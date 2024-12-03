#include "udf.h"
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* dynamic memory allocation for 1D and 2D arrays */
#define DECLARE_MEMORY(name, type) type *name = NULL

#define DECLARE_MEMORY_N(name, type, dim) type *name[dim] = {NULL}

#define RELEASE_MEMORY(name)                                        \
if (NNULLP(name)) {                                                 \
    free(name);                                                     \
    name = NULL;                                                    \
}

#define RELEASE_MEMORY_N(name, dim)                                 \
for (_d = 0; _d < dim; _d++) {                                      \
    RELEASE_MEMORY(name[_d]);                                       \
}

#define ASSIGN_MEMORY(name, size, type)                             \
if (size) {                                                         \
    if (NNULLP(name)) {                                             \
        name = (type *)realloc(name, size * sizeof(type));          \
    } else {                                                        \
        name = (type *)malloc(size * sizeof(type));                 \
    }                                                               \
    if (NULLP(name)) {                                              \
        Error("\nUDF-error: Memory assignment failed for name.");   \
        exit(1);                                                    \
    }                                                               \
}

#define ASSIGN_MEMORY_N(name, size, type, dim)                      \
for (_d = 0; _d < dim; _d++) {                                      \
    ASSIGN_MEMORY(name[_d], size, type);                            \
}

/* send and receive arrays in parallel (see Fluent UDF manual) */
#define PRF_CSEND_INT_N(to, name, n, tag, dim)                      \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CSEND_INT(to, name[_d], n, tag);                            \
}

#define PRF_CSEND_REAL_N(to, name, n, tag, dim)                     \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CSEND_REAL(to, name[_d], n, tag);                           \
}

#define PRF_CRECV_INT_N(from, name, n, tag, dim)                    \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CRECV_INT(from, name[_d], n, tag);                          \
}

#define PRF_CRECV_REAL_N(from, name, n, tag, dim)                   \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CRECV_REAL(from, name[_d], n, tag);                         \
}

/* global variables */
#define mnpf 2
int _d; /* don't use in UDFs! (overwritten by functions above) */
int n_threads;
DECLARE_MEMORY(thread_ids, int);
int timestep = 0;


  /*----------------*/
 /* get_thread_ids */
/*----------------*/

DEFINE_ON_DEMAND(get_thread_ids) {
    /* Read in thread ids, which is the name of the ModelParts in Fluent. Should be called early on. */

#if RP_HOST /* only host process is involved, code not compiled for node */
    int k;
    FILE *file;
    file = fopen("bcs.txt", "r");
    fscanf(file, "%i", &n_threads);
#endif /* RP_HOST */

    host_to_node_int_1(n_threads); /* send n_threads (read from file on host) from host and receive on node(s) */
    ASSIGN_MEMORY(thread_ids, n_threads, int);

#if RP_HOST /* only host process is involved, code not compiled for node */
    for (k = 0; k < n_threads; k++) {
        fscanf(file, "%*s %i", &thread_ids[k]);
    }
    fclose(file);
#endif /* RP_HOST */

    host_to_node_int(thread_ids, n_threads); /* thread_ids (read from file on host) communicated to nodes */
}


  /*----------------------*/
 /* store_coordinates_id */
/*----------------------*/

DEFINE_ON_DEMAND(store_coordinates_id) {
    /* Loop over all ModelParts (in this function denoted as face_thread) and write files nodes_timestep%i_thread%i.dat
    and faces_timestep%i_thread%i.dat. The first contains the coordinates of the nodes and the dynamic mesh node ids,
    which are guaranteed to be unique and consistent. The second contains the face centroid coordinates and the dynamic
    mesh node ids of the nodes (vertices) of that face. The looping over the faces and mesh nodes (vertices) can only be
    done by the compute nodes for their partition of the mesh. Using low-level message passing macros this data is
    collected on node_zero, which communicates it further to the host. The host finally writes the files. */
    if (myid == 0) {printf("\n\nStarted UDF store_coordinates_id.\n"); fflush(stdout);}

    /* declaring variables */
    int thread, n_nodes, n_faces, i_n, i_f, d, compute_node;
    DECLARE_MEMORY_N(node_coords, real, ND_ND);
    DECLARE_MEMORY(node_ids, int);
    DECLARE_MEMORY_N(face_coords, real, ND_ND);
    DECLARE_MEMORY_N(face_ids, int, mnpf);

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    Domain *domain;
    Thread *face_thread; /* points to face thread i.e. number of faces corresponding (to part) of moving interface */
    face_t face; /* variable to keep track of current face in face loop */
    Node *node; /* points to node */
    int node_number, j;
    real centroid[ND_ND]; /* array to store 2 or 3 coordinates */
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
    char file_nodes_name[256];
    char file_faces_name[256];
    FILE *file_nodes = NULL;
    FILE *file_faces = NULL;
    timestep = RP_Get_Integer("udf/timestep"); /* host process reads "udf/timestep" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(timestep); /* host process shares timestep variable with nodes */

    for (thread=0; thread<n_threads; thread++) { /* both host and node execute loop over face_threads (= ModelParts) */

#if RP_HOST /* only host process is involved, code not compiled for node */
        sprintf(file_nodes_name, "nodes_timestep%i_thread%i.dat", timestep, thread_ids[thread]);
        sprintf(file_faces_name, "faces_timestep%i_thread%i.dat", timestep, thread_ids[thread]);

        if (NULLP(file_nodes = fopen(file_nodes_name, "w"))) {
            Error("\nUDF-error: Unable to open %s for writing\n", file_nodes_name);
            exit(1);
        }
        if (NULLP(file_faces = fopen(file_faces_name, "w"))) {
            Error("\nUDF-error: Unable to open %s for writing\n", file_faces_name);
            exit(1);
        }

#if RP_2D
        fprintf(file_nodes, "%27s %27s %10s\n", "x-coordinate", "y-coordinate", "unique-id");
        fprintf(file_faces, "%27s %27s  %10s\n", "x-coordinate", "y-coordinate", "unique-ids");
#else /* RP_2D */
        fprintf(file_nodes, "%27s %27s %27s %10s\n", "x-coordinate", "y-coordinate", "z-coordinate", "unique-id");
        fprintf(file_faces, "%27s %27s %27s  %10s\n", "x-coordinate", "y-coordinate", "z-coordinate", "unique-ids");
#endif /* RP_2D */
#endif /* RP_HOST */

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
        domain = Get_Domain(1);
        face_thread = Lookup_Thread(domain, thread_ids[thread]);

        n_nodes = 0; /* to store number of nodes in current compute node partition of the face_thread */
        /* Loop over all faces (tracked by variable "face") in partition of face_thread to which current compute node
        has access and count number of nodes in each face. Knowing the number of nodes in partition allows assigning
        memory for the coordinate arrays and id arrays */
        begin_f_loop(face, face_thread) {
            n_nodes += F_NNODES(face, face_thread);
        } end_f_loop(face, face_thread)
        n_faces = THREAD_N_ELEMENTS_INT(face_thread); /* get number of faces in this partition of face_thread */

        /* assign memory on compute node that accesses its partition of the face_thread */
        ASSIGN_MEMORY_N(node_coords, n_nodes, real, ND_ND);
        ASSIGN_MEMORY(node_ids, n_nodes, int);
        ASSIGN_MEMORY_N(face_coords, n_faces, real, ND_ND);
        ASSIGN_MEMORY_N(face_ids, n_faces, int, mnpf);

        i_n = 0;
        i_f = 0;
        /* loop over all faces (tracked by variable "face") in current compute node partition of face_thread */
        begin_f_loop(face, face_thread) {
            if (i_f >= n_faces) {Error("\nIndex %i >= array size %i.", i_f, n_faces);}

            F_CENTROID(centroid, face, face_thread); /* store centroid coordinates of face in variable "centroid" */
            for (d = 0; d < ND_ND; d++) {
                face_coords[d][i_f] = centroid[d];
            }

            for (j = 0; j < mnpf; j++) {
                 /* -1 is placeholder, it is usually overwritten, but not always in unstructured grids */
                face_ids[j][i_f] = -1;
            }

            j = 0;
            /* loop over all nodes in the face, "node_number" is the local node index number that keeps track of the
            current node */
            f_node_loop(face, face_thread, node_number) {
                node = F_NODE(face, face_thread, node_number); /* obtain global face node number */

                if (j >= mnpf) {Error("\nIndex %i >= array size %i.", j, mnpf);}
                 /* store dynamic mesh node ids of all nodes of current face, requires
                 (enable-dynamic-mesh-node-ids #t) in journal file */
                face_ids[j][i_f] = NODE_DM_ID(node);
                j++;

                if (i_n >= n_nodes) {Error("\nIndex %i >= array size %i.", i_n, n_nodes);}
                node_ids[i_n] = NODE_DM_ID(node); /* store dynamic mesh node id of current node */
                for (d = 0; d < ND_ND; d++) {
                    node_coords[d][i_n] = NODE_COORD(node)[d];
                }
                i_n++;
            }
            i_f++;
        } end_f_loop(face, face_thread);

        /* assign destination ID compute_node to "node_host" or "node_zero" (these names are known) */
        compute_node = (I_AM_NODE_ZERO_P) ? node_host : node_zero;

        /* send from node to either node_zero, or if the current process is node_zero, to node_host */
        /* the tag argument myid is the ID of the sending node, as the convention is to have the tag argument the same
        as the from argument (that is, the first argument) for receive messages */
        PRF_CSEND_INT(compute_node, &n_nodes, 1, myid);
        PRF_CSEND_INT(compute_node, &n_faces, 1, myid);

        PRF_CSEND_REAL_N(compute_node, node_coords, n_nodes, myid, ND_ND);
        PRF_CSEND_INT(compute_node, node_ids, n_nodes, myid);
        PRF_CSEND_REAL_N(compute_node, face_coords, n_faces, myid, ND_ND);
        PRF_CSEND_INT_N(compute_node, face_ids, n_faces, myid, mnpf);

        /* memory can be freed once the data is sent */
        RELEASE_MEMORY_N(node_coords, ND_ND);
        RELEASE_MEMORY(node_ids);
        RELEASE_MEMORY_N(face_coords, ND_ND);
        RELEASE_MEMORY_N(face_ids, mnpf);

        /* node_zero is the only one that can communicate with host, so it first receives from the other nodes, then
        sends to the host */
        if(I_AM_NODE_ZERO_P){
            compute_node_loop_not_zero(compute_node) { /* loop over all other nodes and receive from each */
                /* the tag argument compute_node is the ID of the sending node, as the convention is to have the tag
                argument the same as the from argument (that is, the first argument) for receive messages */
                PRF_CRECV_INT(compute_node, &n_nodes, 1, compute_node);
                PRF_CRECV_INT(compute_node, &n_faces, 1, compute_node);

                /* Once n_nodes and n_faces has been received, the correct amount of memory can be allocated on compute
                node 0. This depends on the partition assigned to the sending compute node. */
                ASSIGN_MEMORY_N(node_coords, n_nodes, real, ND_ND);
                ASSIGN_MEMORY(node_ids, n_nodes, int);
                ASSIGN_MEMORY_N(face_coords, n_faces, real, ND_ND);
                ASSIGN_MEMORY_N(face_ids, n_faces, int, mnpf);

                /* receive the 2D-arrays from the other nodes on node_zero */
                PRF_CRECV_REAL_N(compute_node, node_coords, n_nodes, compute_node, ND_ND);
                PRF_CRECV_INT(compute_node, node_ids, n_nodes, compute_node);
                PRF_CRECV_REAL_N(compute_node, face_coords, n_faces, compute_node, ND_ND);
                PRF_CRECV_INT_N(compute_node, face_ids, n_faces, compute_node, mnpf);

                /* Send the variables to the host. Deviating from the tag convention, the message tag is now the
                original non-zero compute node on which the mesh data was stored, even though node 0 does the actual
                communication */
                PRF_CSEND_INT(node_host, &n_nodes, 1, compute_node);
                PRF_CSEND_INT(node_host, &n_faces, 1, compute_node);

                PRF_CSEND_REAL_N(node_host, node_coords, n_nodes, compute_node, ND_ND);
                PRF_CSEND_INT(node_host, node_ids, n_nodes, compute_node);
                PRF_CSEND_REAL_N(node_host, face_coords, n_faces, compute_node, ND_ND);
                PRF_CSEND_INT_N(node_host, face_ids, n_faces, compute_node, mnpf);

                /* once all data has been sent to host, memory on the node can be freed */
                RELEASE_MEMORY_N(node_coords, ND_ND);
                RELEASE_MEMORY(node_ids);
                RELEASE_MEMORY_N(face_coords, ND_ND);
                RELEASE_MEMORY_N(face_ids, mnpf);
            }
        }
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
        /* loop over all compute nodes (corresponding to the message tags), receive data and append to file for each */
        compute_node_loop(compute_node) {
            PRF_CRECV_INT(node_zero, &n_nodes, 1, compute_node);
            PRF_CRECV_INT(node_zero, &n_faces, 1, compute_node);

            /* once n_nodes and n_faces has been received, the correct amount of memory can be allocated on the host */
            ASSIGN_MEMORY_N(node_coords, n_nodes, real, ND_ND);
            ASSIGN_MEMORY(node_ids, n_nodes, int);
            ASSIGN_MEMORY_N(face_coords, n_faces, real, ND_ND);
            ASSIGN_MEMORY_N(face_ids, n_faces, int, mnpf);

            /* receive the 2D-arrays from node_zero */
            PRF_CRECV_REAL_N(node_zero, node_coords, n_nodes, compute_node, ND_ND);
            PRF_CRECV_INT(node_zero, node_ids, n_nodes, compute_node);
            PRF_CRECV_REAL_N(node_zero, face_coords, n_faces, compute_node, ND_ND);
            PRF_CRECV_INT_N(node_zero, face_ids, n_faces, compute_node, mnpf);

            for (i_n = 0; i_n < n_nodes; i_n++) {
                for (d = 0; d < ND_ND; d++) {
                    fprintf(file_nodes, "%27.17e ", node_coords[d][i_n]);
                }
                fprintf(file_nodes, "%10d\n", node_ids[i_n]);
            }

            for (i_f = 0; i_f < n_faces; i_f++) {
                for (d = 0; d < ND_ND; d++) {
                    fprintf(file_faces, "%27.17e ", face_coords[d][i_f]);
                }
                for (d = 0; d < mnpf; d++) {
                    fprintf(file_faces, " %10d", face_ids[d][i_f]);
                }
                fprintf(file_faces, "\n");
            }

            /* after files have been appended, the memory can be freed */
            RELEASE_MEMORY_N(node_coords, ND_ND);
            RELEASE_MEMORY(node_ids);
            RELEASE_MEMORY_N(face_coords, ND_ND);
            RELEASE_MEMORY_N(face_ids, mnpf);
        } /* close compute_node_loop */

        fclose(file_nodes);
        fclose(file_faces);
#endif /* RP_HOST */

    } /* close loop over threads */

    if (myid == 0) {printf("\nFinished UDF store_coordinates_id.\n"); fflush(stdout);}
}


  /*-------------------------*/
 /* store_pressure_traction */
/*-------------------------*/


DEFINE_ON_DEMAND(store_pressure_traction) {
    /* Similar as the function store_coordinates_id, but now the pressure and traction acting on the faces are stored in
    a file pressure_traction_timestep%i_thread%i.dat. */
    if (myid == 0) {printf("\nStarted UDF store_pressure_traction.\n"); fflush(stdout);}

    /* declaring variables */
    int thread, n, i, d, compute_node;
    DECLARE_MEMORY_N(array, real, ND_ND + 1);
    DECLARE_MEMORY_N(ids, int, mnpf);

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    Domain *domain;
    Thread *face_thread;
    face_t face;
    Node *node;
    int node_number, j;
    real traction[ND_ND], area[ND_ND]; /* store 2 or 3 traction and up-to-date face area in arrays */
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
    char file_name[256];
    FILE *file = NULL;
    timestep = RP_Get_Integer("udf/timestep"); /* host process reads "udf/timestep" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(timestep); /* host process shares timestep variable with nodes */

    for (thread=0; thread<n_threads; thread++) { /* both host and node execute loop over face_threads (= ModelParts) */

#if RP_HOST /* only host process is involved, code not compiled for node */
        sprintf(file_name, "pressure_traction_timestep%i_thread%i.dat",
                timestep, thread_ids[thread]);

        if (NULLP(file = fopen(file_name, "w"))) {
            Error("\nUDF-error: Unable to open %s for writing\n", file_name);
            exit(1);
        }

#if RP_2D
        fprintf(file, "%27s %27s %27s  %10s\n",
            "x-shear", "y-shear", "pressure", "unique-ids");
#else /* RP_2D */
        fprintf(file, "%27s %27s %27s %27s  %10s\n",
            "x-shear", "y-shear", "z-shear", "pressure", "unique-ids");

#endif /* RP_2D */
#endif /* RP_HOST */

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
        domain = Get_Domain(1);
        face_thread = Lookup_Thread(domain, thread_ids[thread]);

        n = THREAD_N_ELEMENTS_INT(face_thread); /* get number of faces in this partition of face_thread */

        /* assign memory on compute node that accesses its partition of the face_thread */
        ASSIGN_MEMORY_N(array, n, real, ND_ND + 1);
        ASSIGN_MEMORY_N(ids, n, int, mnpf);

        i = 0;
        /* loop over all faces (tracked by variable "face") in current compute node partition of face_thread */
        begin_f_loop(face, face_thread) {
            if (i >= n) {Error("\nIndex %i >= array size %i.", i, n);}

            F_AREA(area, face, face_thread);
            NV_VS(traction, =, F_STORAGE_R_N3V(face, face_thread, SV_WALL_SHEAR), *, -1.0 / NV_MAG(area));
            for (d = 0; d < ND_ND; d++) {
                array[d][i] = traction[d];
            }
            array[ND_ND][i] = F_P(face, face_thread);

            for (j = 0; j < mnpf; j++) {
                /* -1 is placeholder, it is usually overwritten, but not always in unstructured grids */
                ids[j][i] = -1;
            }

            j = 0;
            /* loop over all nodes in the face, "node_number" is the local node index number that keeps track of the
            current node */
            f_node_loop(face, face_thread, node_number) {
                if (j >= mnpf) {Error("\nIndex %i >= array size %i.", j, mnpf);}
                node = F_NODE(face, face_thread, node_number);
                ids[j][i] = NODE_DM_ID(node); /* store dynamic mesh node id of current node */
                j++;
            }
            i++;
        } end_f_loop(face, face_thread);

        /* assign destination ID compute_node to "node_host" or "node_zero" (these names are known) */
        compute_node = (I_AM_NODE_ZERO_P) ? node_host : node_zero;

        /* send from node to either node_zero, or if the current process is node_zero, to node_host */
        /* the tag argument myid is the ID of the sending node, as the convention is to have the tag argument the same
        as the from argument (that is, the first argument) for receive messages */
        PRF_CSEND_INT(compute_node, &n, 1, myid);

        PRF_CSEND_REAL_N(compute_node, array, n, myid, ND_ND + 1);
        PRF_CSEND_INT_N(compute_node, ids, n, myid, mnpf);

        /* memory can be freed once the data is sent */
        RELEASE_MEMORY_N(array, ND_ND + 1);
        RELEASE_MEMORY_N(ids, mnpf);

        /* node_zero is the only one that can communicate with host, so it first receives from the other nodes, then
        sends to the host */
        if(I_AM_NODE_ZERO_P){
            compute_node_loop_not_zero(compute_node) { /* loop over all other nodes and receive from each */
                /* the tag argument compute_node is the ID of the sending node, as the convention is to have the tag
                argument the same as the from argument (that is, the first argument) for receive messages */
                PRF_CRECV_INT(compute_node, &n, 1, compute_node);

                /* Once n has been received, the correct amount of memory can be allocated on compute node 0. This
                depends on the partition assigned to the sending compute node. */
                ASSIGN_MEMORY_N(array, n, real, ND_ND + 1);
                ASSIGN_MEMORY_N(ids, n, int, mnpf);

                /* receive the 2D-arrays from the other nodes on node_zero */
                PRF_CRECV_REAL_N(compute_node, array, n, compute_node, ND_ND + 1);
                PRF_CRECV_INT_N(compute_node, ids, n, compute_node, mnpf);

                /* Send the variables to the host. Deviating from the tag convention, the message tag is now the
                original non-zero compute node on which the mesh data was stored, even though node 0 does the actual
                communication */
                PRF_CSEND_INT(node_host, &n, 1, compute_node);

                PRF_CSEND_REAL_N(node_host, array, n, compute_node, ND_ND + 1);
                PRF_CSEND_INT_N(node_host, ids, n, compute_node, mnpf);

                /* once all data has been sent to host, memory on the node can be freed */
                RELEASE_MEMORY_N(array, ND_ND + 1);
                RELEASE_MEMORY_N(ids, mnpf);
            }
        }
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
        /* loop over all compute nodes (corresponding to the message tags), receive data and append to file for each */
        compute_node_loop(compute_node) {
            PRF_CRECV_INT(node_zero, &n, 1, compute_node);

            /* once n has been received, the correct amount of memory can be allocated on the host */
            ASSIGN_MEMORY_N(array, n, real, ND_ND + 1);
            ASSIGN_MEMORY_N(ids, n, int, mnpf);

            /* receive the 2D-arrays from node_zero */
            PRF_CRECV_REAL_N(node_zero, array, n, compute_node, ND_ND + 1);
            PRF_CRECV_INT_N(node_zero, ids, n, compute_node, mnpf);

            for (i = 0; i < n; i++) {
                for (d = 0; d < ND_ND + 1; d++) {
                    fprintf(file, "%27.17e ", array[d][i]);
                }
                for (d = 0; d < mnpf; d++) {
                    fprintf(file, " %10d", ids[d][i]);
                }
                fprintf(file, "\n");
            }
            /* after files have been appended, the memory can be freed */
            RELEASE_MEMORY_N(array, ND_ND + 1);
            RELEASE_MEMORY_N(ids, mnpf);
        } /* close compute_node_loop */

        fclose(file);
#endif /* RP_HOST */

    } /* close loop over threads */

    if (myid == 0) {printf("\nFinished UDF store_pressure_traction.\n"); fflush(stdout);}
}


  /*------------*/
 /* move_nodes */
/*------------*/

DEFINE_GRID_MOTION(move_nodes, domain, dynamic_thread, time, dtime) {
    /* UDF that can be assigned to a deformable wall in Fluent. It will read the updated positions of the mesh nodes
    (vertices) and apply them. As only the compute nodes have access to their own partition of the mesh, each one is
    responsible for reading the file containing the updated coordinates and selecting the correct row. */
    char file_name[256];
    Thread *face_thread = DT_THREAD(dynamic_thread); /* face_thread to which UDF is assigned in Fluent */
    int thread_id = THREAD_ID(face_thread);

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    face_t face;
    Node *node;
    int i, d, n, node_number;
    DECLARE_MEMORY_N(coords, real, ND_ND);
    DECLARE_MEMORY(ids, int);
    FILE *file = NULL;
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
    timestep = RP_Get_Integer("udf/timestep"); /* host process reads "udf/timestep" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(timestep); /* host process shares timestep variable with nodes */

#if RP_HOST /* only host process is involved, code not compiled for node */
    sprintf(file_name, "nodes_update_timestep%i_thread%i.dat",
            timestep, thread_id);
    host_to_node_sync_file(file_name); /* send file to the compute nodes */
#else
    struct stat st = {0};
    /* create a temporary directory if it does not exist yet, this is needed as multiple machines can be involved */
    if (stat("/tmp/coconut_vsc48608_2211195_fluent", &st) == -1) {
        mkdir("/tmp/coconut_vsc48608_2211195_fluent", 0700);
    }
    sprintf(file_name, "/tmp/coconut_vsc48608_2211195_fluent/nodes_update_timestep%i_thread%i.dat",
            timestep, thread_id);
    host_to_node_sync_file("/tmp/coconut_vsc48608_2211195_fluent");  /* receive file on compute nodes and store in a temporary folder */
#endif /* RP_HOST */

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    if (NULLP(file = fopen(file_name, "r"))) {
        Error("\nUDF-error: Unable to open %s for reading\n", file_name);
        exit(1);
    }

    fscanf(file, "%i", &n); /* first line of file says how many nodes (=lines) there are */

    ASSIGN_MEMORY_N(coords, n, real, ND_ND);
    ASSIGN_MEMORY(ids, n, int);

    for (i=0; i < n; i++) {
        for (d = 0; d < ND_ND; d++) {
            fscanf(file, "%lf", &coords[d][i]); /* read coordinates from file */
        }
        fscanf(file, "%i", &ids[i]); /* read node ids from file */
    }

    fclose(file);

    begin_f_loop(face, face_thread) { /* loop over all faces in face_thread */
        f_node_loop(face, face_thread, node_number) { /* loop over all nodes in current face */
            node = F_NODE(face, face_thread, node_number); /* get global face ndoe index from local node index */
            if NODE_POS_NEED_UPDATE(node) { /* only execute if position has not been updated yet (efficiency) */
                for (i=0; i < n; i++) { /* loop over all lines to find the correct node */
                    if (NODE_DM_ID(node) == ids[i]) { /* correct node has the same dynamic mesh node id */
                        for (d = 0; d < ND_ND; d++) {
                            NODE_COORD(node)[d] = coords[d][i]; /* modify node coordinates */
                        }
                        NODE_POS_UPDATED(node); /* flag node as updated */
                        break; /* break loop for efficiency */
                    }
                    if (i == n - 1) {
                        Error("\nUDF-error: No match for node id %i\n", NODE_DM_ID(node));
                        exit(1);
                    }
                }
            }
        }
    } end_f_loop(face, face_thread);

    RELEASE_MEMORY_N(coords, ND_ND);
    RELEASE_MEMORY(ids);

    if (myid == 0) {
        sprintf(file_name, "/tmp/coconut_vsc48608_2211195_fluent/nodes_update_timestep%i_thread%i.dat",
                timestep-1, thread_id);
        remove(file_name);}
#endif /* RP_NODE */

    if (myid == 0) {printf("\nFinished UDF move_nodes.\n"); fflush(stdout);}
}
