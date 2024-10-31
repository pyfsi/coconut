#include "udf.h"
#include "sg.h"
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

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

/* User defined memory -- not available in 3D! */
#define PR_1_X 0 // previous x-coordinate of first node in face
#define PR_1_Y 1 // previous y-coordinate of first node in face
#define PR_2_X 2 // previous x-coordinate of second node in face
#define PR_2_Y 3 // previous y-coordinate of second node in face
#define NEW_1_X 4 // new x-coordinate of first node in face
#define NEW_1_Y 5 // new y-coordinate of first node in face
#define NEW_2_X 6 // new x-coordinate of second node in face
#define NEW_2_Y 7 // new y-coordinate of second node in face
#define ADJ 8 // flag for cells adjacent to the interface
#define D_VOL 9 // swept volume during move_nodes operation
#define PR_H 10 // Cell liquid fraction in previous converged timestep
#define SIGN 11 // defines wether the cell is on the growing or shrinking side during phase change

/* global variables */
#define mnpf |MAX_NODES_PER_FACE|
real dt = |TIME_STEP_SIZE|;
real LH = |LATENT_HEAT|;
real TM = |MELT_TEMP|;
bool unsteady = |UNSTEADY|;
bool fluid = |FLUID|;
int _d; /* don't use in UDFs! (overwritten by functions above) */
int n_threads;
DECLARE_MEMORY(thread_ids, int);
int timestep = 0;
int iteration = 0;

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


  /*-------------------*/
 /* store_temperature */
/*-------------------*/

DEFINE_ON_DEMAND(store_temperature){
    /* Similar as the function store_coordinates_id and store_pressure_traction, but now the temperature at the faces
    is stored in a file temperature_timestep%i_thread%i.dat. */
    if (myid == 0) {printf("\nStarted UDF store_temperature.\n"); fflush(stdout);}

    /* declaring variables */
    int thread, n, i, d, compute_node;
    DECLARE_MEMORY(temp, real);
    DECLARE_MEMORY_N(ids, int, mnpf);

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    Domain *domain;
    Thread *face_thread;
    face_t face;
    Node *node;
    int node_number, j;
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
    char file_name[256];
    FILE *file = NULL;
    timestep = RP_Get_Integer("udf/timestep"); /* host process reads "udf/timestep" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(timestep); /* host process shares timestep variable with nodes */

    for (thread=0; thread<n_threads; thread++) { /* both host and node execute loop over face_threads (= ModelParts) */

#if RP_HOST /* only host process is involved, code not compiled for node */
        sprintf(file_name, "temperature_timestep%i_thread%i.dat",
                timestep, thread_ids[thread]);

        if (NULLP(file = fopen(file_name, "w"))) {
            Error("\nUDF-error: Unable to open %s for writing\n", file_name);
            exit(1);
        }

        fprintf(file, "%27s %10s\n",
            "temperature", "unique-ids");

#endif /* RP_HOST */

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
        domain = Get_Domain(1);
        face_thread = Lookup_Thread(domain, thread_ids[thread]);

        n = THREAD_N_ELEMENTS_INT(face_thread); /* get number of faces in this partition of face_thread */

        /* assign memory on compute node that accesses its partition of the face_thread */
        ASSIGN_MEMORY(temp, n, real);
        ASSIGN_MEMORY_N(ids, n, int, mnpf);

        i = 0;
        /* loop over all faces (tracked by variable "face") in current compute node partition of face_thread */
        begin_f_loop(face, face_thread) {
            if (i >= n) {Error("\nIndex %i >= array size %i.", i, n);}

            /* Strategy to allow subcooling of the solid domain */
            if (fluid) {
                temp[i] = F_T(face, face_thread);
            } else {
                if (F_T(face, face_thread) >= TM) {
                    temp[i] = TM;
                } else {
                    temp[i] = F_T(face, face_thread);
                }
            }

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

        PRF_CSEND_REAL(compute_node, temp, n, myid);
        PRF_CSEND_INT_N(compute_node, ids, n, myid, mnpf);

        /* memory can be freed once the data is sent */
        RELEASE_MEMORY(temp);
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
                ASSIGN_MEMORY(temp, n, real);
                ASSIGN_MEMORY_N(ids, n, int, mnpf);

                /* receive the 2D-arrays from the other nodes on node_zero */
                PRF_CRECV_REAL(compute_node, temp, n, compute_node);
                PRF_CRECV_INT_N(compute_node, ids, n, compute_node, mnpf);

                /* Send the variables to the host. Deviating from the tag convention, the message tag is now the
                original non-zero compute node on which the mesh data was stored, even though node 0 does the actual
                communication */
                PRF_CSEND_INT(node_host, &n, 1, compute_node);

                PRF_CSEND_REAL(node_host, temp, n, compute_node);
                PRF_CSEND_INT_N(node_host, ids, n, compute_node, mnpf);

                /* once all data has been sent to host, memory on the node can be freed */
                RELEASE_MEMORY(temp);
                RELEASE_MEMORY_N(ids, mnpf);
            }
        }
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
        /* loop over all compute nodes (corresponding to the message tags), receive data and append to file for each */
        compute_node_loop(compute_node) {
            PRF_CRECV_INT(node_zero, &n, 1, compute_node);

            /* once n has been received, the correct amount of memory can be allocated on the host */
            ASSIGN_MEMORY(temp, n, real);
            ASSIGN_MEMORY_N(ids, n, int, mnpf);

            /* receive the 2D-arrays from node_zero */
            PRF_CRECV_REAL(node_zero,temp, n, compute_node);
            PRF_CRECV_INT_N(node_zero, ids, n, compute_node, mnpf);

            for (i = 0; i < n; i++) {

                fprintf(file, "%27.17e ", temp[i]);

                for (d = 0; d < mnpf; d++) {
                    fprintf(file, " %10d", ids[d][i]);
                }
                fprintf(file, "\n");
            }
            /* after files have been appended, the memory can be freed */
            RELEASE_MEMORY(temp);
            RELEASE_MEMORY_N(ids, mnpf);
        } /* close compute_node_loop */

        fclose(file);
#endif /* RP_HOST */

    } /* close loop over threads */

    if (myid == 0) {printf("\nFinished UDF store_temperature.\n"); fflush(stdout);}
}


  /*-----------------*/
 /* store_heat_flux */
/*-----------------*/

DEFINE_ON_DEMAND(store_heat_flux){
    /* Similar as the function store_temperature, but now the heat flux magnitude projected on the normal
    direction outwards of the faces is stored in a file heat_flux_timestep%i_thread%i.dat. */
    if (myid == 0) {printf("\nStarted UDF store_heat_flux.\n"); fflush(stdout);}

    /* declaring variables */
    int thread, n, i, d, compute_node;
    DECLARE_MEMORY(flux, real); /* heat flux magnitude outwards of the faces */
    DECLARE_MEMORY_N(ids, int, mnpf);

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    Domain *domain;
    Thread *face_thread;
    face_t face;
    Node *node;
    int node_number, j;
    real area[ND_ND]; /* store 2 or 3D face area vector in arrays */
    real heat;
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
    char file_name[256];
    FILE *file = NULL;
    timestep = RP_Get_Integer("udf/timestep"); /* host process reads "udf/timestep" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(timestep); /* host process shares timestep variable with nodes */

    for (thread=0; thread<n_threads; thread++) { /* both host and node execute loop over face_threads (= ModelParts) */

#if RP_HOST /* only host process is involved, code not compiled for node */
        sprintf(file_name, "heat_flux_timestep%i_thread%i.dat",
                timestep, thread_ids[thread]);

        if (NULLP(file = fopen(file_name, "w"))) {
            Error("\nUDF-error: Unable to open %s for writing\n", file_name);
            exit(1);
        }

        fprintf(file, "%27s %10s\n",
            "flux-normal", "unique-ids");

#endif /* RP_HOST */

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
        domain = Get_Domain(1);
        face_thread = Lookup_Thread(domain, thread_ids[thread]);

        n = THREAD_N_ELEMENTS_INT(face_thread); /* get number of faces in this partition of face_thread */

        /* assign memory on compute node that accesses its partition of the face_thread */
        ASSIGN_MEMORY(flux, n, real);
        ASSIGN_MEMORY_N(ids, n, int, mnpf);

        i = 0;
        /* loop over all faces (tracked by variable "face") in current compute node partition of face_thread */
        begin_f_loop(face, face_thread) {
            if (i >= n) {Error("\nIndex %i >= array size %i.", i, n);}

            F_AREA(area, face, face_thread);
            heat = BOUNDARY_HEAT_FLUX(face, face_thread);
            flux[i] = heat/NV_MAG(area);

            /* Code currently only for melting: store_heat_flux udf only used in liquid domain & liquid domain will grow */
            C_UDMI(F_C0(face, face_thread),THREAD_T0(face_thread),SIGN) = 1.0;

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

        PRF_CSEND_REAL(compute_node, flux, n, myid);
        PRF_CSEND_INT_N(compute_node, ids, n, myid, mnpf);

        /* memory can be freed once the data is sent */
        RELEASE_MEMORY(flux);
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
                ASSIGN_MEMORY(flux, n, real);
                ASSIGN_MEMORY_N(ids, n, int, mnpf);

                /* receive the 2D-arrays from the other nodes on node_zero */
                PRF_CRECV_REAL(compute_node, flux, n, compute_node);
                PRF_CRECV_INT_N(compute_node, ids, n, compute_node, mnpf);

                /* Send the variables to the host. Deviating from the tag convention, the message tag is now the
                original non-zero compute node on which the mesh data was stored, even though node 0 does the actual
                communication */
                PRF_CSEND_INT(node_host, &n, 1, compute_node);

                PRF_CSEND_REAL(node_host, flux, n, compute_node);
                PRF_CSEND_INT_N(node_host, ids, n, compute_node, mnpf);

                /* once all data has been sent to host, memory on the node can be freed */
                RELEASE_MEMORY(flux);
                RELEASE_MEMORY_N(ids, mnpf);
            }
        }
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
        /* loop over all compute nodes (corresponding to the message tags), receive data and append to file for each */
        compute_node_loop(compute_node) {
            PRF_CRECV_INT(node_zero, &n, 1, compute_node);

            /* once n has been received, the correct amount of memory can be allocated on the host */
            ASSIGN_MEMORY(flux, n, real);
            ASSIGN_MEMORY_N(ids, n, int, mnpf);

            /* receive the 2D-arrays from node_zero */
            PRF_CRECV_REAL(node_zero, flux, n, compute_node);
            PRF_CRECV_INT_N(node_zero, ids, n, compute_node, mnpf);

            for (i = 0; i < n; i++) {

                fprintf(file, "%27.17e ", flux[i]);

                for (d = 0; d < mnpf; d++) {
                    fprintf(file, " %10d", ids[d][i]);
                }
                fprintf(file, "\n");
            }
            /* after files have been appended, the memory can be freed */
            RELEASE_MEMORY(flux);
            RELEASE_MEMORY_N(ids, mnpf);
        } /* close compute_node_loop */

        fclose(file);
#endif /* RP_HOST */

    } /* close loop over threads */

    if (myid == 0) {printf("\nFinished UDF store_heat_flux.\n"); fflush(stdout);}
}


  /*-----------------*/
 /* set_temperature */
/*-----------------*/

DEFINE_PROFILE(set_temperature, face_thread, var) {
    /* UDF that can be used to define a custom temperature boundary profile in Fluent. It will read the updated
    temperature profile from a file and set them as thermal boundary condition. As only the compute nodes have access
    to their own partition of the mesh, each one is responsible for reading the file containing the updated
    temperatures and selecting the correct row. */
    if (myid == 0) {printf("\nStarted UDF set_temperature.\n"); fflush(stdout);}
    char file_name[256];
    int thread_id = THREAD_ID(face_thread);
    int timestep;
    bool skip_search; // Boolean to indicate whether loop should search for corresponding node id's
    skip_search = false;

// DEFINE_PROFILE UDF does not execute HOST processes -> only node associated with given face thread

#if RP_NODE
    face_t face;
    Node *node;
    int i, d, n, node_number, id;
    DECLARE_MEMORY(temp, real);
    DECLARE_MEMORY(flag, bool);
    DECLARE_MEMORY_N(ids, int, mnpf);
    FILE *file = NULL;
#endif /* RP_NODE */

    if (unsteady) {
        timestep = N_TIME;
    } else {
        timestep = N_ITER;
        if (timestep != 0) {
            timestep = 1;
        }
    }

    sprintf(file_name, "temperature_timestep%i_thread%i.dat",
            timestep, thread_id);

#if RP_NODE
    if (NULLP(file = fopen(file_name, "r"))) {
        if (timestep == 0) {
            skip_search = true;
        } else {
            Error("\nUDF-error: Unable to open %s for reading\n", file_name);
            exit(1);
        }
    }

    if (!skip_search) {
        char line[256];
        n = 0;
        while (fgets(line, sizeof(line), file) != NULL) {
                n++;
            }
        n--;
        char line_bis[64];
        fseek(file, 0L, SEEK_SET);
        fgets(line_bis, sizeof(line_bis), file); // Discard the header line

        ASSIGN_MEMORY(temp, n, real);
        ASSIGN_MEMORY(flag, n, bool);
        ASSIGN_MEMORY_N(ids, n, int, mnpf);

        for (i=0; i < n; i++) {
            fscanf(file, "%lf", &temp[i]); /* read temperature from file */
            flag[i] = false;
            for (d = 0; d < mnpf; d++) {
                fscanf(file, "%i", &ids[d][i]); /* read node ids from file */
            }
        }

        fclose(file);
    }

    char error_msg[256] = "UDF-error: No match for face with node ids";
    bool isNodeFound;
    begin_f_loop(face, face_thread) { /* loop over all faces in face_thread */
        if (!skip_search) {
            for (i=0; i < n; i++) { /* loop over all faces in ids array */
                if (flag[i]) { /* skip faces in ids array that have been assigned already */
                    continue;
                } else {
                    isNodeFound = true;
                    f_node_loop(face, face_thread, node_number) { /* loop over all nodes in current face */
                        if (!isNodeFound) { /* break loop if previous node was not found */
                            break;
                        }
                        node = F_NODE(face, face_thread, node_number); /* get global face node index from local node index */
                        if (NNULLP(THREAD_STORAGE(NODE_THREAD(node), SV_DM_ID))) { // checks if node values are already loaded
                            id = NODE_DM_ID(node);
                        } else {
                            skip_search = true;
                            isNodeFound = false;
                            break;
                        }
                        for (d = 0; d < mnpf; d++) { /* loop over all node ids for current face in ids array */
                            if (id == ids[d][i]) {
                                isNodeFound = true;
                                break;
                            } else {
                                isNodeFound = false;
                            }
                        }
                    }
                    if (skip_search) {
                        break;
                    }
                    if (isNodeFound) { /* All nodes have been found, so the faces match */
                        flag[i] = true;
                        F_PROFILE(face, face_thread, var) = temp[i];
                        break;
                    }
                }
            }
        }
        if (skip_search) {
            F_PROFILE(face, face_thread, var) = 300.0; // assign arbitrary value when faces cannot be identified by node id's
        }
        if (!isNodeFound && !skip_search) {
            for (d = 0; d < mnpf; d++) {
                char nodeID[11];
                sprintf(nodeID, " %d", ids[d][i]);
                strcat(error_msg, nodeID);
                if (d < mnpf - 1) {
                    strcat(error_msg, ",");
                }
            }
            Error("\n%s\n", error_msg);
            exit(1);
        }
    } end_f_loop(face, face_thread);

    RELEASE_MEMORY(temp);
    RELEASE_MEMORY(flag);
    RELEASE_MEMORY_N(ids, mnpf);
#endif /* RP_NODE */

    if (myid == 0) {printf("\nFinished UDF set_temperature.\n"); fflush(stdout);}
}


  /*---------------*/
 /* set_heat_flux */
/*---------------*/

DEFINE_PROFILE(set_heat_flux, face_thread, var) {
    /* UDF that can be used to define a custom heat flux boundary profile in Fluent, similar to set_temperature.
    It will read the updated heat flux profile from a file and set them as thermal boundary condition.
    As only the compute nodes have access to their own partition of the mesh, each one is responsible for reading
    the file containing the updated heat fluxes and selecting the correct row. */
    if (myid == 0) {printf("\nStarted UDF set_heat_flux.\n"); fflush(stdout);}
    char file_name[256];
    int thread_id = THREAD_ID(face_thread);
    int timestep;
    bool skip_search; // Boolean to indicate whether loop should search for corresponding node id's
    skip_search = false;

#if RP_NODE
    face_t face;
    Node *node;
    int i, d, n, node_number, id, id_bis;
    DECLARE_MEMORY(heat_flux, real);
    DECLARE_MEMORY(flag, bool);
    DECLARE_MEMORY_N(ids, int, mnpf);
    FILE *file = NULL;
#endif /* RP_NODE */

    if (unsteady) {
        timestep = N_TIME;
    } else {
        timestep = N_ITER;
        if (timestep != 0) {
            timestep = 1;
        }
    }

    sprintf(file_name, "heat_flux_timestep%i_thread%i.dat",
            timestep, thread_id);

#if RP_NODE
    if (NULLP(file = fopen(file_name, "r"))) {
        if (timestep == 0) {
            skip_search = true;
        } else {
            Error("\nUDF-error: Unable to open %s for reading\n", file_name);
            exit(1);
        }
    }

    if (!skip_search) {
        char line[256];
        n = 0;
        while (fgets(line, sizeof(line), file) != NULL) {
                n++;
            }
        n--;
        char line_bis[64];
        fseek(file, 0L, SEEK_SET);
        fgets(line_bis, sizeof(line_bis), file); // Discard the header line

        ASSIGN_MEMORY(heat_flux, n, real);
        ASSIGN_MEMORY(flag, n, bool);
        ASSIGN_MEMORY_N(ids, n, int, mnpf);

        for (i = 0; i < n; i++) {
            fscanf(file, "%lf", &heat_flux[i]); /* read normal incoming heat flux from file */
            flag[i] = false;
            for (d = 0; d < mnpf; d++) {
                fscanf(file, "%i", &ids[d][i]); /* read node ids from file */
            }
        }

        fclose(file);
    }

    char error_msg[256] = "UDF-error: No match for face with node ids";
    bool isNodeFound;
    begin_f_loop(face, face_thread) { /* loop over all faces in face_thread */
        if (!skip_search) {
            for (i=0; i < n; i++) { /* loop over all faces in ids array */
                if (flag[i]) { /* skip faces in ids array that have been assigned already */
                    continue;
                } else {
                    isNodeFound = true;
                    f_node_loop(face, face_thread, node_number) { /* loop over all nodes in current face */
                        if (!isNodeFound) { /* break loop if previous node was not found */
                            break;
                        }
                        node = F_NODE(face, face_thread, node_number); /* get global face node index from local node index */
                        if (NNULLP(THREAD_STORAGE(NODE_THREAD(node), SV_DM_ID))) { // checks if node values are already loaded
                            id = NODE_DM_ID(node);
                        } else {
                            skip_search = true;
                            isNodeFound = false;
                            break;
                        }
                        for (d = 0; d < mnpf; d++) { /* loop over all node ids for current face in ids array */
                            if (id == ids[d][i]) {
                                isNodeFound = true;
                                break;
                            } else {
                                isNodeFound = false;
                            }
                        }
                    }
                    if (skip_search) {
                        break;
                    }
                    if (isNodeFound) { /* All nodes have been found, so the faces match */
                        flag[i] = true;
                        F_PROFILE(face, face_thread, var) = -1*heat_flux[i];
                        /* Code currently only for melting: set_heat_flux udf only used in solid domain & solid domain will shrink */
                        C_UDMI(F_C0(face, face_thread),THREAD_T0(face_thread),SIGN) = -1.0;
                        break;
                    }
                }
            }
        }
        if (skip_search) {
            F_PROFILE(face, face_thread, var) = 0.0; // assign zero value when faces cannot be identified by node id's
        }
        if (!isNodeFound && !skip_search) {
            for (d = 0; d < mnpf; d++) {
                char nodeID[11];
                sprintf(nodeID, " %d", ids[d][i]);
                strcat(error_msg, nodeID);
                if (d < mnpf - 1) {
                    strcat(error_msg, ",");
                }
            }
            Error("\n%s\n", error_msg);
            exit(1);
        }
    } end_f_loop(face, face_thread);

    RELEASE_MEMORY(heat_flux);
    RELEASE_MEMORY(flag);
    RELEASE_MEMORY_N(ids, mnpf);
#endif /* RP_NODE */

    if (myid == 0) {printf("\nFinished UDF set_heat_flux.\n"); fflush(stdout);}
}


  /*------------*/
 /* move_nodes */
/*------------*/

DEFINE_GRID_MOTION(move_nodes, domain, dynamic_thread, time, dtime) {
    /* UDF that can be assigned to a deformable wall in Fluent. It will read the updated positions of the mesh nodes
    (vertices) and apply them. As only the compute nodes have access to their own partition of the mesh, each one is
    responsible for reading the file containing the updated coordinates and selecting the correct row. */
    if (myid == 0) {printf("\nStarted UDF move_nodes.\n"); fflush(stdout);}
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
    if (stat("|TMP_DIRECTORY_NAME|", &st) == -1) {
        mkdir("|TMP_DIRECTORY_NAME|", 0700);
    }
    sprintf(file_name, "|TMP_DIRECTORY_NAME|/nodes_update_timestep%i_thread%i.dat",
            timestep, thread_id);
    host_to_node_sync_file("|TMP_DIRECTORY_NAME|");  /* receive file on compute nodes and store in a temporary folder */
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
        sprintf(file_name, "|TMP_DIRECTORY_NAME|/nodes_update_timestep%i_thread%i.dat",
                timestep-1, thread_id);
        remove(file_name);}
#endif /* RP_NODE */

    if (myid == 0) {printf("\nFinished UDF move_nodes.\n"); fflush(stdout);}
}

  /*--------------*/
 /* set_adjacent */
/*--------------*/

DEFINE_ON_DEMAND(set_adjacent) {
/* Store previous and new interface node coordinates in UDM's and flags cells adjacent to the boundary.*/

#if RP_HOST /* only host process is involved, code not compiled for node */
    iteration = RP_Get_Integer("udf/iteration"); /* host process reads "udf/iteration" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(iteration); /* host process shares iteration variable with nodes */

#if RP_NODE
    Domain *domain;
    Thread *cell_thread, *face_thread;
    cell_t cell;
    face_t face;
    Node *node;
    int face_number, node_number, surface, i;

    domain = Get_Domain(1);
    thread_loop_c(cell_thread,domain)
    {
        begin_c_loop(cell,cell_thread)
        {
            /* Initialize all UDMI's at zero value */
            C_UDMI(cell,cell_thread,ADJ) = 0.0;
            c_face_loop(cell,cell_thread,face_number)
            {
                face_thread = C_FACE_THREAD(cell,cell_thread,face_number);
                for (surface = 0; surface < n_threads; surface++)
                {
                    if (THREAD_ID(face_thread) == thread_ids[surface])
                    {
                        C_UDMI(cell,cell_thread,ADJ) = 1.0;
                        face = C_FACE(cell,cell_thread,face_number);
                        i = 0;
                        f_node_loop(face, face_thread, node_number) {
                            node = F_NODE(face, face_thread, node_number);
                            if (iteration == 0) { // only executed once at start of journal file
                                if (i==0) {
                                    C_UDMI(cell,cell_thread,NEW_1_X) = NODE_X(node);
                                    C_UDMI(cell,cell_thread,NEW_1_Y) = NODE_Y(node);
                                    C_UDMI(cell,cell_thread,PR_1_X) = NODE_X(node);
                                    C_UDMI(cell,cell_thread,PR_1_Y) = NODE_Y(node);
                                    i ++;
                                }
                                else {
                                    if (i==1) {
                                        C_UDMI(cell,cell_thread,NEW_2_X) = NODE_X(node);
                                        C_UDMI(cell,cell_thread,NEW_2_Y) = NODE_Y(node);
                                        C_UDMI(cell,cell_thread,PR_2_X) = NODE_X(node);
                                        C_UDMI(cell,cell_thread,PR_2_Y) = NODE_Y(node);
                                        i --;
                                    }
                                }
                            }
                            else { // iteration != 0
                                if (iteration == 1) {
                                    if (i==0) {
                                        C_UDMI(cell,cell_thread,PR_1_X) = C_UDMI(cell,cell_thread,NEW_1_X);
                                        C_UDMI(cell,cell_thread,PR_1_Y) = C_UDMI(cell,cell_thread,NEW_1_Y);
                                        C_UDMI(cell,cell_thread,NEW_1_X) = NODE_X(node);
                                        C_UDMI(cell,cell_thread,NEW_1_Y) = NODE_Y(node);
                                        i ++;
                                    }
                                    else {
                                        if (i==1) {
                                            C_UDMI(cell,cell_thread,PR_2_X) = C_UDMI(cell,cell_thread,NEW_2_X);
                                            C_UDMI(cell,cell_thread,PR_2_Y) = C_UDMI(cell,cell_thread,NEW_2_Y);
                                            C_UDMI(cell,cell_thread,NEW_2_X) = NODE_X(node);
                                            C_UDMI(cell,cell_thread,NEW_2_Y) = NODE_Y(node);
                                            i --;
                                        }
                                    }
                                }
                                else { // iteration != 0 & 1
                                    if (i==0) {
                                        C_UDMI(cell,cell_thread,NEW_1_X) = NODE_X(node);
                                        C_UDMI(cell,cell_thread,NEW_1_Y) = NODE_Y(node);
                                        i ++;
                                    }
                                    else {
                                        if (i==1) {
                                            C_UDMI(cell,cell_thread,NEW_2_X) = NODE_X(node);
                                            C_UDMI(cell,cell_thread,NEW_2_Y) = NODE_Y(node);
                                            i --;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } end_c_loop(c,cell_thread)
    }
#endif /* RP_NODE */
}


  /*--------------------*/
 /* calc_volume_change */
/*--------------------*/

DEFINE_ON_DEMAND(calc_volume_change)
{
/* WARNING: only 2D compatibility, no 3D */
    bool zero_vol; // Boolean to indicate whether swept volume is zero or not
    zero_vol = false;

#if RP_NODE
    Domain *domain;
    domain = Get_Domain(1);
    Thread *cell_thread;
    cell_t cell;
    real x_ref, y_ref, dummy;
    int i, j, d, dum_i;
    DECLARE_MEMORY(theta, real); // will be sorted
    DECLARE_MEMORY_N(vertices, real, ND_ND); // NOT sorted
    DECLARE_MEMORY_N(sort_vert, real, ND_ND); // sorted
    DECLARE_MEMORY(index, int); // will be sorted

    ASSIGN_MEMORY(theta, 4, real);
    ASSIGN_MEMORY_N(vertices, 4, real, ND_ND);
    ASSIGN_MEMORY_N(sort_vert, 4, real, ND_ND);
    ASSIGN_MEMORY(index, 4, int);

    thread_loop_c(cell_thread,domain)
    {
        begin_c_loop(cell,cell_thread) //loop over all cells
        {
            zero_vol = false;
            if (C_UDMI(cell,cell_thread,ADJ) == 1.0) {
                vertices[0][0] = C_UDMI(cell,cell_thread,PR_1_X);
                vertices[0][1] = C_UDMI(cell,cell_thread,PR_2_X);
                vertices[0][2] = C_UDMI(cell,cell_thread,NEW_1_X);
                vertices[0][3] = C_UDMI(cell,cell_thread,NEW_2_X);
                vertices[1][0] = C_UDMI(cell,cell_thread,PR_1_Y);
                vertices[1][1] = C_UDMI(cell,cell_thread,PR_2_Y);
                vertices[1][2] = C_UDMI(cell,cell_thread,NEW_1_Y);
                vertices[1][3] = C_UDMI(cell,cell_thread,NEW_2_Y);

                // No volume change when nodes don't move
                if ((vertices[0][0] == vertices[0][2]) && (vertices[1][0] == vertices[1][2])) {
                    if ((vertices[0][1] == vertices[0][3]) && (vertices[1][1] == vertices[1][3])) {
                        zero_vol = true; // PR_1 == NEW_1 & PR_2 == NEW_2
                    }
                }
                if ((vertices[0][0] == vertices[0][3]) && (vertices[1][0] == vertices[1][3])) {
                    if ((vertices[0][1] == vertices[0][2]) && (vertices[1][1] == vertices[1][2])) {
                        zero_vol = true; // PR_1 == NEW_2 & PR_2 == NEW_1
                    }
                }

                if (zero_vol) {
                    C_UDMI(cell,cell_thread,D_VOL) = 0.0;
                }
                else {
                    // determine reference coordinates
                    x_ref = (vertices[0][0]+vertices[0][1]+vertices[0][2]+vertices[0][3])/4;
                    y_ref = (vertices[1][0]+vertices[1][1]+vertices[1][2]+vertices[1][3])/4;
                    for (i=0; i < 4; i++) { // determine angles wrt reference point
                        theta[i] = atan2((vertices[1][i]-y_ref),(vertices[0][i])-x_ref);
                        index[i] = i;
                    }
                    for (i=0; i < 3; i++) { // bubble sort indices in counter-clockwise manner (necessary for shoelace theorem)
                        for (j=0; j < 3-i; j++) {
                            if (theta[j] > theta[j+1]) {
                                dummy = theta[j];
                                theta[j] = theta[j+1];
                                theta[j+1] = dummy;
                                dum_i = index[j];
                                index[j] = index[j+1];
                                index[j+1] = dum_i;
                            }
                        }
                    }
                    for (i=0; i < 4; i++) { // sort coordinates according sorted indices
                        for (d = 0; d < ND_ND; d++) {
                            sort_vert[d][i] = vertices[d][index[i]];
                        }
                    } // calculate swept area with shoelace theorem
                    C_UDMI(cell,cell_thread,D_VOL) = 0.5*((sort_vert[0][0]*sort_vert[1][1]+sort_vert[0][1]*sort_vert[1][2]+sort_vert[0][2]*sort_vert[1][3]+sort_vert[0][3]*sort_vert[1][0])-(sort_vert[0][1]*sort_vert[1][0]+sort_vert[0][2]*sort_vert[1][1]+sort_vert[0][3]*sort_vert[1][2]+sort_vert[0][0]*sort_vert[1][3]));
                }
            }
        } end_c_loop(cell,cell_thread)
    }

    RELEASE_MEMORY(theta);
    RELEASE_MEMORY_N(vertices, ND_ND);
    RELEASE_MEMORY_N(sort_vert, ND_ND);
    RELEASE_MEMORY(index);
#endif /* RP_NODE */
}


  /*--------------------*/
 /* write_displacement */
/*--------------------*/

DEFINE_ON_DEMAND(write_displacement) {
/* Store previous and new interface node coordinates in UDM's and flags cells adjacent to the boundary.*/
    if (myid == 0) {printf("\nStarted UDF write_displacement.\n"); fflush(stdout);}

    /* declaring variables */
    int thread, n, i, d, compute_node;
    DECLARE_MEMORY_N(disp, real, ND_ND); /* displacement vector */
    DECLARE_MEMORY_N(ids, int, mnpf);
    DECLARE_MEMORY(face_area, real);

#if RP_NODE
    Domain *domain;
    Thread *t, *face_thread, *itf_thread;
    cell_t c;
    face_t face, cell_face;
    Node *node;
    int face_number, node_number, j;
    real heat, vel, src, ds, A_by_es;
    real normal[ND_ND], area[ND_ND], A[ND_ND], es[ND_ND], dr0[ND_ND], dr1[ND_ND], avg_flux[ND_ND], v0[ND_ND];
#endif /* RP_NODE */
    
#if RP_HOST
    char file_name[256];
    FILE *file = NULL;
    timestep = RP_Get_Integer("udf/timestep"); /* host process reads "udf/timestep" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(timestep); /* host process shares timestep variable with nodes */

    for (thread=0; thread<n_threads; thread++) { /* both host and node execute loop over face_threads (= ModelParts) */

#if RP_HOST /* only host process is involved, code not compiled for node */
        sprintf(file_name, "displacement_timestep%i_thread%i.dat",
                timestep, thread_ids[thread]);

        if (NULLP(file = fopen(file_name, "w"))) {
            Error("\nUDF-error: Unable to open %s for writing\n", file_name);
            exit(1);
        }

#if RP_2D
        fprintf(file, "%27s %27s %10s\n",
            "x-disp", "y-disp", "area", "unique-ids");
#else /* RP_2D */
        fprintf(file, "%27s %27s %27s %10s\n",
            "x-disp", "y-disp", "z-disp", "area", "unique-ids");

#endif /* RP_2D */
#endif /* RP_HOST */

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
        domain = Get_Domain(1);
        itf_thread = Lookup_Thread(domain, thread_ids[thread]);

        n = THREAD_N_ELEMENTS_INT(itf_thread); /* get number of faces in this partition of face_thread */

        /* assign memory on compute node that accesses its partition of the face_thread */
        ASSIGN_MEMORY_N(disp, n, real, ND_ND);
        ASSIGN_MEMORY_N(ids, n, int, mnpf);
        ASSIGN_MEMORY(face_area, n, real);

        i = 0;
        /* loop over all faces (tracked by variable "face") in current compute node partition of face_thread */
        begin_f_loop(face, itf_thread) {
            if (i >= n) {Error("\nIndex %i >= array size %i.", i, n);}

            c = F_C0(face, itf_thread); /* adjacent cell */
            t = THREAD_T0(itf_thread); /* adjacent cell thread (F_C0_THREAD(face, itf_thread)) */
            heat = 0;

            /* Allocates memory for gradient and reconstruction gradient macros */
            Alloc_Storage_Vars(domain, SV_T_RG, SV_T_G, SV_NULL);
            /* (Re-)calculate the gradients and store them in the memory allocated before. Have to be called in the correct order. */
            Scalar_Reconstruction(domain, SV_T, -1, SV_T_RG, NULL);
            // Scalar_Derivatives(domain, SV_T, -1, SV_T_G, SV_T_RG, NULL); -> does not work multi-core... & not needed atm

            c_face_loop(c,t,face_number)
            {
                face_thread = C_FACE_THREAD(c,t,face_number);
                cell_face = C_FACE(c,t,face_number);

                /* for cell faces at coupling interface */
                if (THREAD_ID(face_thread) == THREAD_ID(itf_thread)) {

                    /* Incoming heat from liquid side */
                    heat += BOUNDARY_HEAT_FLUX(cell_face, face_thread); // [W]
                    F_AREA(area, cell_face, face_thread);
                    NV_VS(normal, =, area, *, 1.0 / NV_MAG(area)); // normal vector

                    /* store face area of cell */
                    face_area[i] = NV_MAG(area); // [m^2]

                    for (j = 0; j < mnpf; j++) {
                        /* -1 is placeholder, it is usually overwritten, but not always in unstructured grids */
                        ids[j][i] = -1;
                    }

                    /* loop over all nodes in the face */
                    j = 0;
                    f_node_loop(cell_face, face_thread, node_number) {
                        if (j >= mnpf) {Error("\nIndex %i >= array size %i.", j, mnpf);}
                        node = F_NODE(cell_face, face_thread, node_number);
                        ids[j][i] = NODE_DM_ID(node); /* store dynamic mesh node id of current node */
                        j++;
                    }
                }
                else {
                    // Outgoing heat at interior cell faces -- currently not used!
                    INTERIOR_FACE_GEOMETRY(cell_face,face_thread,A,ds,es,A_by_es,dr0,dr1);
                    if (F_C1(cell_face, face_thread) != -1) {
                        NV_VS_VS(avg_flux, =, C_T_RG(c,t), *, 0.5, +, C_T_RG(F_C0(cell_face, face_thread),F_C0_THREAD(cell_face, face_thread)), *, 0.5);
                        NV_VS(v0, =, es, *, A_by_es);
                        heat += C_K_L(c,t)*(C_T(F_C0(cell_face, face_thread),F_C0_THREAD(cell_face, face_thread))-C_T(c,t))*A_by_es/ds + C_K_L(c,t)*(NV_DOT(avg_flux, A)-NV_DOT(avg_flux, v0));
                    }
                }
            }

            if (C_UDMI(c,t,PR_H) >= 0.0) { // during melting
                vel = heat/(NV_MAG(area)*LH*C_R(c,t));
            } else {
                if (C_H(c,t) >= (C_CP(c,t) * (TM - 298.15))) { // transition to melting
                    vel = (C_UDMI(c,t,PR_H) * C_VOLUME_M1(c,t) / dt + heat / C_R(c,t)) / (NV_MAG(area) * LH);
                } else { // subcooled --> no interface velocity
                    vel = 0.0;
                }
            }

            j = 0;
            for (j = 0; j < ND_ND; j++) {
                disp[j][i] = -1.0*normal[j]*vel*dt; // Positive flux results in interface motion opposite to the face normals
            }
            i ++;

            /* Free memory of reconstruction gradients and gradients */
            Free_Storage_Vars(domain, SV_T_RG, SV_T_G, SV_NULL);
        } end_f_loop(face, itf_thread);

        /* assign destination ID compute_node to "node_host" or "node_zero" (these names are known) */
        compute_node = (I_AM_NODE_ZERO_P) ? node_host : node_zero;

        /* send from node to either node_zero, or if the current process is node_zero, to node_host */
        /* the tag argument myid is the ID of the sending node, as the convention is to have the tag argument the same
        as the from argument (that is, the first argument) for receive messages */
        PRF_CSEND_INT(compute_node, &n, 1, myid);

        PRF_CSEND_REAL_N(compute_node, disp, n, myid, ND_ND);
        PRF_CSEND_INT_N(compute_node, ids, n, myid, mnpf);
        PRF_CSEND_REAL(compute_node, face_area, n, myid);

        /* memory can be freed once the data is sent */
        RELEASE_MEMORY_N(disp, ND_ND);
        RELEASE_MEMORY_N(ids, mnpf);
        RELEASE_MEMORY(face_area);

        /* node_zero is the only one that can communicate with host, so it first receives from the other nodes, then
        sends to the host */
        if(I_AM_NODE_ZERO_P){
            compute_node_loop_not_zero(compute_node) { /* loop over all other nodes and receive from each */
                /* the tag argument compute_node is the ID of the sending node, as the convention is to have the tag
                argument the same as the from argument (that is, the first argument) for receive messages */
                PRF_CRECV_INT(compute_node, &n, 1, compute_node);

                /* Once n has been received, the correct amount of memory can be allocated on compute node 0. This
                depends on the partition assigned to the sending compute node. */
                ASSIGN_MEMORY_N(disp, n, real, ND_ND);
                ASSIGN_MEMORY_N(ids, n, int, mnpf);
                ASSIGN_MEMORY(face_area, n, real);

                /* receive the 2D-arrays from the other nodes on node_zero */
                PRF_CRECV_REAL_N(compute_node, disp, n, compute_node, ND_ND);
                PRF_CRECV_INT_N(compute_node, ids, n, compute_node, mnpf);
                PRF_CRECV_REAL(compute_node, face_area, n, compute_node);

                /* Send the variables to the host. Deviating from the tag convention, the message tag is now the
                original non-zero compute node on which the mesh data was stored, even though node 0 does the actual
                communication */
                PRF_CSEND_INT(node_host, &n, 1, compute_node);

                PRF_CSEND_REAL_N(node_host, disp, n, compute_node, ND_ND);
                PRF_CSEND_INT_N(node_host, ids, n, compute_node, mnpf);
                PRF_CSEND_REAL(node_host, face_area, n, compute_node);

                /* once all data has been sent to host, memory on the node can be freed */
                RELEASE_MEMORY_N(disp, ND_ND);
                RELEASE_MEMORY_N(ids, mnpf);
                RELEASE_MEMORY(face_area);
            }
        }
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
        /* loop over all compute nodes (corresponding to the message tags), receive data and append to file for each */
        compute_node_loop(compute_node) {
            PRF_CRECV_INT(node_zero, &n, 1, compute_node);

            /* once n has been received, the correct amount of memory can be allocated on the host */
            ASSIGN_MEMORY_N(disp, n, real, ND_ND);
            ASSIGN_MEMORY_N(ids, n, int, mnpf);
            ASSIGN_MEMORY(face_area, n, real);

            /* receive the 2D-arrays from node_zero */
            PRF_CRECV_REAL_N(node_zero, disp, n, compute_node, ND_ND);
            PRF_CRECV_INT_N(node_zero, ids, n, compute_node, mnpf);
            PRF_CRECV_REAL(node_zero, face_area, n, compute_node);

            for (i = 0; i < n; i++) {
                for (d = 0; d < ND_ND; d++) {
                    fprintf(file, "%27.17e ", disp[d][i]);
                }
                fprintf(file, "%27.17e ", face_area[i]);
                for (d = 0; d < mnpf; d++) {
                    fprintf(file, " %10d", ids[d][i]);
                }
                fprintf(file, "\n");
            }
            /* after files have been appended, the memory can be freed */
            RELEASE_MEMORY_N(disp, ND_ND);
            RELEASE_MEMORY_N(ids, mnpf);
            RELEASE_MEMORY(face_area);
        } /* close compute_node_loop */

        fclose(file);
#endif /* RP_HOST */

    } /* close loop over threads */

    if (myid == 0) {printf("\nFinished UDF write_displacement.\n"); fflush(stdout);}
}


  /*----------------------*/
 /* update_cell_enthalpy */
/*----------------------*/

DEFINE_ON_DEMAND(update_cell_enthalpy)
{
#if RP_NODE
    Domain *d;
    d = Get_Domain(1);
    Thread *t;
    cell_t c;
    real H;

    // loop over all cells
    thread_loop_c(t,d) {
        begin_c_loop(c,t) {
            C_UDMI(c,t,PR_H) = C_H(c,t) - C_CP(c,t)*(TM - 298.15);
        } end_c_loop(c,t)
    }
#endif /* RP_NODE */
}


  /*----------*/
 /* ini_sign */
/*----------*/

DEFINE_ON_DEMAND(ini_sign)
{
#if RP_NODE
    Domain *d;
    d = Get_Domain(1);
    Thread *t;
    cell_t c;
    thread_loop_c(t,d)
    {
        begin_c_loop(c,t) // loop over all cells
        {
            C_UDMI(c,t,SIGN) = 0.0;
        } end_c_loop(c,t)
    }
#endif /* RP_NODE */
}


  /*-----------------*/
 /* udf_mass_source */
/*-----------------*/

DEFINE_SOURCE(udf_mass_source,c,t,dS,eqn)
{
/*Source term for continuity equation, to compensate mass loss or mass gain.*/
real source;

source = C_UDMI(c,t,SIGN)*C_R(c,t)*C_UDMI(c,t,ADJ)*C_UDMI(c,t,D_VOL)/(C_VOLUME(c,t)*dt);
dS[eqn] = C_UDMI(c,t,SIGN)*C_UDMI(c,t,ADJ)*C_UDMI(c,t,D_VOL)/(C_VOLUME(c,t)*dt);
return source;
}


  /*-------------------*/
 /* udf_energy_source */
/*-------------------*/

DEFINE_SOURCE(udf_energy_source,c,t,dS,eqn)
{
/*Source term for energy equation, to account for latent and sensible heat.*/
real source, test;

if (C_UDMI(c,t,ADJ) == 1.0) {
    if (fluid) {
        /* This definition is only valid for fluid during melting! */
        source = (C_UDMI(c,t,SIGN)*C_R(c,t)*(C_CP(c,t)*(TM - 298.15))*C_UDMI(c,t,D_VOL))/(C_VOLUME(c,t)*dt);
        dS[eqn] = 0.0;
    } else {
        /* This definition is only valid for solid during melting! */
        if (C_UDMI(c,t,PR_H) >= 0.0) {
            source = C_UDMI(c,t,SIGN)*C_R(c,t)*C_UDMI(c,t,PR_H)/dt;
        } else {
            source = 0.0;
        }
        // Maybe need to correct with correct cell volumes
        dS[eqn] = 0.0;
    }
} else {
    source = 0.0;
    dS[eqn] = 0.0;
}

return source;
}