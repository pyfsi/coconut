#include "udf.h"
#include "sg.h"
#include "dynamesh_tools.h"
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

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/* User defined memory -- not available in 3D! */
#define ADJ 0 // flag for cells adjacent to the interface
#define D_VOL 1 // swept volume during move_nodes operation
#define SIGN 2 // defines wether the cell is on the growing or shrinking side during phase change

/* User defined node memory -- not available in 3D! */
#define N_ID 0 // u-d node memory where initial node ids are saved
#define PR_X 1 // previous x-coordinate of first node in face
#define PR_Y 2 // previous y-coordinate of first node in face
#define NEW_X 3 // new x-coordinate of first node in face
#define NEW_Y 4 // new y-coordinate of first node in face
#define OS_PR_X 5 // new x-coordinate of first node in face
#define OS_PR_Y 6 // new y-coordinate of first node in face

/* global variables */
#define mnpf |MAX_NODES_PER_FACE|
real dt = |TIME_STEP_SIZE|;
real TM = |MELT_TEMP|;
real HM = |MELT_ENTHALPY|;
real rho_s = |SOLID_DENSITY|;
int TS_START = |TIME_STEP_START|;
bool unsteady = |UNSTEADY|;
int _d; /* don't use in UDFs! (overwritten by functions above) */
int n_threads, n_os_threads;
DECLARE_MEMORY(thread_ids, int);
DECLARE_MEMORY(overset_thread_ids, int);
int timestep = 0;
int iteration = 0;
int rb_it = 0;


  /*----------------*/
 /* get_thread_ids */
/*----------------*/

DEFINE_ON_DEMAND(get_thread_ids) {
    /* Read in thread ids, which is the name of the ModelParts in Fluent. Should be called early on. */

#if RP_HOST /* only host process is involved, code not compiled for node */
    int k;
    FILE *file;
    FILE *file2;
    file = fopen("bcs.txt", "r");
    file2 = fopen("bcs_overset.txt", "r");

    if (file == NULL || file2 == NULL) {
        Error("\nUDF-error: could not open bcs.txt or bcs_overset.txt\n");
        exit(1);
    }

    fscanf(file, "%i", &n_threads);
    fscanf(file2, "%i", &n_os_threads);

    if (n_threads != n_os_threads) {
        Error("\nUDF-error: mismatch in number of boundary and overset threads (%i vs %i)\n", n_threads, n_os_threads);
        fclose(file);
        fclose(file2);
        exit(1);
    }
#endif /* RP_HOST */

    /* Broadcast sizes to nodes */
    host_to_node_int_1(n_threads);
    host_to_node_int_1(n_os_threads);

    /* Allocate arrays */
    ASSIGN_MEMORY(thread_ids, n_threads, int);
    ASSIGN_MEMORY(overset_thread_ids, n_os_threads, int);

#if RP_HOST /* only host process is involved, code not compiled for node */
    for (k = 0; k < n_threads; k++) {
        fscanf(file, "%*s %i", &thread_ids[k]);
        fscanf(file2, "%*s %i", &overset_thread_ids[k]);
    }
    fclose(file);
    fclose(file2);
#endif /* RP_HOST */

    /* Broadcast arrays to nodes */
    host_to_node_int(thread_ids, n_threads); /* thread_ids (read from file on host) communicated to nodes */
    host_to_node_int(overset_thread_ids, n_os_threads);
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
                if (timestep == TS_START) {
                    N_UDMI(node,N_ID) = NODE_DM_ID(node); /* store initial dynamic mesh node id in a node UDM */
                    for (d = 0; d < ND_ND; d++) {
                        node_coords[d][i_n] = NODE_COORD(node)[d];
                        /* initialise node UDM's with melting only node locations */
                        N_UDMI(node,PR_X) = NODE_X(node);
                        N_UDMI(node,NEW_X) = NODE_X(node);
                        N_UDMI(node,PR_Y) = NODE_Y(node);
                        N_UDMI(node,NEW_Y) = NODE_Y(node);
                    }
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


  /*------------------*/
 /* store_rigid_body */
/*------------------*/

DEFINE_ON_DEMAND(store_rigid_body) {
    /* The rigid body properties of the solid phase are derived from the outser surface as seen from the liquid solver.
    In a fist step, the volume, center of mass (x,y) and moment of inertia (z) are calculated.
    These are used in a second step to calculate the force- and moment balance.
    Finally, these values are stored in a file rigid_body.dat. */
    if (myid == 0) {printf("\nStarted UDF store_rigid_body.\n"); fflush(stdout);}

    /* declaring variables */
    int thread, n, k, i, compute_node;
    /* these memory locations will be used to communicate between host- and node processes */
    DECLARE_MEMORY(vol, real);
    DECLARE_MEMORY(com, real);
    DECLARE_MEMORY(moi, real);
    DECLARE_MEMORY(COM_host, real);
    DECLARE_MEMORY(F_node, real);
    DECLARE_MEMORY(M_node, real);

#if RP_HOST /* only host process is involved, code not compiled for node */
    char file_name[256];
    FILE *file = NULL;
    timestep = RP_Get_Integer("udf/timestep"); /* host process reads "udf/timestep" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(timestep); /* host process shares timestep variable with nodes */

#if RP_HOST /* only host process is involved, code not compiled for node */
    sprintf(file_name, "rigid_body_timestep%i.dat", timestep);

    if (NULLP(file = fopen(file_name, "w"))) {
        Error("\nUDF-error: Unable to open %s for writing\n", file_name);
        exit(1);
    }

    fprintf(file, "Order: volume, COM-xyz, MOI-xyz, force-xyz, moment-xyz\n");
#endif /* RP_HOST */

#if RP_NODE /* only compute_nodes are involved, code not compiled for host */
    Domain *domain;
    Thread *face_thread;
    face_t face;

    real face_area;
    real area[ND_ND];
    real face_centroid[ND_ND];
    real unit_normal[ND_ND];
    real x, y, nx, ny;

    /* Contribution of each compute node to the total volume, center of mass & moment of inertia */
    ASSIGN_MEMORY(vol, 1, real); /* Volume of solid phase */
    ASSIGN_MEMORY(com, 3, real); /* Numerator for center of mass */
    ASSIGN_MEMORY(moi, 3, real); /* Numerator for moment of inertia */

    /* Initialize accumulators and relevant components */
    vol[0] = 0.0;
    for (k = 0; k < 3; k++) {
        com[k] = 0.0;
        moi[k] = 0.0;
    }

    for (thread=0; thread<n_threads; thread++) { /* only node executes loop over face_threads (= ModelParts) */
        domain = Get_Domain(1);
        face_thread = Lookup_Thread(domain, thread_ids[thread]);

        n = THREAD_N_ELEMENTS_INT(face_thread); /* get number of faces in this partition of face_thread */

        i = 0;
        /* loop over all faces (tracked by variable "face") in current compute_node partition of face_thread */
        begin_f_loop(face, face_thread) {
            if (i >= n) {Error("\nIndex %i >= array size %i.", i, n);}

            F_CENTROID(face_centroid, face, face_thread);
            F_AREA(area, face, face_thread);
            NV_VS(unit_normal, =, area, *, 1.0 / NV_MAG(area)); // normal vector
            face_area = NV_MAG(area); // [m^2]

            x = face_centroid[0];
            y = face_centroid[1];
            nx = -1.0 * unit_normal[0]; // normal vectors points inwards the solid volume --> has to be inverted
            ny = -1.0 * unit_normal[1]; // normal vectors points inwards the solid volume --> has to be inverted

            // Volume calculation using divergence theorem (in 2D)
            vol[0] += (1.0/2.0) * (x * nx + y * ny) * face_area;

            // Numerator for center of mass (x, y)
            com[0] += 0.5 * x * x * nx * face_area;
            com[1] += 0.5 * y * y * ny * face_area;

            // Numerator for Moment of Inertia about Z-axis (about origin)
            moi[2] += (rho_s * ( (1.0/3.0) * x * x * x * nx + (1.0/3.0) * y * y * y * ny ) * face_area);

            i++;
        } end_f_loop(face, face_thread);
    } /* close loop over threads */

    /* assign destination ID compute_node to "node_host" or "node_zero" (these names are known) */
    compute_node = (I_AM_NODE_ZERO_P) ? node_host : node_zero;

    /* send from node to either node_zero, or if the current process is node_zero, to node_host */
    PRF_CSEND_REAL(compute_node, vol, 1, myid);
    PRF_CSEND_REAL(compute_node, com, 3, myid);
    PRF_CSEND_REAL(compute_node, moi, 3, myid);

    /* memory can be freed once the data is sent */
    RELEASE_MEMORY(vol);
    RELEASE_MEMORY(com);
    RELEASE_MEMORY(moi);

    /* node_zero is the only one that can communicate with host, so it first receives from the other nodes, then
    sends to the host */
    if(I_AM_NODE_ZERO_P){
        compute_node_loop_not_zero(compute_node) { /* loop over all other nodes and receive from each */

            ASSIGN_MEMORY(vol, 1, real);
            ASSIGN_MEMORY(com, 3, real);
            ASSIGN_MEMORY(moi, 3, real);

            /* receive the 1D-arrays from the other nodes on node_zero */
            PRF_CRECV_REAL(compute_node, vol, 1, compute_node);
            PRF_CRECV_REAL(compute_node, com, 3, compute_node);
            PRF_CRECV_REAL(compute_node, moi, 3, compute_node);

            /* Send the variables to the host. Deviating from the tag convention, the message tag is now the
            original non-zero compute_node on which the mesh data was stored, even though node 0 does the actual
            communication */
            PRF_CSEND_REAL(node_host, vol, 1, compute_node);
            PRF_CSEND_REAL(node_host, com, 3, compute_node);
            PRF_CSEND_REAL(node_host, moi, 3, compute_node);

            /* once all data has been sent to host, memory on the node can be freed */
            RELEASE_MEMORY(vol);
            RELEASE_MEMORY(com);
            RELEASE_MEMORY(moi);
        }
    }
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
    /* Total volume, center of mass & moment of inertia is a summation of each compute node's contribution */
    real volume = 0;
    real MOI[3], COM[3];
    ASSIGN_MEMORY(COM_host, 3, real);

    for (k = 0; k < 3; k++) {
        MOI[k] = 0.0;
        COM[k] = 0.0;
        COM_host[k] = 0.0;
    }

    /* loop over all compute_nodes (corresponding to the message tags), receive data and append to file for each */
    compute_node_loop(compute_node) {
        ASSIGN_MEMORY(vol, 1, real);
        ASSIGN_MEMORY(com, 3, real);
        ASSIGN_MEMORY(moi, 3, real);
        if (!vol || !com || !moi) {
            Error("HOST_UDF_ERROR: Failed to assign memory for receiving buffers for node %d.\n", compute_node);
            exit(1);
        }

        PRF_CRECV_REAL(node_zero, vol, 1, compute_node);
        PRF_CRECV_REAL(node_zero, com, 3, compute_node);
        PRF_CRECV_REAL(node_zero, moi, 3, compute_node);

        /* sum the contributions of each compute node's part of the partitioned interface */
        volume += vol[0];
        for (k = 0; k < 3; k++) {
            COM[k] += com[k];
            MOI[k] += moi[k];
        }

        /* after files have been appended, the memory can be freed */
        RELEASE_MEMORY(vol);
        RELEASE_MEMORY(com);
        RELEASE_MEMORY(moi);
    } /* close compute_node_loop */

    if (fabs(volume) > 1e-9) { // Check for non-zero volume
        // Divide by volume to get center of mass coordinates
        COM[0] = COM[0] / volume;
        COM[1] = COM[1] / volume;
        // Apply Parallel Axis Theorem for moment of inertia through center of mass
        MOI[2] = MOI[2] - rho_s * volume * (COM[0] * COM[0] + COM[1] * COM[1]);
    } else {
        printf("HOST_UDF_WARNING: Total volume is near zero, COM and MOI about COM might be inaccurate.\n"); fflush(stdout);
    }

    for (k = 0; k < 3; k++) {
        COM_host[k] = COM[k];
    }

    // send to node zero to calculate moments around te correct center of mass
    PRF_CSEND_REAL(node_zero, COM_host, 3, node_host);
    RELEASE_MEMORY(COM_host);

    fprintf(file, "%27.17e\n", volume);
    fprintf(file, "%27.17e %27.17e %27.17e\n", COM[0], COM[1], COM[2]);
    fprintf(file, "%27.17e %27.17e %27.17e\n", MOI[0], MOI[1], MOI[2]);
#endif /* RP_HOST */

#if RP_NODE /* Force and moment balance is only performed at node zero */
    Domain *domain_2;
    Thread *tf;
    face_t f;

    /* face_area, area[ND_ND], unit_normal[ND_ND] & face_centroid[ND_ND] already earlier declared for nodes */
    real traction[ND_ND];

    real force[3], moment[3]; /* Force and moment per cell face */
    real r_rel[3]; /* (x,y,z) vector of face centroid relative to center of mass */

    ASSIGN_MEMORY(F_node, 3, real); /* Total force per node (summed over all faces) */
    ASSIGN_MEMORY(M_node, 3, real); /* Total moment per node (summed over all faces) */

    if(I_AM_NODE_ZERO_P){ /* only node zero communicates with the host */
        /* receive center of mass from host process */
        ASSIGN_MEMORY(COM_host, 3, real);
        PRF_CRECV_REAL(node_host, COM_host, 3, node_host);

        compute_node_loop_not_zero(compute_node) {
            PRF_CSEND_REAL(compute_node, COM_host, 3, node_zero);
        }
    }

    if(! I_AM_NODE_ZERO_P){ /* other compute node receive from node zero */
        ASSIGN_MEMORY(COM_host, 3, real);
        PRF_CRECV_REAL(node_zero, COM_host, 3, node_zero);
    }

    for (k = 0; k < 3; k++) {
        force[k] = 0.0;
        moment[k] = 0.0;
        F_node[k] = 0.0;
        M_node[k] = 0.0;
        r_rel[k] = 0.0;
    }

    for (thread=0; thread<n_threads; thread++) { /* only node executes loop over face_threads (= ModelParts) */
        domain_2 = Get_Domain(1);
        tf = Lookup_Thread(domain_2, thread_ids[thread]);

        n = THREAD_N_ELEMENTS_INT(tf); /* get number of faces in this partition of face_thread */

        i = 0;
        /* loop over all faces (tracked by variable "face") in current compute_node partition of face thread 'tf' */
        begin_f_loop(f, tf) {
            if (i >= n) {Error("\nIndex %i >= array size %i.", i, n);}

            F_AREA(area, f, tf);
            F_CENTROID(face_centroid, f, tf);
            face_area = NV_MAG(area); // [m^2]
            NV_VS(traction, =, F_STORAGE_R_N3V(f, tf, SV_WALL_SHEAR), *, -1.0 / face_area);
            NV_VS(unit_normal, =, area, *, 1.0 / face_area); // normal vector

            // force integral (pressure and traction)
            for (k = 0; k < ND_ND; k++) {
                force[k] = F_P(f, tf) * area[k] + traction[k] * face_area;
            }

            /* Calculate relative position vector (r - r_com) */
            NV_VV(r_rel, =, face_centroid, -, COM_host);

            /* Moment as cross product of force and relative position to center of mass */
            moment[2] = NV_CROSS_Z(r_rel, force);
            // moment[2] = r_rel[0] * force[1] - r_rel[1] * force[0];

            F_node[0] += force[0];
            F_node[1] += force[1];
            M_node[2] += moment[2];

            for (k = 0; k < 3; k++) {
                force[k] = 0.0;
                moment[k] = 0.0;
                r_rel[k] = 0.0;
            }

            i++;
        } end_f_loop(face, face_thread); /* close face loop */
    } /* close loop over threads */

    /* COM memory can be freed once the data is used */
    RELEASE_MEMORY(COM_host);

    /* assign destination ID compute_node to "node_host" or "node_zero" (these names are known) */
    compute_node = (I_AM_NODE_ZERO_P) ? node_host : node_zero;

    /* send from node to either node_zero, or if the current process is node_zero, to node_host */
    PRF_CSEND_REAL(compute_node, F_node, 3, myid);
    PRF_CSEND_REAL(compute_node, M_node, 3, myid);

    /* memory can be freed once the data is sent */
    RELEASE_MEMORY(F_node);
    RELEASE_MEMORY(M_node);

    /* node_zero is the only one that can communicate with host, so it first receives from the other nodes, then
    sends to the host */
    if(I_AM_NODE_ZERO_P){
        compute_node_loop_not_zero(compute_node) { /* loop over all other nodes and receive from each */

            ASSIGN_MEMORY(F_node, 3, real);
            ASSIGN_MEMORY(M_node, 3, real);

            /* receive the 1D-arrays from the other nodes on node_zero */
            PRF_CRECV_REAL(compute_node, F_node, 3, compute_node);
            PRF_CRECV_REAL(compute_node, M_node, 3, compute_node);

            /* Send the variables to the host. Deviating from the tag convention, the message tag is now the
            original non-zero compute_node on which the mesh data was stored, even though node 0 does the actual
            communication */
            PRF_CSEND_REAL(node_host, F_node, 3, compute_node);
            PRF_CSEND_REAL(node_host, M_node, 3, compute_node);

            /* once all data has been sent to host, memory on the node can be freed */
            RELEASE_MEMORY(F_node);
            RELEASE_MEMORY(M_node);
        }
    }
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
    real force_tot[3], moment_tot[3]; /* Force and moment summed over each compute node */

    for (k = 0; k < 3; k++) {
        force_tot[k] = 0.0;
        moment_tot[k] = 0.0;
    }

    compute_node_loop(compute_node) {
        ASSIGN_MEMORY(F_node, 3, real);
        ASSIGN_MEMORY(M_node, 3, real);
        if (!F_node || !M_node) {
            Error("HOST_UDF_ERROR: Failed to assign memory for receiving buffers for node %d.\n", compute_node);
            exit(1);
        }

        PRF_CRECV_REAL(node_zero, F_node, 3, compute_node);
        PRF_CRECV_REAL(node_zero, M_node, 3, compute_node);

        /* sum the contributions of each compute node's part of the partitioned interface */
        for (k = 0; k < 3; k++) {
            force_tot[k] += F_node[k];
            moment_tot[k] += M_node[k];
        }

        /* after files have been appended, the memory can be freed */
        RELEASE_MEMORY(F_node);
        RELEASE_MEMORY(M_node);
    } /* close compute_node_loop */

    fprintf(file, "%27.17e %27.17e %27.17e\n", force_tot[0], force_tot[1], force_tot[2]);
    fprintf(file, "%27.17e %27.17e %27.17e\n", moment_tot[0], moment_tot[1], moment_tot[2]);
    fclose(file);
#endif /* RP_HOST */

    if (myid == 0) {printf("\nFinished UDF store_rigid_body.\n"); fflush(stdout);}
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
            // C_UDMI(F_C0(face, face_thread),THREAD_T0(face_thread),SIGN) = 1.0;

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
                ids[j][i] = N_UDMI(node,N_ID); /* store dynamic mesh node id of current node */
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
        if (timestep == TS_START) {
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
                            id = N_UDMI(node,N_ID);
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


  /*------------*/
 /* move_nodes */
/*------------*/

DEFINE_GRID_MOTION(move_nodes, domain, dynamic_thread, time, dtime) {
    /* UDF that can be assigned to a deformable wall in Fluent. It will read the updated positions of the mesh nodes
    (vertices) and apply them. As only the compute nodes have access to their own partition of the mesh, each one is
    responsible for reading the file containing the updated coordinates and selecting the correct row. Furthermore,
    the node coordinates in the solid domain, which result from phase change only, are written to node UDM's to allow
    calculation of volume change due to phase change */
    if (myid == 0) {printf("\nStarted UDF move_nodes.\n"); fflush(stdout);}
    char file_name[256];
    char file_name_2[256];
    Thread *face_thread = DT_THREAD(dynamic_thread); /* face_thread to which UDF is assigned in Fluent */
    int thread_id = THREAD_ID(face_thread);

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    face_t face;
    Node *node;
    int i, d, n, n_pc, node_number;
    DECLARE_MEMORY_N(coords, real, ND_ND);
    DECLARE_MEMORY(ids, int);
    DECLARE_MEMORY_N(coords_pc, real, ND_ND);
    DECLARE_MEMORY(ids_pc, int);
    FILE *file = NULL;
    FILE *file2 = NULL;
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
    timestep = RP_Get_Integer("udf/timestep"); /* host process reads "udf/timestep" from Fluent (nodes cannot) */
    iteration = RP_Get_Integer("udf/iteration"); /* host process reads "udf/iteration" from Fluent (nodes cannot) */
    rb_it = RP_Get_Integer("udf/rb_it"); /* host process reads "udf/rb_it" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(timestep); /* host process shares timestep variable with nodes */
    host_to_node_int_1(iteration); /* host process shares iteration variable with nodes */
    host_to_node_int_1(rb_it); /* host process shares rb_it variable with nodes */

#if RP_HOST /* only host process is involved, code not compiled for node */
    /* File with real node locations */
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

#if RP_HOST /* only host process is involved, code not compiled for node */
    /* File with node locations due to phase change only */
    sprintf(file_name_2, "nodes_pc_timestep%i_thread%i.dat",
            timestep, thread_id);
    host_to_node_sync_file(file_name_2); /* send file to the compute nodes */
#else
    sprintf(file_name_2, "|TMP_DIRECTORY_NAME|/nodes_pc_timestep%i_thread%i.dat",
            timestep, thread_id);
    host_to_node_sync_file("|TMP_DIRECTORY_NAME|"); /* receive file on compute nodes and store in a temporary folder */
#endif /* RP_HOST */

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    /* Read and assign the nodes_update file for the current position of the nodes */
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

    /* Read and assign the nodes_pc file for the position of the nodes in the solid solver */
    if (NULLP(file2 = fopen(file_name_2, "r"))) {
        Error("\nUDF-error: Unable to open %s for reading\n", file_name_2);
        exit(1);
    }

    fscanf(file2, "%i", &n_pc); /* number of nodes in second file */

    if (n_pc != n) {
        Error("\nUDF-error: Node count mismatch (%d vs %d) between %s and %s\n",
              n, n_pc, file_name, file_name_2);
        exit(1);
    }

    ASSIGN_MEMORY_N(coords_pc, n, real, ND_ND);
    ASSIGN_MEMORY(ids_pc, n, int);

    for (i=0; i < n; i++) {
        for (d = 0; d < ND_ND; d++) {
            fscanf(file2, "%lf", &coords_pc[d][i]); /* read coordinates from file */
        }
        fscanf(file2, "%i", &ids_pc[i]); /* read node ids from file */
    }

    fclose(file2);

    /* Compare IDs */
    for (i = 0; i < n; i++) {
        if (ids[i] != ids_pc[i]) {
            Error("\nUDF-error: Node ID mismatch at index %d (%d != %d)\n",
                  i, ids[i], ids_pc[i]);
            exit(1);
        }
    }

    begin_f_loop(face, face_thread) { /* loop over all faces in face_thread */
        f_node_loop(face, face_thread, node_number) { /* loop over all nodes in current face */
            node = F_NODE(face, face_thread, node_number); /* get global face ndoe index from local node index */
            if NODE_POS_NEED_UPDATE(node) { /* only execute if position has not been updated yet (efficiency) */
                int found_node = 0;
                for (i=0; i < n; i++) { /* loop over all lines to find the correct node */
                    if (N_UDMI(node,N_ID) == ids[i]) { /* correct node has the same dynamic mesh node id */
                        for (d = 0; d < ND_ND; d++) {
                            NODE_COORD(node)[d] = coords[d][i]; /* modify node coordinates */
                            /* Update melting only node locations */
                            if (d == 0) {
                                if ((iteration == 1) && (rb_it == 1)) {
                                    N_UDMI(node,PR_X) = N_UDMI(node,NEW_X); /* Update previous node position only when new time step occurs */
                                }
                                N_UDMI(node,NEW_X) = coords_pc[d][i];
                            } else if (d == 1) {
                                if ((iteration == 1) && (rb_it == 1)) {
                                    N_UDMI(node,PR_Y) = N_UDMI(node,NEW_Y); /* Update previous node position only when new time step occurs */
                                }
                                N_UDMI(node,NEW_Y) = coords_pc[d][i];
                            }
                        }
                        NODE_POS_UPDATED(node); /* flag node as updated */
                        found_node = 1;
                        break; /* break loop for efficiency */
                    }
                }
                if (!found_node) {
                    printf("\nUDF-warning: No match for node id %i\n", N_UDMI(node,N_ID));
                    Error("\nUDF-error: No match for node id %i\n", N_UDMI(node,N_ID));
                    exit(1);
                }
            }
        }
    } end_f_loop(face, face_thread);

    RELEASE_MEMORY_N(coords, ND_ND);
    RELEASE_MEMORY(ids);
    RELEASE_MEMORY_N(coords_pc, ND_ND);
    RELEASE_MEMORY(ids_pc);

    if (myid == 0) {
        sprintf(file_name, "|TMP_DIRECTORY_NAME|/nodes_update_timestep%i_thread%i.dat",
                timestep-1, thread_id);
        remove(file_name);
        sprintf(file_name_2, "|TMP_DIRECTORY_NAME|/nodes_pc_timestep%i_thread%i.dat",
                timestep-1, thread_id);
        remove(file_name_2);
    }
#endif /* RP_NODE */

    if (myid == 0) {printf("\nFinished UDF move_nodes.\n"); fflush(stdout);}
}


  /*-----------------------*/
 /* move_overset_boundary */
/*-----------------------*/

DEFINE_GRID_MOTION(move_overset_boundary, domain, dynamic_thread, time, dtime) {
    /* UDF that can be assigned to the overset boundary in a component mesh in Fluent.
    It will read the updated (translational and rotational) bulk velocity and update the nodes accordingly. */
    if (myid == 0) {printf("\nStarted UDF move_overset_boundary.\n"); fflush(stdout);}
    char file_name[256];
    Thread *face_thread = DT_THREAD(dynamic_thread); /* face_thread to which UDF is assigned in Fluent */
    int thread_id = THREAD_ID(face_thread);

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    face_t face;
    Node *node;
    int i, d, n, node_number;
    DECLARE_MEMORY(v_trans, real);    /* Translational velocity */
    DECLARE_MEMORY(v_rot, real);      /* Rotational velocity */
    DECLARE_MEMORY(com, real);        /* Center of mass */
    DECLARE_MEMORY(r_new, real);      /* New node coordinates */
    DECLARE_MEMORY(r_prev, real);      /* New node coordinates */
    real disp[3], cross[3], r_rel[3]; /* helper arrays */
    FILE *file = NULL;
    int skip_motion = 0;   /* flag to skip node update */
#endif /* RP_NODE */

#if RP_HOST /* only host process is involved, code not compiled for node */
    timestep = RP_Get_Integer("udf/timestep"); /* host process reads "udf/timestep" from Fluent (nodes cannot) */
    iteration = RP_Get_Integer("udf/iteration"); /* host process reads "udf/iteration" from Fluent (nodes cannot) */
    rb_it = RP_Get_Integer("udf/rb_it"); /* host process reads "udf/rb_it" from Fluent (nodes cannot) */
#endif /* RP_HOST */

    host_to_node_int_1(timestep); /* host process shares timestep variable with nodes */
    host_to_node_int_1(iteration); /* host process shares iteration variable with nodes */
    host_to_node_int_1(rb_it); /* host process shares rb_it variable with nodes */

#if RP_HOST /* only host process is involved, code not compiled for node */
    sprintf(file_name, "RB_update_timestep%i_thread%i.dat",
            timestep, thread_id);
    host_to_node_sync_file(file_name); /* send file to the compute nodes */
#else
    struct stat st = {0};
    /* create a temporary directory if it does not exist yet, this is needed as multiple machines can be involved */
    if (stat("|TMP_DIRECTORY_NAME|", &st) == -1) {
        mkdir("|TMP_DIRECTORY_NAME|", 0700);
    }
    sprintf(file_name, "|TMP_DIRECTORY_NAME|/RB_update_timestep%i_thread%i.dat",
            timestep, thread_id);
    host_to_node_sync_file("|TMP_DIRECTORY_NAME|");  /* receive file on compute nodes and store in a temporary folder */
#endif /* RP_HOST */

#if RP_NODE /* only compute nodes are involved, code not compiled for host */
    if (NULLP(file = fopen(file_name, "r"))) {
        Error("\nUDF-error: Unable to open %s for reading\n", file_name);
        exit(1);
    }

    ASSIGN_MEMORY(v_trans, 3, real);
    ASSIGN_MEMORY(v_rot, 3, real);
    ASSIGN_MEMORY(com, 3, real);

     /* read from file */
    for (d = 0; d < 3; d++) {fscanf(file, "%lf", &v_trans[d]);}
    for (d = 0; d < 3; d++) {fscanf(file, "%lf", &v_rot[d]);}
    for (d = 0; d < 3; d++) {fscanf(file, "%lf", &com[d]);}
    fclose(file);

    /* check if both v_trans and v_rot are zero vectors (within tolerance) */
    {
        real tol = 1.0e-15;
        real sumsq = 0.0;
        for (d = 0; d < 3; d++) {
            sumsq += v_trans[d]*v_trans[d] + v_rot[d]*v_rot[d];
        }
        if (sumsq < tol) {
            skip_motion = 1;
        }
    }

    if (!skip_motion) {
        ASSIGN_MEMORY(r_new, 3, real);
        ASSIGN_MEMORY(r_prev, 3, real);
        begin_f_loop(face, face_thread) { /* loop over all faces in face_thread */
            f_node_loop(face, face_thread, node_number) { /* loop over all nodes in current face */
                node = F_NODE(face, face_thread, node_number); /* get global face ndoe index from local node index */
                if NODE_POS_NEED_UPDATE(node) { /* only execute if position has not been updated yet (efficiency) */
                    /* Initialise vectors */
                    for (d = 0; d < 3; d++) {
                        r_new[d] = 0.0;
                        disp[d] = 0.0;
                        cross[d] = 0.0;
                        r_rel[d] = 0.0;
                        r_prev[d] = 0.0;
                        if ((iteration == 1) && (rb_it == 1)) {
                            if (d == 0) N_UDMI(node, OS_PR_X) = NODE_X(node);
                            else if (d == 1) N_UDMI(node, OS_PR_Y) = NODE_Y(node);
                        }
                    }

                    r_prev[0] = N_UDMI(node, OS_PR_X);
                    r_prev[1] = N_UDMI(node, OS_PR_Y);

                    /* rigid body kinematics */
                    NV_VV(r_rel, =, r_prev, -, com);
                    NV_CROSS(cross, v_rot, r_rel);
                    NV_VS_VS(disp, =, v_trans, *, dt, +, cross, *, dt);
                    NV_VV(r_new, =, r_prev, +, disp);

                    /*
                    if ((NODE_COORD(node)[0] < 0.00001) && (NODE_COORD(node)[0] > -0.00001)) {
                        if (NODE_COORD(node)[1] > 0.04) {
                            for (d = 0; d < 3; d++) {
                                printf("\nr_rel[%i] = %lf\n", d, r_rel[d]); fflush(stdout);
                                printf("\nr_new[%i] = %lf\n", d, r_new[d]); fflush(stdout);
                                printf("\ndisp[%i] = %lf\n", d, disp[d]); fflush(stdout);
                                printf("\ncross[%i] = %lf\n", d, cross[d]); fflush(stdout);
                                printf("\ncom[%i] = %lf\n", d, com[d]); fflush(stdout);
                                printf("\nv_rot[%i] = %lf\n", d, v_rot[d]); fflush(stdout);
                                printf("\nv_trans[%i] = %lf\n", d, v_trans[d]); fflush(stdout);
                            }
                        }
                    }
                    */

                    for (d = 0; d < ND_ND; d++) {
                        NODE_COORD(node)[d] = r_new[d]; /* update node coordinates */
                    }
                    NODE_POS_UPDATED(node); /* flag node as updated */
                }
            }
        } end_f_loop(face, face_thread);
        RELEASE_MEMORY(r_new);
        RELEASE_MEMORY(r_prev);
    }

    RELEASE_MEMORY(v_trans);
    RELEASE_MEMORY(v_rot);
    RELEASE_MEMORY(com);

    if (myid == 0) {
        sprintf(file_name, "|TMP_DIRECTORY_NAME|/RB_update_timestep%i_thread%i.dat",
                timestep-1, thread_id);
        remove(file_name);}
#endif /* RP_NODE */

    if (myid == 0) {printf("\nFinished UDF move_overset_boundary.\n"); fflush(stdout);}
}


  /*--------------*/
 /* set_adjacent */
/*--------------*/

DEFINE_ON_DEMAND(set_adjacent) {
/* Flags cells adjacent to the boundary.*/

#if RP_NODE
    Domain *domain;
    Thread *cell_thread, *face_thread;
    cell_t cell;
    int face_number, node_number, surface;

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
    Thread *cell_thread, *face_thread;
    cell_t cell;
    face_t face;
    Node *node;
    int face_number, node_number, surface;
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

    domain = Get_Domain(1);
    thread_loop_c(cell_thread,domain)
    {
        begin_c_loop(cell,cell_thread) //loop over all cells
        {
            zero_vol = false;
            if (C_UDMI(cell,cell_thread,ADJ) == 1.0) {
                c_face_loop(cell,cell_thread,face_number)
                    {
                        face_thread = C_FACE_THREAD(cell,cell_thread,face_number);
                        for (surface = 0; surface < n_threads; surface++)
                        {
                            if (THREAD_ID(face_thread) == thread_ids[surface])
                            {
                                face = C_FACE(cell,cell_thread,face_number);
                                i = 0;
                                f_node_loop(face, face_thread, node_number) {
                                    if (i > 1) {
                                        Error("\nUDF-error: A face with more than two nodes is impossible in 2D\n");
                                        exit(1);
                                    }

                                    node = F_NODE(face, face_thread, node_number);
                                    vertices[0][i] = N_UDMI(node,PR_X);
                                    vertices[0][i+2] = N_UDMI(node,NEW_X);
                                    vertices[1][i] = N_UDMI(node,PR_Y);
                                    vertices[1][i+2] = N_UDMI(node,NEW_Y);
                                    i ++;
                                }
                            }
                        }
                    }

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


  /*---------*/
 /* ini_udm */
/*---------*/

DEFINE_ON_DEMAND(ini_udm)
{
#if RP_NODE
    Domain *d;
    d = Get_Domain(1);
    Thread *t;
    cell_t c;
    thread_loop_c(t,d) {
        begin_c_loop(c,t) { // loop over all cells
            C_UDMI(c,t,SIGN) = 1.0; // currently only melting assumed
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
real source, rho;

if (rho_s == 0.0) {
    rho = C_R(c,t);
} else {
    rho = rho_s;
}

source = C_UDMI(c,t,SIGN)*C_UDMI(c,t,ADJ)*rho*C_UDMI(c,t,D_VOL)/(C_VOLUME(c,t)*dt);
dS[eqn] = 0.0;

return source;
}


  /*-------------------*/
 /* udf_energy_source */
/*-------------------*/

DEFINE_SOURCE(udf_energy_source,c,t,dS,eqn)
{
/*Source term for energy equation, to account for sensible heat.*/
real source, rho;

if (rho_s == 0.0) {
    rho = C_R(c,t);
} else {
    rho = rho_s;
}

// This definition is only valid for fluid during melting!
if (HM != 0.0) {
    source = (C_UDMI(c,t,ADJ)*C_UDMI(c,t,SIGN)*rho*HM*C_UDMI(c,t,D_VOL))/(C_VOLUME(c,t)*dt);
} else {
    source = (C_UDMI(c,t,ADJ)*C_UDMI(c,t,SIGN)*rho*C_CP(c,t)*(TM - 298.15)*C_UDMI(c,t,D_VOL))/(C_VOLUME(c,t)*dt);
}
dS[eqn] = 0.0;

return source;
}