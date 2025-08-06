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

int _d; /* don't use in UDFs! (overwritten by functions above) */
int n_threads;
DECLARE_MEMORY(thread_ids, int);
real BODY_DENSITY = 870.0;


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

    sprintf(file_name, "rigid_body.dat");

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
            moi[2] += (BODY_DENSITY * ( (1.0/3.0) * x * x * x * nx + (1.0/3.0) * y * y * y * ny ) * face_area);

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
        MOI[2] = MOI[2] - BODY_DENSITY * volume * (COM[0] * COM[0] + COM[1] * COM[1]);
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

            for (k = 0; k < ND_ND; k++) {
                printf("\nforce[%i] = %lf\n", k, force[k]); fflush(stdout);
                printf("\nr_rel[%i] = %lf\n", k, r_rel[k]); fflush(stdout);
            }

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

    printf("\nSummed at compute node:\n"); fflush(stdout);
    printf("\nF_node_x = %lf\n", F_node[0]); fflush(stdout);
    printf("\nF_node_y = %lf\n", F_node[1]); fflush(stdout);
    printf("\nM_node_z = %lf\n", M_node[2]); fflush(stdout);

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

    printf("\nHOST: Summing from all compute nodes\n"); fflush(stdout);
    printf("\nF_host_x = %lf\n", force_tot[0]); fflush(stdout);
    printf("\nF_host_y = %lf\n", force_tot[1]); fflush(stdout);
    printf("\nM_host_z = %lf\n", moment_tot[2]); fflush(stdout);

    fprintf(file, "%27.17e %27.17e %27.17e\n", force_tot[0], force_tot[1], force_tot[2]);
    fprintf(file, "%27.17e %27.17e %27.17e\n", moment_tot[0], moment_tot[1], moment_tot[2]);
    fclose(file);
#endif /* RP_HOST */

    if (myid == 0) {printf("\nFinished UDF store_rigid_body.\n"); fflush(stdout);}
}