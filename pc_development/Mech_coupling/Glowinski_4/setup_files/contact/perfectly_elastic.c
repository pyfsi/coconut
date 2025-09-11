#include "udf.h"
#include "dpm.h"
#include "sg_mphase.h"
#include "mem.h"
#include "dynamesh_tools.h"

/************************************************************\
* 2-degree of freedom equation of motion compiled UDF *
\************************************************************/
DEFINE_CONTACT(contact_props, dt, contacts)
{
    Objp *o;
    face_t face;
    Thread *thread;
    Domain *domain = NULL;
    Dynamic_Thread *ndt = NULL;
    int tid, nid, n_faces;

    real v0dotn1, v1dotn0;
    real nc_mag, norm0_mag, norm1_mag;
    real N3V_VEC (vel_rel);
    real N3V_VEC (nc), N3V_VEC (nctmp);
    real N3V_VEC (xc), N3V_VEC (xctmp);
    real N3V_VEC (vel0), N3V_VEC (omega0), N3V_VEC (theta0), N3V_VEC (norm0);
    real N3V_VEC (vel1), N3V_VEC (omega1), N3V_VEC (theta1), N3V_VEC (norm1);

    if (!Data_Valid_P())
    {
        return;
    }

    /* Define a common contact point / plane */
    N3V_S (nc, =, 0.0);
    N3V_S (xc, =, 0.0);

    /* Fetch current thread ID */
    tid = THREAD_ID (DT_THREAD (dt));
    nid = -1;
    n_faces = 0;
    loop (o, contacts)
    {
        face = O_F (o);
        thread = O_F_THREAD (o);

        /* Skip faces on current thread */
        if (THREAD_ID (thread) == tid)
        {
            continue;
        }

        /* Note ID of the other thread for posterity */
        if (nid == -1)
        {
            nid = THREAD_ID (thread);
        }

        /* Initialize to zero */
        N3V_S (nctmp, =, 0.0);
        N3V_S (xctmp, =, 0.0);
        F_AREA (nctmp, face, thread);
        F_CENTROID (xctmp, face, thread);
        /**
        * Negative sum because wall normals
        * point out of the fluid domain
        */
        N3V_V (nc, -=, nctmp);
        N3V_V (xc, +=, xctmp);
        n_faces++;
    }
    # if RP_NODE
    {
        /* Reduce in parallel */
        nid = PRF_GIHIGH1 (nid);
        n_faces = PRF_GISUM1 (n_faces);
        PRF_GRSUM3 (nc[0], nc[1], nc[2]);
        PRF_GRSUM3 (xc[0], xc[1], xc[2]);
    }
    # endif

    /* Propagate to host */
    node_to_host_int_2 (nid, n_faces);
    node_to_host_real_3 (nc[0], nc[1], nc[2]);
    node_to_host_real_3 (xc[0], xc[1], xc[2]);
    if (n_faces > 0)
    {
        nc_mag = N3V_MAG (nc) + REAL_MIN;
        N3V_S (nc, /=, nc_mag);
        N3V_S (xc, /=, n_faces);
    }
    else
    {
        return;
    }

    Message
    (
    "\nContact:: tid: %d nid: %d n_faces: %d "
    "Point: (%f %f %f) Normal: (%f %f %f)",
    tid, nid, n_faces,
    xc[0], xc[1], xc[2],
    nc[0], nc[1], nc[2]
    );

    /* Fetch thread for opposite body */
    domain = THREAD_DOMAIN (DT_THREAD (dt));
    thread = Lookup_Thread (domain, nid);
    if (NULLP (thread))
    {
        Message ("\nWarning: No thread for nid: %d ", nid);
        return;
    }
    else
    {
        ndt = THREAD_DT (thread);
    }

    /* Fetch body parameters */
    SDOF_Get_Motion (dt, vel0, omega0, theta0);

    /* Compute difference vectors and normalize */
    N3V_VV (norm0, =, xc, -, DT_CG (dt));
    norm0_mag = N3V_MAG (norm0) + REAL_MIN;
    N3V_S (norm0, /=, norm0_mag);
    if (NULLP (ndt))
    {
        /* Stationary body / wall. Use contact normal */
        N3V_V (norm1, =, nc);
        /* Compute relative velocity */
        N3V_S (vel1, =, 0.0);
        N3V_V (vel_rel, =, vel0);
    }
    else
    {
        /* Fetch body parameters */
        SDOF_Get_Motion (ndt, vel1, omega1, theta1);
        /* Compute relative velocity */
        N3V_VV (vel_rel, =, vel0, -, vel1);
        /* Compute difference vectors and normalize */
        N3V_VV (norm1, =, xc, -, DT_CG (ndt));
        norm1_mag = N3V_MAG (norm1) + REAL_MIN;
        N3V_S (norm1, /=, norm1_mag);

        /* Check if velocity needs to be reversed */
        if (N3V_DOT (vel_rel, nc) < 0.0)
        {
            /* Reflect velocity across the normal */
            v1dotn0 = 2.0 * N3V_DOT (vel1, norm0);
            N3V_S (norm0, *=, v1dotn0);
            N3V_V (vel1, -=, norm0);
            /* Override body velocity */
            SDOF_Overwrite_Motion (ndt, vel1, omega1, theta1);
        }
    }

    /* Check if velocity needs to be reversed */
    if (N3V_DOT (vel_rel, nc) < 0.0)
    {
        /* Reflect velocity across the normal */
        v0dotn1 = 2.0 * N3V_DOT (vel0, norm1);
        N3V_S (norm1, *=, v0dotn1);
        N3V_V (vel0, -=, norm1);
        /* Override body velocity */
        SDOF_Overwrite_Motion (dt, vel0, omega0, theta0);
    }

    Message
    (
    "\ncontact_props: Updated :: vel0 = (%f %f %f) vel1 = (%f %f %f)",
    vel0[0], vel0[1], vel0[2], vel1[0], vel1[1], vel1[2]
    );
}