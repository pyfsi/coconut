/***********************************************************************
* DEFINE_CONTACT UDF for a Short-Range Repulsive Force Model
*
* This UDF implements the collision model from Glowinski et al. by
* applying a repulsive force when an object comes within close range
* of a wall.
*
* It is assumed that this UDF is called by the solver when the object
* is within a "contact detection" range, *before* geometric overlap.
*
* The calculated force is applied by converting it to a change in
* velocity over the current time step (Impulse-Momentum).
***********************************************************************/

#include "udf.h"
#include "dynamesh_tools.h"

/* --- USER-DEFINED PARAMETERS FOR REPULSIVE FORCE MODEL --- */

/* Define the properties of the moving object and the force model */
#define OBJECT_RADIUS   0.00125   /* Radius of the moving object (spherical/circular) */
#define RHO             0.001   /* Force range, as defined in the paper */
#define EPSILON         1.0e-5 /* "Small" positive stiffness parameter */

/* --- END OF USER-DEFINED PARAMETERS --- */


DEFINE_CONTACT(repulsive_force_contact, dt, contacts)
{
    Objp *o;
    face_t face;
    Thread *thread;
    int tid, nid, n_faces;

    real mass, dt_mag;
    real N3V_VEC(xc), N3V_VEC(nc);
    real N3V_VEC(xctmp), N3V_VEC(nctmp);
    real N3V_VEC(vel0), N3V_VEC(omega0), N3V_VEC(theta0);
    real N3V_VEC(force_vec), N3V_VEC(delta_v);

    if (!Data_Valid_P())
    {
        return;
    }

    /* --- Step 1: Calculate Average "Near-Contact" Point and Normal --- */
    N3V_S(nc, =, 0.0);
    N3V_S(xc, =, 0.0);
    tid = THREAD_ID(DT_THREAD(dt));
    nid = -1;
    n_faces = 0;
    loop(o, contacts)
    {
        face = O_F(o);
        thread = O_F_THREAD(o);
        if (THREAD_ID(thread) == tid) continue;
        if (nid == -1) nid = THREAD_ID(thread);

        /* Initialize to zero */
        N3V_S (nctmp, =, 0.0);
        N3V_S (xctmp, =, 0.0);

        F_AREA(nctmp, face, thread);
        F_CENTROID(xctmp, face, thread);

        N3V_V(nc, -=, nctmp);
        N3V_V(xc, +=, xctmp);
        n_faces++;
    }

#if RP_NODE
    {
        /* Reduce in parallel */
        n_faces = PRF_GISUM1(n_faces);
        PRF_GRSUM3(nc[0], nc[1], nc[2]);
        PRF_GRSUM3(xc[0], xc[1], xc[2]);
    }
#endif

    node_to_host_int_1(n_faces);
    node_to_host_real_3(nc[0], nc[1], nc[2]);
    node_to_host_real_3(xc[0], xc[1], xc[2]);

    if (n_faces > 0)
    {
        dt_mag = N3V_MAG(nc);
        N3V_S(nc, /=, dt_mag); /* Normalize the contact normal */
        N3V_S(xc, /=, n_faces); /* Average the near-contact point */
    }
    else
    {
        return;
    }

    /* --- Step 2: Calculate Repulsive Force based on Separation Distance --- */
    real vec_cg_to_contact[3];
    N3V_VV(vec_cg_to_contact, =, xc, -, DT_CG(dt));
    real dist_cg_to_plane = N3V_DOT(vec_cg_to_contact, nc);
    real surface_dist = dist_cg_to_plane - OBJECT_RADIUS;

    real force_mag = 0.0;

    /* The repulsive force is only active if the surface distance is within the range RHO [cite: 12] */
    if (surface_dist < RHO)
    {
        /* This section implements the force model from Eq. [cite_start]5.2 of the document [cite: 27] */
        /* F_ij = (c_ij / ε) * ( ((d_ij - R_i - R_j - ρ) / ρ)^- )^2 * (unit_vector) */
        /* where ξ^- = max(0, -ξ) [cite: 29] */

        real arg = (surface_dist - RHO) / RHO;
        real penalty_term = pow(fmax(0.0, -arg), 2.0);

        /* A natural choice for the scaling factor is mass * gravity [cite: 31] */
        real force_scaling_factor = DT_MASS(dt) * 9.81;

        force_mag = (force_scaling_factor / EPSILON) * penalty_term;

        /* The force vector acts along the contact normal, pushing away from the wall */
        N3V_S(force_vec, =, nc, *, force_mag);

        /* --- Step 3: Apply Force as an Impulse to Change Velocity --- */
        SDOF_Get_Motion(dt, vel0, omega0, theta0);
        mass = DT_MASS(dt);
        real current_dt = CURRENT_TIMESTEP;

        /* Calculate change in velocity: delta_v = (F * dt) / mass */
        N3V_S(delta_v, =, force_vec, *, (current_dt / mass));

        /* Calculate the new velocity by adding the velocity change */
        N3V_V(vel0, +=, delta_v);

        /* --- Step 4: Override the Object's Motion in the Solver --- */
        SDOF_Overwrite_Motion(dt, vel0, omega0, theta0);

#if RP_HOST
        Message("\nRepulsive Force Active: t=%.4f, Surface Dist: %.6f m, Force: %.4f N, New Vel: (%.3f, %.3f, %.3f)\n",
                CURRENT_TIME, surface_dist, force_mag, vel0[0], vel0[1], vel0[2]);
#endif
    }
}