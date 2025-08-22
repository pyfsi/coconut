/************************************************************\
* Contact UDF: moving object vs steady, flat wall            *
* - Inelastic normal collision (restitution e_n)             *
* - Tangential Coulomb friction (mu) with stick/slip         *
* - Based on collision handling discussed in the document    *
\************************************************************/

#include "udf.h"
#include "dpm.h"
#include "sg_mphase.h"
#include "mem.h"
#include "dynamesh_tools.h"

/* ====== User parameters (tune as needed) ====== */
#define EN_RESTITUTION   0.50   /* e_n: 0 <= e_n <= 1 (normal) */
#define MU_TANGENTIAL    0.30   /* Coulomb friction coefficient */
#define V_T_STICK        1.0e-6 /* tangential stick threshold (m/s) */

/* ====== Vector helpers come from Fluent's N3V_* macros ====== */

DEFINE_CONTACT(contact_wall_inelastic, dt, contacts)
{
  Objp *o;
  face_t face;
  Thread *thread;
  Domain *domain = NULL;

  int tid, nid, n_faces;

  real nc_mag, norm0_mag;
  real vn, vt_mag, alpha;

  /* Vectors */
  real N3V_VEC (vel_rel);
  real N3V_VEC (nc), N3V_VEC (nctmp);
  real N3V_VEC (xc), N3V_VEC (xctmp);
  real N3V_VEC (vel0), N3V_VEC (omega0), N3V_VEC (theta0);
  real N3V_VEC (vn_vec), N3V_VEC (vt_vec);
  real N3V_VEC (norm0);

  if (!Data_Valid_P())
    return;

  /* Build a common contact plane/point from faces on the other thread */
  N3V_S (nc, =, 0.0);
  N3V_S (xc, =, 0.0);

  /* Current dynamic-thread id (moving body) */
  tid = THREAD_ID (DT_THREAD (dt));
  nid = -1;
  n_faces = 0;

  loop (o, contacts)
  {
    face   = O_F (o);
    thread = O_F_THREAD (o);

    /* Skip faces belonging to this same moving thread */
    if (THREAD_ID (thread) == tid)
      continue;

    if (nid == -1)
      nid = THREAD_ID (thread); /* record “other” thread id (the wall) */

    N3V_S (nctmp, =, 0.0);
    N3V_S (xctmp, =, 0.0);

    F_AREA     (nctmp, face, thread);   /* face area vector = (|A| n)   */
    F_CENTROID (xctmp, face, thread);   /* face centroid                 */

    /* Wall normals point out of the fluid domain => flip sign for contact */
    N3V_V (nc, -=, nctmp);
    N3V_V (xc, +=, xctmp);

    n_faces++;
  }

# if RP_NODE
  {
    /* Parallel reduction */
    nid = PRF_GIHIGH1 (nid);
    n_faces = PRF_GISUM1  (n_faces);
    PRF_GRSUM3 (nc[0], nc[1], nc[2]);
    PRF_GRSUM3 (xc[0], xc[1], xc[2]);
  }
# endif

  /* Propagate to host (required for parallel-safe output) */
  node_to_host_int_2  (nid, n_faces);
  node_to_host_real_3 (nc[0], nc[1], nc[2]);
  node_to_host_real_3 (xc[0], xc[1], xc[2]);

  if (n_faces <= 0)
    return;

  /* Normalize contact normal; average contact point */
  nc_mag = N3V_MAG (nc) + REAL_MIN;
  N3V_S (nc, /=, nc_mag);
  N3V_S (xc, /=, n_faces);

  Message("\nContact(Wall):: tid:%d nid:%d faces:%d   Xc:(%g %g %g)  n:(%g %g %g)",
          tid, nid, n_faces, xc[0], xc[1], xc[2], nc[0], nc[1], nc[2]);

  /* --- Fetch motion for the moving body --- */
  SDOF_Get_Motion (dt, vel0, omega0, theta0);

  /* Build approach direction from body CG to contact, for info */
  N3V_VV (norm0, =, xc, -, DT_CG (dt));
  norm0_mag = N3V_MAG (norm0) + REAL_MIN;
  N3V_S (norm0, /=, norm0_mag);

  /* Wall is steady: relative velocity at contact is simply body velocity */
  N3V_V (vel_rel, =, vel0);

  /* Decompose velocity into normal and tangential components wrt contact plane */
  vn = N3V_DOT (vel_rel, nc);               /* scalar normal component (+ away from wall) */
  N3V_VV (vn_vec, =, nc, *, vn);            /* vector normal component */
  N3V_VV (vt_vec, =, vel_rel, -, vn_vec);   /* tangential component */
  vt_mag = N3V_MAG (vt_vec);

  /* Only handle collisions when approaching the wall: vn < 0 */
  if (vn < 0.0)
  {
    /* ===== Normal, inelastic impact: v_n,out = -e_n * v_n,in ===== */
    {
      real vn_out = -EN_RESTITUTION * vn;
      N3V_VV (vn_vec, =, nc, *, vn_out);
    }

    /* ===== Tangential, Coulomb friction with stick/slip =====
       Use an impulse-like limiter: |Δv_t| ≤ α |v_t|, where
       α = min(1, μ (1+e_n) |v_n,in| / (|v_t| + ε))
       If |v_t| < threshold => sticking (vt_out = 0)
    */
    if (vt_mag <= V_T_STICK)
    {
      N3V_S (vt_vec, =, 0.0); /* sticking */
    }
    else
    {
      alpha = MU_TANGENTIAL * (1.0 + EN_RESTITUTION) * fabs(vn) / (vt_mag + REAL_MIN);
      if (alpha > 1.0) alpha = 1.0;

      /* Reduce tangential speed by factor alpha in the opposite direction */
      /* vt_out = (1 - alpha) * vt_in */
      N3V_S (vt_vec, *=, (1.0 - alpha));
    }

    /* Recompose post-impact velocity (wall is steady) */
    N3V_VV (vel0, =, vn_vec, +, vt_vec);

    /* Overwrite body translational velocity; keep angular as-is */
    SDOF_Overwrite_Motion (dt, vel0, omega0, theta0);

    Message("\nWallContact:: Updated vel = (%g %g %g)  (vn<0 handled)", vel0[0], vel0[1], vel0[2]);
  }
  else
  {
    /* Separating / sliding away: no change */
    /* (Optional) small tangential damping against grazing contact could be added here */
  }
}
