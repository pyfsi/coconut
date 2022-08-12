#include "udf.h"

DEFINE_PROFILE(x_velocity_unsteady,face_thread,nv) {
#if !RP_HOST
    real x[ND_ND]; /* this will hold the position vector */
    real time, velocity, y, h;
    face_t face;
    h = 1.0; /* domain height */
    time = CURRENT_TIME;
    velocity = 1.0-cos(2.0*M_PI*time/5.0);

	begin_f_loop(face,face_thread) {
        F_CENTROID(x, face, face_thread);
        y = x[1];
        F_PROFILE(face,face_thread,nv) = velocity*(y/h-0.875)/(0.125);
    } end_f_loop(face,face_thread)
#endif /* !RP_HOST */
}

DEFINE_GEOM(left,domain,dt,position)
{
    position[0] = 0.0; /* motion along x = 0 */
}

DEFINE_GEOM(right,domain,dt,position)
{
    position[0] = 1.0; /* motion along x = 1 */
}

DEFINE_GRID_MOTION(bottom,domain,dt,time,dtime)
{
    Thread *tf = DT_THREAD(dt);
    face_t f;
    Node *v;
    real NV_VEC(dx);
    real tm;
    int n;
    /* set deforming flag on adjacent cell zone */
    SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));
    tm = (sin (0.4 * 3.1415 * time) + 1) / 2;
    begin_f_loop(f,tf)
    {
        f_node_loop(f,tf,n)
        {
            v = F_NODE(f,tf,n);
            /* update node if the current node has not been previously
            visited when looping through previous faces */
            if (NODE_POS_NEED_UPDATE (v))
            {
                /* indicate that node position has been update
                so that it’s not updated more than once */
                NODE_POS_UPDATED(v);
                dx[0] = 0;
                dx[1] = NODE_Y(v);
                dx[1] -= tm * 0.4 * 4.0 * NODE_X(v) * (NODE_X(v) - 1.0);
                dx[2] = 0;
                NV_V(NODE_COORD(v), +=, dx);
            }
        }
    }
    end_f_loop(f,tf);
}
