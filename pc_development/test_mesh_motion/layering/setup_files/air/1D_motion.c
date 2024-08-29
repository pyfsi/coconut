#include "udf.h"

DEFINE_CG_MOTION(x_mov,dt,vel,omega,time,dtime)
{
    vel[0] = 0.001; // m/s --> 1 mm/s
}