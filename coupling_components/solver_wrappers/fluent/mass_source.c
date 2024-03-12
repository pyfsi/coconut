# include "udf.h"
#define dt |TIME_STEP_SIZE|

/* User defined memory ------------------------------------------------------------------------------------------------
index 0 = flags the cells adjacent to the interface
index 1 = cell volume previous time step converged
-------------------------------------------------------------------------------------------------------------------- */

DEFINE_ON_DEMAND(update_volumes_ts)
{
Domain *d;
d = Get_Domain(1);
Thread *t;
cell_t c;
thread_loop_c(t,d) {
    begin_c_loop(c,t) //loop over all cells
        C_UDMI(c,t,1) = C_VOLUME(c,t);
    end_c_loop(c,t)
    }
printf("\nFinished UDF update_volumes_ts.\n");
fflush(stdout);
}

/*Source term for energy equation, to include latent heat.*/
DEFINE_SOURCE(udf_mass_source,c,t,dS,eqn)
{
real source = -C_R(c,t)*C_UDMI(c,t,0)*(C_VOLUME(c,t) - C_UDMI(c,t,1))/(C_VOLUME(c,t)*dt);
dS[eqn] = -C_UDMI(c,t,0)*(C_VOLUME(c,t) - C_UDMI(c,t,1))/(C_VOLUME(c,t)*dt);
return source;
}
