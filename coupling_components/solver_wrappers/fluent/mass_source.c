# include "udf.h"

/* User defined memory ------------------------------------------------------------------------------------------------
index 0 = previous iteration cell volume                -- timestep m+1 ; iteration i
index 1 = cell volume previous time step converged 	    -- timestep m ; converged solution
-------------------------------------------------------------------------------------------------------------------- */

DEFINE_INIT(init_volumes,d) //use in setup-files
{
Thread *t;
cell_t c;
thread_loop_c(t,d) {
    begin_c_loop(c,t) //loop over all cells
        C_UDMI(c,t,0) = C_VOLUME(c,t);
        C_UDMI(c,t,1) = C_VOLUME(c,t);
    end_c_loop(c,t)
    }
}

DEFINE_ON_DEMAND(update_volumes_it)
{
Domain *d;
d = Get_Domain(1);
Thread *t;
cell_t c;
thread_loop_c(t,d) {
    begin_c_loop(c,t) //loop over all cells
        C_UDMI(c,t,0) = C_VOLUME(c,t);
    end_c_loop(c,t)
    }
}

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
}

/*Source term for energy equation, to include latent heat.*/
DEFINE_SOURCE(udf_mass_source,c,t,dS,eqn)
{
real source = -C_R(c,t)*(C_VOLUME(c,t) - C_UDMI(c,t,1))/(C_VOLUME(c,t)*(CURRENT_TIME-PREVIOUS_TIME));
dS[eqn] = -(C_VOLUME(c,t) - C_UDMI(c,t,1))/(C_VOLUME(c,t)*(CURRENT_TIME-PREVIOUS_TIME));
return source;
}
