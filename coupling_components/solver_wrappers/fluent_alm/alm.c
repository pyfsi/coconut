#include "udf.h"
#include "cxndsearch.h"
#include <math.h>
#include <stdio.h>

/* Interpolate density, yarn orientation, dynamic viscosity and velocity components */
#define SCALAR_R 0
#define SCALAR_THETA 1
#define SCALAR_MU 2
#define SCALAR_U 3
#define SCALAR_V 4
#if RP_3D
#define SCALAR_W 5
#endif /* RP_3D */
#define NSCALARS (SCALAR_U+ND_ND)

/* User-defined memory */
#define UDMI_FMAG 0
#define UDMI_SX 1
#define UDMI_SY 2
#if RP_3D
#define UDMI_SZ 3
#endif /* RP_3D */

/* Declare data fitting coefficients of the force coefficients */
real c_axial[5];
real c_drag[4];
real c_moment[6];

/* Parameters */
#define NYARNPOINTS |NYARNPOINTS|
#define YARN_DIAMETER |YARN_DIAMETER|
#define G_EPS |G_EPS|
#define N_CIRC_S |N_CIRC_S|
#define DELTA_T |DELTA_T|
#define UNSTEADY |UNSTEADY|

typedef struct yarn_point_struct {
	int partition;
	int highest_partition;
	cell_t cell;
	Thread *cell_thread;
	cxboolean found;
	real vs_w[N_CIRC_S+1];  // velocity sampling: sum of weights
} Yarn_Point;

/* Global variables */
/* These global variables are not allocated in a dynamic manner, as it is hard to do so. 
Moreover, they are needed during the entire timestep, so there is not much gain in allocating them dynamically.
At the moment, only 1 yarn is defined, but by adding an extra dimension to these arrays, this can be changed rather easily. */
Yarn_Point yarn_points[NYARNPOINTS];
real yarn_coordinates[NYARNPOINTS][ND_ND];
real yarn_coordinates_old[NYARNPOINTS][ND_ND];
real yarn_tangent[NYARNPOINTS][ND_ND];
real yarn_velocities[NYARNPOINTS][ND_ND];
real yarn_forces[NYARNPOINTS][ND_ND];

real local_theta[NYARNPOINTS];
real sample_coordinates[NYARNPOINTS][N_CIRC_S+1][ND_ND];
real sampled_velocity[NYARNPOINTS][2][ND_ND];
real air_values[NYARNPOINTS][NSCALARS];

int timestep = 0;


  /*------------------*/
 /* helper_functions */
/*------------------*/


void read_yarn_points(char file_name[256])
{
#if RP_HOST
    FILE *input;
    int n_points, point, d;

    if (NULLP(input = fopen(file_name, "r"))) {
        Error("\nUDF-error: Unable to open %s for reading\n", file_name);
        exit(1);
    }

    fscanf(input, "%i", &n_points);
    if (NYARNPOINTS != n_points)
    {
        Error("UDF-error: Number of nodes (%i) not as expected (%i)\n", n_points, NYARNPOINTS);
        exit(1);
    }

    for (point = 0; point < NYARNPOINTS; point++)
    {
        for (d=0; d<2*ND_ND; d++)
        {
            if (d<ND_ND) {
                fscanf(input, "%lf", &yarn_coordinates[point][d]);
            }
            else {
                fscanf(input, "%lf", &yarn_tangent[point][d-ND_ND]);
            }
        }
    }
    fclose(input);

#endif /* RP_HOST */

    host_to_node_real(&yarn_coordinates[0][0], ND_ND*NYARNPOINTS);
    host_to_node_real(&yarn_tangent[0][0], ND_ND*NYARNPOINTS);
}

void find_cell_at_yarn_points(Yarn_Point point_struct[NYARNPOINTS], real coordinate_list[NYARNPOINTS][ND_ND])
{
#if RP_NODE
    CX_Cell_Id *cx_cell;
    static ND_Search *domain_table = NULL;
    int point, nwithforce;
    cell_t cell;
	Thread *cell_thread;

    /* initialize all yarn points as not found */
	for (point = 0; point < NYARNPOINTS; point++)
	{
		point_struct[point].found = FALSE;
	}

    /* search for points */
    domain_table = CX_Start_ND_Point_Search(domain_table, TRUE, -1);
    for (point = 0; point < NYARNPOINTS; point++)
    {
        cx_cell = CX_Find_Cell_With_Point(domain_table, coordinate_list[point], 0.0);
        if (cx_cell)
        {
            cell = RP_CELL(cx_cell);
            cell_thread = RP_THREAD(cx_cell);
            point_struct[point].partition = C_PART(cell, cell_thread);  // Partition id of process that has 'cell' as internal cell
            point_struct[point].cell = cell;
            point_struct[point].cell_thread = cell_thread;
            point_struct[point].found = TRUE;
        }
    }
    domain_table = CX_End_ND_Point_Search(domain_table);

    nwithforce = 0;
    for (point = 0; point < NYARNPOINTS; point ++)
    {
        point_struct[point].highest_partition = PRF_GIHIGH1(point_struct[point].partition);
        point_struct[point].found = PRF_GIHIGH1(point_struct[point].found);
        if (point_struct[point].found && (myid == point_struct[point].highest_partition))
        {
            nwithforce++;
        }
    }

    nwithforce = PRF_GISUM1(nwithforce);
    /*if (myid == 0) {printf("\n%i of %i points found in fluid domain\n", nwithforce, NYARNPOINTS); fflush(stdout);}*/

#endif /* RP_NODE */
}

void set_yarn_velocities()
{
#if RP_NODE
    int point;
    if ((UNSTEADY) && (timestep > 0))
    {
        for (point = 0; point < NYARNPOINTS; point++)
        {
            NV_VS_VS(yarn_velocities[point], =, yarn_coordinates[point], /, DELTA_T, -, yarn_coordinates_old[point], /, DELTA_T);
        }
    }
    else
    {
        for (point = 0; point < NYARNPOINTS; point++)
        {
            NV_S(yarn_velocities[point], =, 0.0);
            NV_V(yarn_coordinates_old[point], =, yarn_coordinates[point]);
        }
    }
#endif /* RP_NODE */
}

void calculate_air_velocity_density_yarn_orientation()
{
#if RP_NODE
    int point, n, d;
    real dot_prod, norm_tan, norm_vel, r, g, omega, mag;
    real e_v[ND_ND], e_drag[ND_ND], centroid[ND_ND], dist[ND_ND];
    cell_t cell;
    Thread *cell_thread;
    Domain *domain;

    real work[NYARNPOINTS*2*ND_ND];
    for (n = 0; n < NYARNPOINTS*2*ND_ND; n++)
    {
        work[n] = 0.0;
    }
    real iwork[N_CIRC_S+1];
    for (n = 0; n < N_CIRC_S+1; n ++)
    {
        iwork[n] = 0.0;
    }
    real iiwork[NYARNPOINTS*NSCALARS];
    for (n = 0; n < NYARNPOINTS*NSCALARS; n++)
    {
        iiwork[n] = 0.0;
    }

    /* Initialize values */
    for (point = 0; point < NYARNPOINTS; point++)
    {
        for (n = 0; n < NSCALARS; n++)
        {
            air_values[point][n] = 0.0;
        }
        for (n = 0; n < N_CIRC_S+1; n ++)
        {
            yarn_points[point].vs_w[n] = 0.0;
            NV_S(sample_coordinates[point][n], =, 0.0);
        }
        NV_S(sampled_velocity[point][0], =, 0.0);
        NV_S(sampled_velocity[point][1], =, 0.0);
    }

    /* Integral velocity sampling at actuator point, no normalization needed  */
    domain = Get_Domain(1);
    thread_loop_c(cell_thread, domain)
	{
		begin_c_loop_int(cell, cell_thread) // Only loop over interior cells to avoid double counting!
		{
		    C_CENTROID(centroid, cell, cell_thread);
		    for(point = 0; point < NYARNPOINTS; point++)
		    {
		        if (yarn_points[point].found)
		        {
		            NV_VV(dist, =, yarn_coordinates[point], -, centroid);
                    g = pow(G_EPS, -3.0) * pow(M_PI, -1.5) * exp(-NV_MAG2(dist)/pow(G_EPS, 2.0));

		            sampled_velocity[point][0][0] += g * C_VOLUME(cell, cell_thread) * C_U(cell, cell_thread);
		            sampled_velocity[point][0][1] += g * C_VOLUME(cell, cell_thread) * C_V(cell, cell_thread);
#if RP_3D
		            sampled_velocity[point][0][2] += g * C_VOLUME(cell, cell_thread) * C_W(cell, cell_thread);
#endif /* RP_3D */
		        }
		    }
		} end_c_loop_int(cell, cell_thread)
    }

    /* Synchronize velocities over different processors, as this is needed to compute the new sampling point */
    PRF_GRSUM(&sampled_velocity[0][0][0], NYARNPOINTS*2*ND_ND, work);

    /* Get sampling point from velocity at actuator point*/
    for (point = 0; point < NYARNPOINTS; point++)
    {
        NV_VV(e_v, =, sampled_velocity[point][0], -, yarn_velocities[point]);  // Compute relative velocity
        NV_S(sampled_velocity[point][0], =, 0.0);   // Clear stored velocity
        mag = NV_MAG(e_v);
        if (mag > DBL_EPSILON)
        {
            NV_S(e_v, /=, mag);
        }
        else
        {
            NV_S(e_v, *=, 0.0);
        }
        mag = NV_MAG(yarn_tangent[point]);
        NV_S(yarn_tangent[point], /=, mag);
//        r = G_EPS * sqrt(-log(0.01));  // Sample at radius where force influence reaches 1% of peak value
        r = YARN_DIAMETER;  // Sample at radius equal to yarn diameter
        local_theta[point] = acos(NV_DOT(e_v, yarn_tangent[point]));

        /* Cross-flow component */
        NV_V_VS(e_drag, =, e_v, -, yarn_tangent[point], *, NV_DOT(e_v, yarn_tangent[point]));
        mag = NV_MAG(e_drag);
        if (mag > DBL_EPSILON)
        {
            NV_S(e_drag, /=, mag);
        }
        else
        {
            NV_S(e_drag, *=, 0.0);
        }
        NV_V_VS(sample_coordinates[point][0], =, yarn_coordinates[point], -, e_drag, *, r);

        /* Axial flow component, only implemented for 3D! */
        real base_vec[3], u[3], e_z[3];
        real alpha_inc, beta;

        NV_S(e_z, =, 0.0);
        e_z[2] = 1.0;
        alpha_inc = 2*M_PI/N_CIRC_S;

        /* Construct rotation vector */
        NV_CROSS(u, e_z, yarn_tangent[point]);
        mag = NV_MAG(u);
        beta = asin(mag);
        if (mag > DBL_EPSILON)
        {
            NV_S(u, /=, mag);
        }
        else
        {
            NV_S(u, *=, 0.0);
        }

        for (n = 0; n < N_CIRC_S; n++)
        {
            /* Construct circle points in xy-plane around origin */
            base_vec[0] = r*cos(alpha_inc*n);
            base_vec[1] = r*sin(alpha_inc*n);
            base_vec[2] = 0.0;

            /* Construct rotated circle points */
            if (abs(beta) > DBL_EPSILON)
            {
                sample_coordinates[point][n+1][0] = yarn_coordinates[point][0] + (cos(beta)+u[0]*u[0]*(1-cos(beta)))*base_vec[0] + (u[0]*u[1]*(1-cos(beta))-u[2]*sin(beta))*base_vec[1];
                sample_coordinates[point][n+1][1] = yarn_coordinates[point][1] + (u[1]*u[0]*(1-cos(beta))+u[2]*sin(beta))*base_vec[0] + (cos(beta)+u[1]*u[1]*(1-cos(beta)))*base_vec[1];
                sample_coordinates[point][n+1][2] = yarn_coordinates[point][2] + (u[2]*u[0]*(1-cos(beta))-u[1]*sin(beta))*base_vec[0] + (u[2]*u[1]*(1-cos(beta))+u[0]*sin(beta))*base_vec[1];
            }
            else  // yarn_tangent is already aligned with e_z: no rotation possible/needed
            {
                NV_VV(sample_coordinates[point][n+1], =, yarn_coordinates[point], +, base_vec);
            }
        }
    }

    /* Compute sum of interpolation weights */
    thread_loop_c(cell_thread, domain)
    {
        begin_c_loop_int(cell, cell_thread)
        {
            C_CENTROID(centroid, cell, cell_thread);
            for (point = 0; point < NYARNPOINTS; point ++)
            {
                if (yarn_points[point].found)
		        {
		            for (n = 0; n < N_CIRC_S+1; n ++)
                    {
                        NV_VV(dist, =, sample_coordinates[point][n], -, centroid);
                        g = pow(G_EPS, -3.0) * pow(M_PI, -1.5) * exp(-NV_MAG2(dist)/pow(G_EPS, 2.0));
                        yarn_points[point].vs_w[n] += g * C_VOLUME(cell, cell_thread);
		            }
                }
            }
        } end_c_loop_int(cell, cell_thread)
    }

    /* Synchronize sum of weights across processors */
    for (point = 0; point < NYARNPOINTS; point ++)
    {
        PRF_GRSUM(&yarn_points[point].vs_w[0], N_CIRC_S+1, iwork);
    }

    /* Integral velocity sampling at new sample points */
    thread_loop_c(cell_thread, domain)
	{
		begin_c_loop_int(cell, cell_thread) // Only loop over interior cells to avoid double counting!
		{
		    C_CENTROID(centroid, cell, cell_thread);
		    for(point = 0; point < NYARNPOINTS; point++)
		    {
		        if (yarn_points[point].found)
		        {
		            /* Cross-flow sampling */
		            NV_VV(dist, =, sample_coordinates[point][0], -, centroid);
                    g = pow(G_EPS, -3.0) * pow(M_PI, -1.5) * exp(-NV_MAG2(dist)/pow(G_EPS, 2.0));

		            sampled_velocity[point][0][0] += g * C_VOLUME(cell, cell_thread) / yarn_points[point].vs_w[0] * C_U(cell, cell_thread);
		            sampled_velocity[point][0][1] += g * C_VOLUME(cell, cell_thread) / yarn_points[point].vs_w[0] * C_V(cell, cell_thread);
#if RP_3D
		            sampled_velocity[point][0][2] += g * C_VOLUME(cell, cell_thread) / yarn_points[point].vs_w[0] * C_W(cell, cell_thread);
#endif /* RP_3D */

                    /* Axial flow sampling */
                    for (n = 0; n < N_CIRC_S; n ++)
                    {
                        NV_VV(dist, =, sample_coordinates[point][n], -, centroid);
                        g = pow(G_EPS, -3.0) * pow(M_PI, -1.5) * exp(-NV_MAG2(dist)/pow(G_EPS, 2.0));

                        sampled_velocity[point][1][0] += g * C_VOLUME(cell, cell_thread) / yarn_points[point].vs_w[n] * C_U(cell, cell_thread) / N_CIRC_S;
                        sampled_velocity[point][1][1] += g * C_VOLUME(cell, cell_thread) / yarn_points[point].vs_w[n] * C_V(cell, cell_thread) / N_CIRC_S;
#if RP_3D
                        sampled_velocity[point][1][2] += g * C_VOLUME(cell, cell_thread) / yarn_points[point].vs_w[n] * C_W(cell, cell_thread) / N_CIRC_S;
#endif /* RP_3D */
                    }
		        }
		    }
		} end_c_loop_int(cell, cell_thread)
    }

    /* Synchronize velocities over different processors, as this is needed to compute the correct flow angle */
    PRF_GRSUM(&sampled_velocity[0][0][0], NYARNPOINTS*2*ND_ND, work);

    /* Now get the other flow variables (density, viscosity) and flow orientation */
    for (point = 0; point < NYARNPOINTS; point++)
    {
        if (yarn_points[point].found && (myid == yarn_points[point].highest_partition))
        {
            cell = yarn_points[point].cell;
            cell_thread = yarn_points[point].cell_thread;
            air_values[point][SCALAR_R] = C_R(cell, cell_thread);
            air_values[point][SCALAR_MU] = C_MU_L(cell, cell_thread);

            /* Use blending function to distinguish between cross-flow sampling and axial flow sampling and make velocity relative to yarn velocity*/
            omega = pow(cos(local_theta[point]), 10);
            for (d = 0; d < ND_ND; d ++)
            {
                air_values[point][SCALAR_U+d] = (1 - omega) * sampled_velocity[point][0][d] + omega * sampled_velocity[point][1][d] - yarn_velocities[point][d];
            }

            /* Compute the orientation of the thread with respect to the flow */
            /* For now: yarn is stationary, so don't take into account yarn velocity yet */
            dot_prod = 0.0;
            norm_tan = 0.0;
            norm_vel = 0.0;
            for (n = 0; n < ND_ND; n++)
            {
                dot_prod += yarn_tangent[point][n] * air_values[point][SCALAR_U+n];
                norm_tan += pow(yarn_tangent[point][n], 2);
                norm_vel += pow(air_values[point][SCALAR_U+n], 2);
            }
            if ((sqrt(norm_tan) < DBL_EPSILON) || (sqrt(norm_vel) < DBL_EPSILON))
            {
                dot_prod = 0.0;
            }
            else
            {
                dot_prod = dot_prod / (sqrt(norm_tan) * sqrt(norm_vel));
            }
            air_values[point][SCALAR_THETA] = acos(dot_prod);
        }
    }

    /* Synchronize flow variables */
    PRF_GRSUM(&air_values[0][0], NYARNPOINTS*NSCALARS, iiwork);

    /*if (myid == 0)
    {
        for (point = 0; point < NYARNPOINTS; point++)
        {
            printf("Point %i: v = [%5.2f, %5.2f, %5.2f] m/s; v_yarn = [%5.2f, %5.2f, %5.2f] m/s; rho = %5.3f kg/m^3; mu = %5.3e N.s/m^2; theta = %5.3f°;\n", point, air_values[point][SCALAR_U], air_values[point][SCALAR_V], air_values[point][SCALAR_W], NV_LIST(yarn_velocities[point]), air_values[point][SCALAR_R], air_values[point][SCALAR_MU], air_values[point][SCALAR_THETA]*180.0/M_PI);
            fflush(stdout);
        }
    }*/
#endif /* RP_NODE */
}

void initialize_force_coefficients()
{
    c_axial[0] = 0.900;
    c_axial[1] = 0.316;
    c_axial[2] = 1.355;
    c_axial[3] = 0.078;
    c_axial[4] = 1.120;

    c_drag[0] = 0.811;
    c_drag[1] = 1.395;
    c_drag[2] = 10.593;
    c_drag[3] = 0.847;

    c_moment[0] = 6.554;
    c_moment[1] = 1.012;
    c_moment[2] = -0.008;
    c_moment[3] = 1.863;
    c_moment[4] = -0.014;
    c_moment[5] = 0.815;
}

void calculate_yarn_forces()
{
#if RP_NODE
    int point, d;
    real e_drag[ND_ND], e_v[ND_ND], velocity[ND_ND], f_a[ND_ND], f_drag[ND_ND];
    real Re, p_dyn, c_a, c_td, multiplier, mag;

    initialize_force_coefficients();

    for (point = 0; point < NYARNPOINTS; point++){

        /* Initialize with zero force */
        NV_S(yarn_forces[point], =, 0.0);

        if (yarn_points[point].found)  // Only calculate yarn force if yarn point is actually in domain
        {
            /* Calculate local reference frame */
            for (d = 0; d < ND_ND; d++)
            {
                velocity[d] = air_values[point][SCALAR_U+d];
            }
            /* Avoid NaN values at zero velocity */
            mag = NV_MAG(velocity);
            if (mag > DBL_EPSILON)
            {
                NV_V(e_v, =, velocity);
                NV_S(e_v, /=, mag);
                NV_V_VS(e_drag, =, e_v, -, yarn_tangent[point], *, NV_DOT(e_v, yarn_tangent[point]));
                mag = NV_MAG(e_drag);
                if (mag < DBL_EPSILON)
                {
                    NV_S(e_drag, *=, 0.);
                }
                else
                {
                    NV_S(e_drag, /=, mag);
                }

                /* Check orientation of flow w.r.t. yarn tangent vector */
                if (NV_DOT(e_v, yarn_tangent[point]) >= 0.)
                {
                    multiplier = 1.;
                }
                else
                {
                    multiplier = -1.;
                }

                /* Calculate flow values */
                Re = air_values[point][SCALAR_R] * NV_MAG(velocity) * YARN_DIAMETER / air_values[point][SCALAR_MU];
                p_dyn = 0.5 * air_values[point][SCALAR_R] * pow(NV_MAG(velocity), 2) * YARN_DIAMETER;

                /* Calculate axial force on yarn */
                c_a = c_axial[0] * pow(Re, -c_axial[1]) * pow(fabs(cos(air_values[point][SCALAR_THETA])), c_axial[2]) + c_axial[3] * pow(fabs(sin(2*air_values[point][SCALAR_THETA])), c_axial[4]);
                NV_V(f_a, =, yarn_tangent[point]);
                NV_S(f_a, *=, (c_a * p_dyn * multiplier));

                /* Calculate transversal drag force on yarn */
                c_td = c_drag[0] * pow(fabs(sin(air_values[point][SCALAR_THETA])), c_drag[1]) + c_drag[2] * pow(Re, -c_drag[3]) * pow(fabs(sin(air_values[point][SCALAR_THETA])), (c_drag[1]-0.5));
                NV_V(f_drag, =, e_drag);
                NV_S(f_drag, *=, (c_td * p_dyn));

                /* Store forces acting on the yarn in a global variable */
                NV_VV(yarn_forces[point], =, f_a, +, f_drag);
            }
            else
            {
                NV_S(e_v, =, 0.0);
                NV_S(e_drag, =, 0.0);
                NV_S(f_a, =, 0.0);
                NV_S(f_drag, =, 0.0);
                c_a = 0.0;
                c_td = 0.0;
            }
        }

        /* Print to console */
        /*if (myid == 0){
            printf("\nPoint %i: c_a = %10.5f; c_t,d = %10.5f\n", point, c_a, c_td);
            printf("Coordinates: [%10.5f, %10.5f, %10.5f]\n", NV_LIST(yarn_coordinates[point]));
            printf("Velocity: [%10.5f, %10.5f, %10.5f], magnitude: %10.5f\n", NV_LIST(velocity), NV_MAG(velocity));
            printf("Reynolds number: %10.5f, dynamic pressure: %10.5f, theta: %10.5f\n", Re, p_dyn, air_values[point][SCALAR_THETA]);
            printf("e_a: [%10.5f, %10.5f, %10.5f], magnitude: %10.5f\n", NV_LIST(yarn_tangent[point]), NV_MAG(yarn_tangent[point]));
            printf("e_td: [%10.5f, %10.5f, %10.5f], magnitude: %10.5f\n", NV_LIST(e_drag), NV_MAG(e_drag));
            printf("f_drag: %10.5f; f_a : %10.5f\n", NV_MAG(f_drag), NV_MAG(f_a));
            printf("Force on yarn: [%10.5f, %10.5f, %10.5f] N/m\n", NV_LIST(yarn_forces[point]));
            fflush(stdout);
        }*/
    }
#endif /* RP_NODE */
}


  /*----------------*/
 /* store_traction */
/*----------------*/


DEFINE_ON_DEMAND(store_traction)
{
	/* Send coordinates from compute node 0 */
	node_to_host_real(&yarn_forces[0][0], NYARNPOINTS*ND_ND);
	node_to_host_real(&air_values[0][0], NYARNPOINTS*NSCALARS);

	/* Write coordinates on host */
#if RP_HOST
	int point, d;
    timestep = RP_Get_Integer("udf/timestep");
	char output_file_name[256];
	FILE *output;

	sprintf(output_file_name, "traction_timestep%i.dat", timestep);
	if (NULLP(output = fopen(output_file_name, "w"))) {
        Error("\nUDF-error: Unable to open %s for writing\n", output_file_name);
        exit(1);
    }
    fprintf(output, "%27s %27s %27s  %10s\n","x-force [N/m]", "y-force [N/m]", "z-force [N/m]", "unique-ids");
	for (point = 0; point < NYARNPOINTS; point++){
	    for (d = 0; d < ND_ND; d++) {
	        fprintf(output, "%27.17e ", yarn_forces[point][d]);
	    }
		fprintf(output, " %10d\n", point);
	}
	fclose(output);

	sprintf(output_file_name, "air_values_timestep%i.dat", timestep);
	if (NULLP(output = fopen(output_file_name, "w"))) {
        Error("\nUDF-error: Unable to open %s for writing\n", output_file_name);
        exit(1);
    }
    fprintf(output, "%27s %27s %27s  %10s\n","x-velocity [m/s]", "y-velocity [m/s]", "z-velocity [m/s]", "unique-ids");
	for (point = 0; point < NYARNPOINTS; point++)
	{
	    for (d = 0; d < ND_ND; d++)
	    {
		    fprintf(output, "%27.17e ", air_values[point][SCALAR_U+d]);
	    }
		fprintf(output, " %10d\n", point);
	}
	fclose(output);
#endif /* RP_HOST */

    if (myid == 0) {printf("\nFinished UDF store_traction\n"); fflush(stdout);}

}


  /*-------------*/
 /* set_sources */
/*-------------*/


DEFINE_ADJUST(set_source, domain)
{
#if RP_NODE
    Thread *cell_thread;
    cell_t cell;
    int point;
    real centroid[ND_ND], r[ND_ND], r_n[ND_ND], t[ND_ND], yt[ND_ND], yt_n[ND_ND], dist[ND_ND], force[ND_ND];
    real g, l_elem, p_s, p_n;

    /* Calculate new yarn forces */
    calculate_air_velocity_density_yarn_orientation();
    calculate_yarn_forces();

    /* Calculate momentum sources and store in UDM's */
    thread_loop_c(cell_thread, domain)
    {
        begin_c_loop_int(cell, cell_thread)
        {
            /* Initialize all UDMI's at zero value */
            C_UDMI(cell, cell_thread, UDMI_FMAG) = 0.0;
            C_UDMI(cell, cell_thread, UDMI_SX) = 0.0;
            C_UDMI(cell, cell_thread, UDMI_SY) = 0.0;
#if RP_3D
            C_UDMI(cell, cell_thread, UDMI_SZ) = 0.0;
#endif /* RP_3D */

            C_CENTROID(centroid, cell, cell_thread);

            /* For each yarn point (actuator element), calculate distance to current cell and apply Gaussian smoothing */
            for (point = 0; point < NYARNPOINTS-1; point++)
            {
                /* Calculate tangential coordinate of current actuator element and element length */
                NV_VV(t, =, yarn_coordinates[point+1], -, yarn_coordinates[point]);
                NV_VV(r, =, centroid, -, yarn_coordinates[point]);
                NV_VV(r_n, =, centroid, -, yarn_coordinates[point+1]);
                l_elem = NV_MAG(t);
                p_s = NV_DOT(r, t) / pow(l_elem, 2);

                /* Get local tangents to yarn points and align with actuator element orientation */
                NV_V(yt, =, yarn_tangent[point]);
                if (NV_DOT(yt, t) < 0.0)
                {
                    NV_S(yt, *=, -1);
                }
                NV_V(yt_n, =, yarn_tangent[point+1]);
                if (NV_DOT(yt_n, t) < 0.0)
                {
                    NV_S(yt_n, *=, -1);
                }

                /* Cell is influenced by current actuator element */
                if ((NV_DOT(yt, r) >= 0.0) && (NV_DOT(yt_n, r_n) < 0.0))
                {
                    if ((p_s >= 0.0) && (p_s <= 1.0))  // Concave case automatically covered due to outer if-loop
                    {
                        /* Calculate radial distance to current actuator element */
                        NV_CROSS(dist, r, t);
                        p_n = NV_MAG(dist) / l_elem;
                    }
                    /* Treat convex angles */
                    else if (p_s < 0.0)
                    {
                        p_s = 0.0;
                        p_n = NV_MAG(r);
                    }
                    else if (p_s > 1.0)
                    {
                        p_s = 1.0;
                        p_n = NV_MAG(r_n);
                    }

                    /* Interpolate forces, calculate 2D Gaussian spreading factor, apply to cell UDMI's */
                    NV_VS_VS(force, =, yarn_forces[point], *, (1-p_s), +, yarn_forces[point+1], *, p_s);
                    g = 1 / (G_EPS * G_EPS * M_PI) * exp(-(p_n*p_n)/(G_EPS*G_EPS));

                    C_UDMI(cell, cell_thread, UDMI_FMAG) += g * NV_MAG(force);
                    C_UDMI(cell, cell_thread, UDMI_SX) -= g * force[0];
                    C_UDMI(cell, cell_thread, UDMI_SY) -= g * force[1];
#if RP_3D
                    C_UDMI(cell, cell_thread, UDMI_SZ) -= g * force[2];
#endif /* RP_3D */
                }
            }
        } end_c_loop_int(cell, cell_thread)
    }
#endif /* RP_NODE */
}

DEFINE_SOURCE(xmom_source, cell, cell_thread, dS, eqn)
{
	real source;
	source = C_UDMI(cell, cell_thread, UDMI_SX);
	dS[eqn] = 0.0;

	return source;
}

DEFINE_SOURCE(ymom_source, cell, cell_thread, dS, eqn)
{
	real source;
	source = C_UDMI(cell, cell_thread, UDMI_SY);
	dS[eqn] = 0.0;

	return source;
}

#if RP_3D
DEFINE_SOURCE(zmom_source, cell, cell_thread, dS, eqn)
{
	real source;
	source = C_UDMI(cell, cell_thread, UDMI_SZ);
	dS[eqn] = 0.0;

	return source;
}
#endif /* RP_3D */


  /*------------------*/
 /* read_coordinates */
/*------------------*/


DEFINE_ON_DEMAND(read_coordinates)
{
#if RP_HOST
    timestep = RP_Get_Integer("udf/timestep");
#endif /* RP_HOST */
    host_to_node_int_1(timestep);

    /* Read yarn coordinates, find in domain, set yarn velocities */
    char filename[256];
    sprintf(filename, "coordinates_update_timestep%i.dat", timestep);
    read_yarn_points(filename);
    find_cell_at_yarn_points(yarn_points, yarn_coordinates);
    set_yarn_velocities();
  
//#if RP_NODE
    /* Calculate cell size at actuator points */
/*    for (point = 0; point < NYARNPOINTS; point ++)
    {
        if ((yarn_points[point].found) && (myid == yarn_points[point].highest_partition))
        {
            cell = yarn_points[point].cell;
            cell_thread = yarn_points[point].cell_thread;
            C_CENTROID(centroid, cell, cell_thread);
            printf("Point %i; coordinates [%10.5e, %10.5e, %10.5e]; cell centroid [%10.5e, %10.5e, %10.5e]; cell size %10.5e; epsilon/dh %10.5e \n", point, NV_LIST(yarn_coordinates[point]), NV_LIST(centroid), pow(C_VOLUME(cell, cell_thread), 1./3.), G_EPS/pow(C_VOLUME(cell, cell_thread), 1./3.));
            fflush(stdout);
        }
    }*/
//#endif /* RP_NODE */

    if (myid == 0) {printf("\nFinished UDF read_coordinates\n"); fflush(stdout);}
}


  /*---------------*/
 /* increase_time */
/*---------------*/


DEFINE_ON_DEMAND(increase_time)
{
/* Update yarn coordinates from the previous timestep before progressing to the new timestep */
#if RP_NODE
    int point;
    for (point = 0; point < NYARNPOINTS; point++)
    {
        NV_V(yarn_coordinates_old[point], =, yarn_coordinates[point]);
    }
#endif /* RP_NODE */

    /* Send coordinates from compute node 0 */
	node_to_host_real(&yarn_coordinates_old[0][0], NYARNPOINTS*ND_ND);

    if (myid == 0) {printf("\nFinished UDF increase_time\n"); fflush(stdout);}
}


  /*----------------------*/
 /* store_coordinates_id */
/*----------------------*/

DEFINE_ON_DEMAND(store_coordinates_id)
{
	/* Write coordinates and node id's on host */
#if RP_HOST
	int point, d;
    timestep = RP_Get_Integer("udf/timestep");
	char output_file_name[256];
	FILE *output;

	sprintf(output_file_name, "nodes_timestep%i.dat", timestep);

	if (NULLP(output = fopen(output_file_name, "w"))) {
        Error("\nUDF-error: Unable to open %s for writing\n", output_file_name);
        exit(1);
    }

#if RP_2D
	fprintf(output,  "%27s %27s  %10s\n", "x-coordinate", "y-coordinate", "unique-ids");
#else /* RP_2D */
	fprintf(output,  "%27s %27s %27s  %10s\n", "x-coordinate", "y-coordinate", "z-coordinate", "unique-ids");
#endif /* RP_2D */

	for (point = 0; point < NYARNPOINTS; point++)
	{
	    for (d = 0; d < ND_ND; d++)
	    {
                    fprintf(output, "%27.17e ", yarn_coordinates[point][d]);
        }
        fprintf(output, "%10d\n", point);
	}
	fclose(output);
#endif /* RP_HOST */

	if (myid == 0) {printf("\nFinished UDF store_coordinates_id\n"); fflush(stdout);}
}
