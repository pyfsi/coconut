#include "udf.h"

DEFINE_PROFILE(x_velocity_unsteady,face_thread,nv) {
#if !RP_HOST
	real y, coordinates[2], time, time_factor;
	face_t face;

	time = CURRENT_TIME;
	time_factor = ((time < 2.0) ? (1.0-cos(M_PI/2.0*time))/2.0 : 1.0);
	begin_f_loop(face,face_thread) {
		F_CENTROID(coordinates,face,face_thread);
		y = coordinates[1];
		F_PROFILE(face,face_thread,nv) = 1.5*1.0*y*(0.41-y)/pow((0.41/2.0), 2.0)*time_factor;
	} end_f_loop(face,face_thread)
#endif /* !RP_HOST */
}
