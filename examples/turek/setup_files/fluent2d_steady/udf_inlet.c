#include "udf.h"

DEFINE_PROFILE(x_velocity_steady,face_thread,nv) {
#if !RP_HOST
	real y, coordinates[2];
	face_t face;

	begin_f_loop(face,face_thread) {
		F_CENTROID(coordinates,face,face_thread);
		y = coordinates[1];
		F_PROFILE(face,face_thread,nv) = 1.5*0.2*y*(0.41-y)/pow((0.41/2.0), 2.0);
	} end_f_loop(face,face_thread)
#endif /* !RP_HOST */
}
