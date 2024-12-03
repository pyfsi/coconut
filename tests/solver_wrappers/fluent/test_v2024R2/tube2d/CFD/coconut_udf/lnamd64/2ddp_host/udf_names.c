/* This file generated automatically. */
/*          Do not modify.            */
#include "udf.h"
#include "prop.h"
#include "dpm.h"
extern DEFINE_ON_DEMAND(get_thread_ids);
extern DEFINE_ON_DEMAND(store_coordinates_id);
extern DEFINE_ON_DEMAND(store_pressure_traction);
extern DEFINE_GRID_MOTION(move_nodes, domain, dynamic_thread, time, dtime);
UDF_Data udf_data[] = {
{"get_thread_ids", (void (*)(void))get_thread_ids, UDF_TYPE_ON_DEMAND},
{"store_coordinates_id", (void (*)(void))store_coordinates_id, UDF_TYPE_ON_DEMAND},
{"store_pressure_traction", (void (*)(void))store_pressure_traction, UDF_TYPE_ON_DEMAND},
{"move_nodes", (void (*)(void))move_nodes, UDF_TYPE_GRID_MOTION},
};
int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data);
#include "version.h"
void UDF_Inquire_Release(int *major, int *minor, int *revision)
{
  *major = RampantReleaseMajor;
  *minor = RampantReleaseMinor;
  *revision = RampantReleaseRevision;
}
