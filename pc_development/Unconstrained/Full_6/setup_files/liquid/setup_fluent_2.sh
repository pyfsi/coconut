#!/bin/bash

# make fluent case with 20 cores for initialisation calculations
fluent 2ddp -t20 -gu -i case_2_flow.jou > setup_fluent.log 2>&1

# delete log file (fluent.log is sufficient) Only in case a UDF is loaded in the journal file
#rm log
