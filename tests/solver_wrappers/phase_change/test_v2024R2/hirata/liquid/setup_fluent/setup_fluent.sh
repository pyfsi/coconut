#!/bin/bash

# make fluent case with 4 cores for initialisation calculations
fluent 2ddp -t4 -gu -i case.jou > setup_fluent.log 2>&1

# delete log file (fluent.log is sufficient) Only in case a UDF is loaded in the journal file
#rm log
