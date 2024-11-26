#!/bin/bash

# make fluent case (12 processor cores for steady stae calculation)
fluent 2ddp -t12 -g -i case_2.jou > setup_fluent.log 2>&1

# delete log file (fluent.log is sufficient) Only in case a UDF is loaded in the journal file
#rm log
