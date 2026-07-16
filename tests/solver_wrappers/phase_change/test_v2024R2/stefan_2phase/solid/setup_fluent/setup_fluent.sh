#!/bin/bash

# make fluent case
fluent 2ddp -g -i case.jou > setup_fluent.log 2>&1

# delete log file (fluent.log is sufficient) Only in case a UDF is loaded in the journal file
#rm log
