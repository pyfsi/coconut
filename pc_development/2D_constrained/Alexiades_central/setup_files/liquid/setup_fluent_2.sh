#!/bin/bash

# make fluent case
fluent 2ddp -gu -i case_2.jou > setup_fluent.log 2>&1
#fluent 2ddp -i case_2.jou > setup_fluent.log 2>&1


# delete log file (fluent.log is sufficient) Only in case a UDF is loaded in the journal file
#rm log
