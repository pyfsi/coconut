#!/bin/bash

# make fluent case
fluent 3ddp -t40 -gu -i case.jou > setup_fluent.log 2>&1

# delete log file (fluent.log is sufficient)
rm report-steady.out
