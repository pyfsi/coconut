#!/bin/bash

# make fluent case
fluent 3ddp -g -i case.jou > setup_fluent.log 2>&1

# delete log file (fluent.log is sufficient)
rm log

