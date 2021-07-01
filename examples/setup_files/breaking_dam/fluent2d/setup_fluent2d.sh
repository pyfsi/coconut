#!/bin/bash

# make fluent case
fluent 2ddp -g -i case.jou > setup_fluent.log 2>&1
fluent 2ddp -g -i case_expl.jou > setup_fluent_expl.log 2>&1
