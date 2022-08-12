#!/bin/bash

cp movie.jou ../
cd ..

ml ANSYS_CFD/2021R2
fluent 2ddp -t1 -gu -i movie.jou > movie.log 2>&1
ml -ANSYS_CFD

rm movie.jou
mv movie.log create_movie/movie.log
