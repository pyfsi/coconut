cd CFD
#TODO remove module load command
ml ANSYS_CFD/2019R3
fluent 2ddp -gu -i movie.jou > movie.log 2>&1
ml -ANSYS_CFD
cd ..