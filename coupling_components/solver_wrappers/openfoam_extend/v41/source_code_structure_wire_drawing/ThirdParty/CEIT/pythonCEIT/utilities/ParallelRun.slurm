#!/bin/bash -l
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH --job-name=RVE_drawing
#SBATCH -A ucd01

# Change to the submission directory
cd $SLURM_SUBMIT_DIR


SOURCE=./compression_appended_paths
TARGET=./compression_drawing_5_cell
NCELLS=9
RELAXATION=0.9
NINCREMENTS=1000


LOGS=$TARGET/logs
RESULTS=$TARGET/results

mkdir -p $TARGET
mkdir $LOGS
mkdir $RESULTS

module load conda
source activate python3
module load taskfarm

rm ./tasks
for FILE in $SOURCE/*
do
    BASENAME="$(basename -- $FILE)"
    echo "python ./mesoscopic_python/ceit_fft.py --nCells $NCELLS --avgFSeries $FILE --nIncrements $NINCREMENTS --outputFile $RESULTS/$BASENAME.output --relaxationFactor $RELAXATION &> $LOGS/log.$BASENAME" >> tasks
done

taskfarm tasks
