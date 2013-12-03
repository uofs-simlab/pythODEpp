#!/bin/sh
NP=$(sed -e 's/^.*max_slots=\([0-9]\+\).*$/\1/' < ~/hostfile | awk '{s+=$1} END {print s}')
mpirun --np $NP --hostfile ~/hostfile --bynode ./ivprunner.py $1
./analysis.py $1 newest

