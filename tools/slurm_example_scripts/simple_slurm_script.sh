#!/bin/bash
#SBATCH --partition=x-men
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=950
#SBATCH --time=0-2:30:00

for i in $(seq 1 32)
do
./sample.exe DecT100.dat sam$i 1 &
sleep 1
done
wait
