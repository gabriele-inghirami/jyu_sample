#!/bin/sh
#SBATCH --partition=general1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --time=0-12:30:00
#SBATCH --mail-type=ALL
#SBATCH --array=0-79:40

module load comp/gcc/8.2.0
#module load comp/intel/2019.0.117

offset=0

for i in $(seq 1 40)
do
  K=$(($offset+$SLURM_ARRAY_TASK_ID+$i))
  id=$(printf "%03d" $K)
  ./sample_gcc.exe DecT100.dat gnures_$id 1 &
done
wait

