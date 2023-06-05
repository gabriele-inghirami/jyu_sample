#!/bin/sh
#SBATCH --partition=test
#SBATCH --cpus-per-task=1
#SBATCH --time=0-0:59:59
#SBATCH --mail-type=ALL
 
in_prefix=gnures_
out_prefix=gccout_
procs=6
num=$(ls $in_prefix* | wc -l)
lowel=$(echo "$num/$procs" | bc)
higel=$(($lowel+1))
rem=$(echo "$num%$procs" | bc)

if (( $num == 0 ))
then
   echo "No usable file found with prefix: "$in_prefix
   exit
fi

function examine() {
  istart=$1
  iend=$(( $istart + $2 -1))
  echo "parameters are: "$istart $iend
  for i in $(seq $istart $iend)
  do
      id=$(printf "%03d" $i)
      echo "launching python script with arguments:"$in_prefix$id $out_prefix$id
      python3 compute_v2_cumulants.py $in_prefix$id $out_prefix$id
  done
}

for i in $(seq 1 $procs)
do
  if (( $i <= $rem ))
  then
    startpoint=$(echo "($i-1)*$higel + 1" | bc)
    echo "launching examine with:"$startpoint $higel
    examine $startpoint $higel&
  else
    startpoint=$(echo "$rem*$higel + ($i-1-$rem)*$lowel + 1" | bc)
    echo "launching examine with:"$startpoint $lowel
    examine $startpoint $lowel&
  fi
done
wait

