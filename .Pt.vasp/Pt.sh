#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -pe mpi 4
#$ -j y
#$ -cwd
#$ -o log/

t0=$(date +%s.%N)
t0_string=$(date)

./run.py --clean
./run.py --pw

t1=$(date +%s.%N)
t1_string=$(date)

t=$(echo "$t1 - $t0"|bc)
h=$(echo "($t/3600)"|bc)
m=$(echo "($t%3600)/60"|bc)
s=$(echo "($t%3600)%60"|bc)

echo ""
echo "# Time Start   : $t0_string"
echo "# Time End     : $t1_string"
echo "# Time Elapsed : ${h}h ${m}m ${s}s"
echo ""
