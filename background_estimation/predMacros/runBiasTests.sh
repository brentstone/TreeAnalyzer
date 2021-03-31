#!/bin/bash
# Usage: . runBiasTest.sh <path_to_doBiasTest.py> <n 10s of toys> <r_inj> <mass>

scriptPath=$1
ccPath=$2
n=$3
rinj=$4
mass=$5

r2=$(awk -v r=$rinj -v mult=2 'BEGIN{rf=(r*mult); print rf;}')
r5=$(awk -v r=$rinj -v mult=5 'BEGIN{rf=(r*mult); print rf;}')

python ${scriptPath}/doBiasTest.py -i ${ccPath}/combinedCard.txt -N ${n} -r $rinj -m ${mass} -l Postfit_r1 -f 1 -s 0 -o 0 -c "combine" > out${mass}_r1 &
python ${scriptPath}/doBiasTest.py -i ${ccPath}/combinedCard.txt -N ${n} -r $r2 -m ${mass} -l Postfit_r2 -f 1 -s 0 -o 0 -c "combine" > out${mass}_r2 &
python ${scriptPath}/doBiasTest.py -i ${ccPath}/combinedCard.txt -N ${n} -r $r5 -m ${mass} -l Postfit_r5 -f 1 -s 0 -o 0 -c "combine" > out${mass}_r5 &

echo $rinj
echo $r2
echo $r5

