#!/bin/sh
scriptname=./locate_all_vs_all_variableslovecs_DPRK.sh
if [ $# != 1 ]
then
  echo
  echo "USAGE: "
  echo "$scriptname   runID    "
  echo "$scriptname   000001    "
  echo
  exit 1
fi
#
runID=$1
exefile=../VIScERaL_GD.py
if test ! -x ${exefile}
then
  echo Script ${exefile} not found ...
  exit 1
fi
#
# Solve for relative locations of DPRK events relative to DPRK2
# and, at the same time, modify the slowness vectors.
# The idea is that the run is given a unique runID and we write
# output for this runID and run the process multiple (many hundreds of)
# times, randomizing the starting values of the event locations and
# the slowness vectors on each run.
#
python ${exefile} \
         --reflat  41.29500   \
         --reflon 129.08000   \
         --maxiter  200000     \
         --maxsdist      0.040      \
         --randomizeslo  0.0050      \
         --maxldistkm    0.600      \
         --randomizeloc  0.200      \
         --maxpositiveslopes     500        \
         --dtfile  ./DPRK_CC_times.txt \
         --abslocflagsfile   new_rellocs_init.txt  \
         --slovecsflagsfile  new_slovecs_init.txt \
         --slooutfile  slovecs_${runID}.txt \
         --locoutfile  rellocs_${runID}.txt \
         --nrmoutfile  nrmvals_${runID}.txt \
         --allowmissing
#
