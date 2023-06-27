#!/bin/sh
scriptname=./locate_all_vs_all_variableslovecs.sh
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
# We copy the file slovecs_flags.txt from the directory
# Finland_exercise1 - these are the fixed ak135 slowness vectors.
# We also copy the file rellocs_vs_H01.txt
# from Finland_exercise1 - this is the output from the 
# calculations where all events are located one by one
# relative to H01.
# Note that H01 is still fixed (flag F) and all the 
# other events have flag S - so they can move.
#
# We now simply run the script VOWEL.py
# I tried --delslow --delkm 1.0 and this did not converge
#
python ${exefile} \
         --reflat  67.93590   \
         --reflon  25.83491   \
         --maxiter  200000     \
         --maxsdist      0.050      \
         --randomizeslo  0.025      \
         --maxldistkm    0.150      \
         --randomizeloc  0.050      \
         --maxpositiveslopes     500        \
         --dtfile  ../Finland_dataset/Finland_CC_times.txt \
         --abslocflagsfile    new_rellocs_init.txt  \
         --slovecsflagsfile   new_slovecs_init.txt \
         --slooutfile  slovecs_${runID}.txt \
         --locoutfile  rellocs_${runID}.txt \
         --nrmoutfile  nrmvals_${runID}.txt \
         --allowmissing
#
