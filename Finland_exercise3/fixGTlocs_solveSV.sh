#!/bin/sh
#
# Steven J Gibbons 2023/06/26
#
# Using the GT-locations for the explosion as fixed, we try to
# solve for the slowness vecotors that best fit these solutions.
#
scriptname=./fixGTlocs_solveSV.sh
exefile=../VIScERaL_GD.py
if test ! -x ${exefile}
then
  echo Script ${exefile} not found ...
  exit 1
fi
#
# Now we want to use the ak135 vectors but allow
# them to change - it is these that we want to solve for
# with the fixed event locations.
#
awk '{print $0, " S " }' ../Finland_dataset/Finland_ak135_slovecs.txt > var_slovecs_flags.txt
#
python ${exefile} \
         --reflat  67.93590   \
         --reflon  25.83491   \
         --writehistory       \
         --maxiter  10000     \
         --maxpositiveslopes     500        \
         --dtfile  ../Finland_dataset/Finland_CC_times.txt \
         --abslocflagsfile   ../Finland_dataset/Finland_events_GT.txt \
         --slovecsflagsfile  var_slovecs_flags.txt \
         --slooutfile      solveSVonly_forFixedGT.txt  --allowmissing
#
