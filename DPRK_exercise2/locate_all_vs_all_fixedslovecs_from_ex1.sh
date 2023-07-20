#!/bin/sh
scriptname=./locate_all_vs_all_fixedslovecs_from_ex1.sh
exefile=../VIScERaL_GD.py
if test ! -x ${exefile}
then
  echo Script ${exefile} not found ...
  exit 1
fi
#
# We copy the file slovecs_flags.txt from the directory
# DPRK_exercise1 - these are the fixed ak135 slowness vectors.
# We also copy the file rellocs_vs_DPRK2.txt
# from DPRK_exercise1 - this is the output from the 
# calculations where all events are located one by one
# relative to DPRK2.
# Note that DPRK2 is still fixed (flag F) and all the 
# other events have flag S - so they can move.
#
# We now simply run the script VIScERaL_GD.py
#
python ${exefile} \
         --reflat  41.29500   \
         --reflon 129.08000   \
         --writehistory       \
         --maxpositiveslopes     500        \
         --maxldistkm    1.500      \
         --randomizeloc  0.000    \
         --maxiter  100000     \
         --dtfile  ../DPRK_dataset/DPRK_CC_times.txt \
         --abslocflagsfile   ./rellocs_vs_DPRK2_fromex1.txt  \
         --slovecsflagsfile  slovecs_flags.txt \
         --locoutfile      allvsall_fixedslo_rellocs.txt  --allowmissing
#
