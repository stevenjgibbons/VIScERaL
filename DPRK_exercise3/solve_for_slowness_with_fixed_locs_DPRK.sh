#!/bin/sh
scriptname=./solve_for_slowness_with_fixed_locs_DPRK.sh
exefile=../VIScERaL_GD.py
if test ! -x ${exefile}
then
  echo Script ${exefile} not found ...
  exit 1
fi
#
# We copy the file slovecs_flags.txt from the directory
# DPRK_exercise1 - these are the fixed ak135 slowness vectors.
# We also copy the file rellocs_vs_D2009.txt
# from DPRK_exercise1 - this is the output from the 
# calculations where all events are located one by one
# relative to D2009.
# Note that D2009 is still fixed (flag F) and all the 
# other events have flag S - so they can move.
#
# We now simply run the script VOWEL.py
#
python ${exefile} \
         --reflat  41.29500   \
         --reflon 129.08000   \
         --writehistory       \
         --FixEvent0       \
         --maxpositiveslopes     5000       \
         --delkm      0.050      \
         --delslow      0.005      \
         --maxsdist      0.050      \
         --randomizeslo  0.000      \
         --maxldistkm    1.500      \
         --randomizeloc  0.000    \
         --maxiter  100000     \
         --dtfile  ../DPRK_dataset/DPRK_CC_times.txt \
         --abslocflagsfile   rellocs_vs_DPRK2_fromex1_fixed.txt   \
         --slovecsflagsfile  slovecs_flags_solve.txt \
         --locoutfile      allvsall_varslo_fromex1.txt  --allowmissing
#
