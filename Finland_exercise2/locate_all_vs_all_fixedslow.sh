#!/bin/sh
#
# Steven J. Gibbons 
# NGI 2023/06/25
#
# Locates all of the events (except for the fixed event H01)
# by Gradient Descent using all differential times between
# all phases and all events. All slowness vectors fixed
# to theoretical ak135 model.
#
scriptname=./locate_all_vs_all_fixedslow.sh
GT_file=../Finland_dataset/Finland_events_GT.txt
initloc_file=./Rellocs_from_FinEx1.txt
AK135_file=../Finland_dataset/Finland_ak135_slovecs.txt
exefile=../VIScERaL_GD.py
if test ! -r ${initloc_file}
then
  echo Ground truth file ${initloc_file} not found ...
  exit 1
fi
if test ! -r ${GT_file}
then
  echo Ground truth file ${GT_file} not found ...
  exit 1
fi
if test ! -r ${AK135_file}
then
  echo AK135 slowness vector file ${AK135_file} not found ...
  exit 1
fi
awk '{print $0, " F "}' ${AK135_file} > slow_flags.txt
if test ! -x ${exefile}
then
  echo Script ${exefile} not found ...
  exit 1
fi
#
python ${exefile} \
         --reflat  67.93590   \
         --reflon  25.83491   \
         --writehistory       \
         --maxiter  10000     \
         --maxpositiveslopes    500    \
         --dtfile  ../Finland_dataset/Finland_CC_times.txt \
         --abslocflagsfile   ${initloc_file}   \
         --slovecsflagsfile  slow_flags.txt    \
         --nrmoutfile      nrmvalues.txt \
         --locoutfile      rellocs.txt  --allowmissing
#
