#!/bin/sh
#
# Locate one of the events DPRK1, DPRK3, DPRK4 and DPRK5 
# relative to the fixed DPRK2
#
scriptname=./locate_one_event.sh
if [ $# != 1 ]
then
  echo
  echo "USAGE: "
  echo "$scriptname   evt_to_locate "
  echo "$scriptname     DPRK1       "
  echo "$scriptname     DPRK3       "
  echo "$scriptname     DPRK4       "
  echo "$scriptname     DPRK5       "
  echo
  exit 1
fi
#
event=$1
if [ ${event} = "DPRK2" ]
then
  echo You cannot locate event DPRK2 relative to itself.
  exit 1
fi
GT_file=./DPRK_events.txt
times_file=../DPRK_dataset/DPRK_CC_times.txt
AK135_file=../DPRK_dataset/DPRK_ak135_slovecs.txt
exefile=../VIScERaL_GD.py
if test ! -r ${GT_file}
then
  echo Ground truth file ${GT_file} not found ...
  exit 1
fi
if test ! -r ${times_file}
then
  echo times file ${times_file} not found ...
  exit 1
fi
if test ! -r ${AK135_file}
then
  echo AK135 slowness vector file ${AK135_file} not found ...
  exit 1
fi
if test ! -x ${exefile}
then
  echo Script ${exefile} not found ...
  exit 1
fi
#
# First we make the file slovecs_flags.txt from the AK135 file.
# This just involves adding the letter F to the end of each
# line to make sure that the slowness vector is fixed under
# the whole calculation.
#
awk '{print $0, "F"}' ${AK135_file} > ./slovecs_flags.txt
#
# Now we need to make the file abslocs_flags.txt with event DPRK2 fixed
# and event ${event} set to the latitude and longitude of DPRK2 but
# free to move.
# We copy the origin time of event ${event} into the file
# but we note that this time is never referred to in the calculations.
#
awk '$4 == "DPRK2" {print $0}' ${GT_file} > abslocs_flags.txt
numlines=`wc abslocs_flags.txt | awk '{print $1}'`
if [ $numlines -ne 1 ]
then
  echo Error in generating file abslocs_flags.txt
  exit 1
fi
initlat=`head -1 abslocs_flags.txt | awk '{print $2}'`
initlon=`head -1 abslocs_flags.txt | awk '{print $3}'`
awk '$4 == EVENT' EVENT=${event} ${GT_file} > tmp.txt
numlines=`wc tmp.txt | awk '{print $1}'`
if [ $numlines -ne 1 ]
then
  echo Are you sure there is an event ${event} in the dataset?
  exit 1
fi
eventtime=`awk '{print $1}' tmp.txt`
cat << EOF >> abslocs_flags.txt
${eventtime}  ${initlat}  ${initlon}  ${event}   S
EOF
python ${exefile} \
         --reflat  41.29500   \
         --reflon 129.08000   \
         --writehistory       \
         --maxpositiveslopes    2000    \
         --maxiter  20000     \
         --maxldistkm  3.0    \
         --randomizeloc  1.0  \
         --dtfile  ${times_file}       \
         --abslocflagsfile   abslocs_flags.txt   \
         --slovecsflagsfile  slovecs_flags.txt \
         --locoutfile      ${event}_rellocs.txt  --allowmissing
rm tmp.txt
#
# 41.29500  129.08000
