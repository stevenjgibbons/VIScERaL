#!/bin/sh
scriptname=./locate_one_event.sh
if [ $# != 1 ]
then
  echo
  echo "USAGE: "
  echo "$scriptname   evt_to_locate "
  echo "$scriptname       H02       "
  echo
  exit 1
fi
#
event=$1
if [ ${event} = "H01" ]
then
  echo You cannot locate event H01 relative to itself.
  exit 1
fi
GT_file=../Finland_dataset/Finland_events_GT.txt
DT_file=../Finland_dataset/Finland_CC_times.txt
AK135_file=../Finland_dataset/Finland_ak135_slovecs.txt
exefile=../VIScERaL_GD.py
if test ! -r ${DT_file}
then
  echo Ground truth file ${DT_file} not found ...
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
# Now we need to make the file abslocs_flags.txt with event H01 fixed
# and event ${event} set to the latitude and longitude of H01 but
# free to move.
# We copy the origin time of event ${event} into the file
# but we note that this time is never referred to in the calculations.
#
awk '$4 == "H01" {print $0}' ${GT_file} > abslocs_flags.txt
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
         --reflat  67.93590   \
         --reflon  25.83491   \
         --writehistory       \
         --maxiter  10000     \
         --maxpositiveslopes    500    \
         --dtfile  ${DT_file} \
         --abslocflagsfile   abslocs_flags.txt   \
         --slovecsflagsfile  slovecs_flags.txt \
         --nrmoutfile      ${event}_nrmvalues.txt \
         --locoutfile      ${event}_rellocs.txt  --allowmissing
rm tmp.txt
#
