#!/bin/bash

PMI_RANK=$1
JOBIDFILE=$2
LOCKFILE=$3
STATUSDIR=$4
MAXNODE=$5

# Create a rank to note that the current job is running
touch $STATUSDIR/$PMI_RANK
BASESTATUSDIR=$(dirname $STATUSDIR)
if [ $PBS_ARRAYID -eq 1 ]; then
  # This means only when all of the jobs were finished, it would be finished.
  MINIMUMJOB=1
else
  MINIMUMJOB=4
fi

WHICHLINE=1
JOBID=0
LASTNODE=N

cd $TMPDIR

function hms
{
  s=$1
  h=$((s/3600))
  s=$((s-(h*3600)));
  m=$((s/60));
  s=$((s-(m*60)));
  printf "%02d:%02d:%02d\n" $h $m $s
}

function checkIfTheJobSuccessfullyFinished 
{
  LOGFILE=${JOBID##*\/}
  BOOTSTATFILE=input/${LOGFILE}_phyml_boot_stats.txt
  BOOTTREEFILE=input/${LOGFILE}_phyml_boot_trees.txt
  STATFILE=input/${LOGFILE}_phyml_stats.txt  
  TREEFILE=input/${LOGFILE}_phyml_tree.txt
  
  FINISHEDOUTPUT=N
  if [ -s $TREEFILE ]; then
      FINISHEDOUTPUT=Y
  fi
  if [ "$FINISHEDOUTPUT" == "Y" ]; then
    # The job is succesfully finished.
    END_TIME=`date +%s`
    ELAPSED=`expr $END_TIME - $START_TIME`
    # echo end-$JOBID
    hms $ELAPSED
    mv $BOOTSTATFILE output
    mv $BOOTTREEFILE output
    mv $STATFILE output
    mv $TREEFILE output
  else
    while [ "$JOBID" != "" ]; do
      lockfile=$LOCKFILE
      if ( set -o noclobber; echo "$$" > "$lockfile") 2> /dev/null; then
        trap 'rm -f "$lockfile"; exit $?' INT TERM
        cp $JOBIDFILE $JOBIDFILE.temp 
        echo $JOBID >> $JOBIDFILE.temp
        mv $JOBIDFILE.temp $JOBIDFILE
        rm -f "$lockfile"
        trap - INT TERM
        # Let the job be finished.
        JOBID=""
      else
        sleep 5
      fi
    done
    JOBID=0
    rm -f $BOOTSTATFILE $BOOTTREEFILE $STATFILE $TREEFILE
  fi
}

# Keep trying to read in jobidfile until the current node is the last one.
while [ "$JOBID" != "" ]; do
  lockfile=$LOCKFILE
  if ( set -o noclobber; echo "$$" > "$lockfile") 2> /dev/null; then
    trap 'rm -f "$lockfile"; exit $?' INT TERM
    read -r JOBID < ${JOBIDFILE}
    sed '1d' $JOBIDFILE > $JOBIDFILE.temp; 
    mv $JOBIDFILE.temp $JOBIDFILE
    rm -f "$lockfile"
    trap - INT TERM

    if [ "$JOBID" == "" ]; then
      NUMJOBS=$(ps ax|grep phyml|grep -v grep|wc -l)
      if [ $NUMJOBS -lt $MINIMUMJOB ]; then
        JOBID=""
      else
        JOBID=0
      fi

    else
      START_TIME=`date +%s`
      eval $JOBID
      checkIfTheJobSuccessfullyFinished
    fi
  else
    JOBID=0
    sleep 10
  fi
done

# Wait for all of the jobs to be finished.
rm -f $STATUSDIR/$PMI_RANK
NUMSTATUS=$(ls -1 $STATUSDIR|wc -l)
if [ $PBS_ARRAYID -gt 1 ]; then
  while [ $NUMSTATUS -gt 0 ]; do
    # Kill those jobs and put JOBID back to jobidfile.
    NUMJOBS=$(ps ax|grep phyml|grep -v grep|wc -l)
    if [ $NUMJOBS -lt $MINIMUMJOB ]; then
      ps ax | grep phyml | grep -v grep | awk '{print $1}' | xargs -i kill {} 2&>/dev/null
    fi
    NUMSTATUS=$(ls -1 $STATUSDIR|wc -l)
    # sleep 5
    sleep 60
  done
fi
