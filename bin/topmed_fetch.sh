#!/bin/bash
#
#   topmed_fetch.sh -submit destrun rmtfile
#
#   e.g. topmed_fetch.sh -submit 2016jun23 fetch a/NWD548276.bam
#
#	Fetch a file from the remote site. The destrun is the dirname
#   for the file and provides detailed information on how/where
#   to fetch the file when 'topmedaspera.pl manifest  destrun' was run.
#   That results in a list of files to fetch obtained from the Manifest.txt
#   and an invocation of this script.
#
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
topmedaspera=/usr/cluster/monitor/bin/topmedaspera.pl
topmednwd=/usr/cluster/monitor/bin/topmed_nwdid.pl
topmedrename=/usr/cluster/monitor/bin/topmedrename.pl
console=/net/topmed/working/topmed-output
mem=2G                  # Artificially high so not too many on small nodes
markverb=arrived
slurmp=topmed
qos=topmed-ncbi
realhost='--nodelist=topmed'

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test arrive $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  # Run this on node where Manifest lives
  #l=(`$topmedcmd where $1 bam`)
  #if [ "${l[1]}" != "" ]; then
  #  realhost="--nodelist=${l[1]}"
  #  qos="${l[1]}-ncbi"
  #fi

  # This could be run before the Manifest has been discovered
  # which means there might not be a bamid yet.
  dest=$1
  bamid=`$topmedcmd whatbamid $2`
  j=$bamid
  if [ "$bamid" = "" ]; then
    j=`basename $2`
    j="${j%.*}"
  fi

  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos $realhost -J $j-arrive --output=$console/$j-arrive.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos $realhost -J $j-arrive --output=$console/$j-arrive.out $0 $*"
    exit 1
  fi
  if [ "$bamid" != "" ]; then $topmedcmd mark $bamid $markverb submitted; fi 
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` arrive ${l[3]} $slurp $slurmqos $mem >> $console/$j.jobids
  fi
  exit
fi

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] dest rmtfile"
  echo ""
  echo "Fetch a file from the remote site and update the database"
  exit 1
fi
dest=$1
rmtfile=$2

#   This could be run before the Manifest has been discovered
#   which means there might not be a bamid yet.
#   This only works if the file is named NWDxxxxxx.bam
nwdfile=`basename $rmtfile`     # Get NWD from a/NWD123456.bam
bamid=`$topmedcmd whatbamid $nwdfile`

#   If this has a bamid, check if it has already arrived
if [ "$bamid" != "" ]; then
  arrived=`$topmedcmd show $bamid state_arrive`
  if [ "$arrived" != "0" -a "$arrived" != "1" -a "$arrived" != "2" ]; then
    echo "File '$nwdfile' [$bamid] already fetched or in process (state_arrive=$arrived)"
    exit 2
  fi
fi
g
#   Mark this as started
if [ "$bamid" != "" ]; then $topmedcmd mark $bamid $markverb started; fi
d=`date +%Y/%m/%d`
s=`hostname`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid dest=$dest rmtfile=$rmtfile ========="

#   We seem to have problems getting a working connection to the broad
#   ascp connects, but sometimes never actually delivers data.
#   See if waiting a bit will help
sleep `shuf -i 10-60 -n 1`

stime=`date +%s`
$topmedaspera -dest $dest fetch $rmtfile
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Command completed in $etime seconds"
if [ "$rc" != "0" ]; then
  if [ "$bamid" != "" ]; then $topmedcmd mark $bamid $markverb failed; fi
  exit 1
fi
if [ "$bamid" != "" ]; then
    #   Set NWDID and other database fields
    p=(`$topmedcmd where $bamid bam`)       # Dir for bam
    bampath=${p[0]}/$nwdfile                # FQP to bam
    $topmednwd -bamid $bamid $bampath
    if [ "$?" != "0" ]; then
      $topmedcmd mark $bamid arrived failed
      exit 3
  fi
  #   Rename the bamfile to NWD and fix the database
  $topmedrename $bamid $bampath
  if [ "$?" != "0" ]; then
    $topmedcmd mark $bamid arrived failed
    exit 4
  fi
  $topmedcmd mark $bamid $markverb completed
fi
echo `date` arrive $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids

exit
