#!/bin/bash
#
#   topmed_gcepull.sh -submit| bamid
#
#	Copy remapped CRAM for a sample from Google Cloud
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcepull
mem=2G
markverb="${me}ed"
constraint="--constraint eth-10g"
qos="--qos=topmed-$me"
slurmp=topmed
realhost=''
gsutil='gsutil -o GSUtil:parallel_composite_upload_threshold=150M'
incominguri='gs://topmed-recabs'
build=38

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test $me $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  # Run this on node where remapped cram lives
  h=`$topmedpath whathost $1 b$build`
  if [ "$h" != "" ]; then
    realhost="--nodelist=$h"
    #qos="--qos=$h-$me"
  fi

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` $me ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Copy remapped CRAM for a sample from Google Cloud"
  exit 1
fi
bamid=$1

#   Get remapped cram file location
crampath=`$topmedpath wherepath $bamid b$build`
if [ "$crampath" = "" ]; then
  echo "Unable to determine where remapped CRAM file should go for '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
s=`hostname`
p=`pwd`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid cramfile=$cramfile pwd=$p ========="

#   Mark this as started
$topmedcmd mark $bamid $markverb started

#======================================================================
#   Copy remapped CRAM from Google Cloud
#======================================================================
nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Unable to get expt_sampleid for bamid '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
cramflagstat=`$topmedcmd show $bamid cramflagstat`
if [ "$cramflagstat" = "" ]; then
  echo "Unable to get cramflagstat for bamid '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

stime=`date +%s`
export BOTO_CONFIG=/net/topmed/working/shared/tpg_gsutil_config.txt
mkdir -p $crampath
if [ "$?" != "0" ]; then
  echo "Unable to create directory for remapped data '$crampath' for bamid '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
echo "$nwdid data will be copied to $crampath"

echo "Checking if flagstat is as we expect"
$gsutil cp  $incominguri/$nwdid/$nwdid.recab.cram.flagstat $crampath
if [ "$?" != "0" ]; then
  echo "Failed to copy file to Google Cloud"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
#   Get number of interest from flagstat file and check it
n=`grep 'paired in sequencing' $crampath/$nwdid.recab.cram.flagstat | awk '{print $1}'`
if [ "$n" != "$cramflagstat" ]; then
  echo "Flagstat '$n' did not match cramflagstat '$cramflagstat' for bamid '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi

echo "Copying remapped CRAM to local file"
$gsutil cp $incominguri/$nwdid/$nwdid.recab.cram $crampath
if [ "$?" != "0" ]; then
  echo "Failed to copy file from Google Cloud to $crampath"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of CRAM from Google CLoud to $crampath completed in $etime seconds"
$topmedcmd -persist mark $bamid $markverb completed
$topmedcmd -persist mark $bamid gceposted requested
$topmedcmd -persist set $bamid b${build}flagstat $cramflagstat
echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit
