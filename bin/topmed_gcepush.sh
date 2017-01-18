#!/bin/bash
#
#   topmed_gcepush.sh -submit| bamid [cramfile]
#
#	Copy CRAM for a sample to Google Cloud
#
bindir=/usr/cluster/bin
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
topmedpath=/usr/cluster/monitor/bin/topmedpath.pl
medir=`dirname $0`
mem=2G
console=/net/topmed/working/topmed-output
markverb=gcepushed
constraint="--constraint eth-10g"
qos='--qos=topmed-gcepush'
slurmp=topmed
realhost=''
gsutil='gsutil -o GSUtil:parallel_composite_upload_threshold=150M'
incominguri='gs://topmed-incoming'

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test gcepush $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  # Run this on node where cram lives
  h=`$topmedpath whathost $1 cram`
  if [ "$h" != "" ]; then
    realhost="--nodelist=$h"
    #qos="--qos=$h-gcepush"
  fi

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-gcepush --output=$console/$1-gcepush.out $0 $sq $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-gcepush --output=$console/$1-gcepush.out $0 $sq $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` gcepush ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid [cramfile]"
  echo ""
  echo "Copy CRAM for a sample to Google Cloud"
  exit 1
fi
bamid=$1
cramfile=$2

#   Get cram file if it was not provided
if [ "$cramfile" = "" ]; then
  cramfile=`$topmedpath wherefile $bamid cram`
fi
if [ "$cramfile" = "" ]; then
  echo "Unable to determine CRAM file for '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
s=`hostname`
p=`pwd`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid cramfile=$cramfile pwd=$p ========="

#   Mark this as started
$topmedcmd mark $bamid $markverb started
stime=`date +%s`

#======================================================================
#   Copy CRAM to Google Cloud
#======================================================================
center=`$topmedcmd show $bamid center`
if [ "$center" = "" ]; then
  echo "Unable to get center for bamid '$bamid'"
  exit 2
fi
run=`$topmedcmd show $bamid run`
if [ "$run" = "" ]; then
  echo "Unable to get run for bamid '$bamid'"
  exit 2
fi
nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Unable to get expt_sampleid for bamid '$bamid'"
  exit 2
fi

stime=`date +%s`
echo "Copying CRAM to $incominguri/$center/$run/$nwdid.src.cram"
export BOTO_CONFIG=/net/topmed/working/shared/tpg_gsutil_config.txt
$gsutil cp $cramfile $incominguri/$center/$run/$nwdid.src.cram
if [ "$?" != "0" ]; then
  echo "Failed to copy file to Google Cloud"
  exit 3
fi

#   Tell the Google Cloud processes that a new file has shown up
curl -f -i -u "csg:WV9kNT35udEE6B9Q" --insecure --data $nwdid https://104.198.71.226/api/unprocessed-samples
if [ "$?" != "0" ]; then
  echo "Failed to notify Google Cloud that new file arrived"
  exit 3
fi

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of CRAM to Google CLoud completed in $etime seconds"
$topmedcmd -persist mark $bamid $markverb completed
echo `date` gcepush $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit
