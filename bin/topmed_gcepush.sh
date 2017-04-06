#!/bin/bash
#
#   topmed_gcepush.sh -submit| bamid [cramfile]
#
#	Copy CRAM for a sample to Google Cloud
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcepush
mem=2G
markverb=$me
cores="--cpus-per-task=2"           # Cores here should be same as gsutil
gsutil='gsutil -o GSUtil:parallel_composite_upload_threshold=150M -o GSUtil:parallel_process_count=2'
incominguri='gs://topmed-incoming'
qos="--qos=topmed-$me"
realhost=''

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedpermit permit test $me $1
  if [ "$?" = "0" ]; then
    echo "$me $1 not permitted" | tee $console/$1-$me.out
    exit 4
  fi 

  # Run this on node where cram lives, unless it is csgspare and then make up something
  h=`$topmedpath whathost $1 cram`
  if [ "$h" != "" ]; then
    if [ "$h" = "csgspare" ]; then
      n=`perl -e "{ print $1%7 }"`
      h="topmed$n"
      if [ "$h" = "topmed0" ]; then
        h=topmed2
      fi
      if [ "$h" = "topmed1" ]; then
        h=topmed4
      fi
    fi
    realhost="--nodelist=$h"
    #qos="--qos=$h-$me"
  fi

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $cores $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM - $l" > $console/$1-$me.out
    $topmedcmd -emsg "Failed to submit command to SLURM - $l" mark $1 $markverb failed
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $cores $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*" >> $console/$1-$me.out
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
  $topmedcmd -persist -emsg "Unable to determine CRAM file for '$bamid'" mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
s=`hostname`
p=`pwd`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid cramfile=$cramfile pwd=$p ========="

#   Mark this as started
$topmedcmd mark $bamid $markverb started

#======================================================================
#   Copy CRAM to Google Cloud
#======================================================================
center=`$topmedcmd show $bamid center`
if [ "$center" = "" ]; then
  $topmedcmd -persist -emsg "Unable to get center for bamid '$bamid'" mark $bamid $markverb failed
  exit 2
fi
run=`$topmedcmd show $bamid run`
if [ "$run" = "" ]; then
  $topmedcmd -persist -emsg "Unable to get run for bamid '$bamid'" mark $bamid $markverb failed
  exit 2
fi
nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  $topmedcmd -persist -emsg "Unable to get expt_sampleid for bamid '$bamid'" mark $bamid $markverb failed
  exit 2
fi

stime=`date +%s`
echo "Copying CRAM to $incominguri/$center/$run/$nwdid.src.cram"
export BOTO_CONFIG=/net/topmed/working/shared/tpg_gsutil_config.txt
$gsutil cp $cramfile $incominguri/$center/$run/$nwdid.src.cram
if [ "$?" != "0" ]; then
  $topmedcmd -persist -emsg "Failed to copy file to GCE" mark $bamid $markverb failed
  exit 3
fi

#   Tell the Google Cloud processes that a new file has shown up
tmperr=/tmp/$$.curlstderr
tmpout=/tmp/$$.curloutput
curl -f -i -u "csg:WV9kNT35udEE6B9Q" --insecure --data $nwdid --stderr $tmperr --output $tmpout  https://104.198.71.226/api/unprocessed-samples
rc=$?
cat $tmpout $tmperr                         # So curl results are in log
a=`grep 'Unprocessable Entity' $tmperr`     # Ignore errors we do not care about
rm $tmperr $tmpout
if [ "$a" = "" ]; then
  if [ "$rc" != "0" ]; then
    $topmedcmd -persist -emsg "Failed to notify Google Cloud that new file arrived" mark $bamid $markverb failed
    exit 3
  fi
else
  echo "I think an Unprocessable Entity means this file was just being sent by accident"
  echo "If that is not the case, get the sample in GCE marked to be reprocessed and repush"
fi

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of CRAM to Google CLoud completed in $etime seconds"
$topmedcmd -persist mark $bamid $markverb completed
echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit
