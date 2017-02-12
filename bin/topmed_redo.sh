#!/bin/bash
#
#   topmed_redo.sh -submit| bamid
#
#	Do some sort of processing for a sample - usually cause I missed something
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=redo
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
  #$topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` $me ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Do some sort of processing for a sample"
  exit 1
fi
bamid=$1

d=`date +%Y/%m/%d`
s=`hostname`
p=`pwd`

echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid ========="

#======================================================================
#   Post process the remapped CRAM from Google Cloud
#======================================================================
#   Mark this as started
#$topmedcmd mark $bamid $markverb started

cramfile=`$topmedpath wherefile $bamid b$build`
if [ "$?" != "0" ]; then
  echo "Unable to get path to CRAM file '$bamid'"
  #$topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
if [ ! -f $cramfile ]; then
  echo "CRAM file does not exist '$bamid' cramfile=$cramfile"
  #$topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

echo "Calculating MD5 for local file ($cramfile)"
md5=(`md5sum $cramfile`)
md5=${md5[0]}
if [ "$md5" = "" ]; then
  echo "Unable to calculate MD5 for remapped '$bamid' [$nwdid] cramfile=$cramfile"
  #$topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
$topmedcmd -persist set $bamid b38cramchecksum $md5
echo "Set checksum for b$build file $bamid ($md5)"

echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit
