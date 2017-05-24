#!/bin/bash
#
#   topmed_redo.sh -submit| bamid
#
#	Do some sort of processing for a sample - usually cause I missed something
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=redo
markverb="${me}ed"
incominguri='gs://topmed-recabs'

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  #MayIRun $me  $bamid
  RandomRealHost $bamid
  SubmitJob $bamid "$realhost-verify" '2G' "$0 $*"
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

Started

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
SetDB $bamid b38cramchecksum $md5
echo "Set checksum for b$build file $bamid ($md5)"

echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit
