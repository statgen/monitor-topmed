#!/bin/bash
#
#   topmed_gcepush.sh -submit| bamid [cramfile]
#
#	Copy CRAM for a sample to Google Cloud
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcepush
markverb=$me
cores="--cpus-per-task=2"           # Cores here should be same as gsutil

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me  $bamid
  RandomRealHost $bamid
  SubmitJob $bamid "topmed-$me" '2G' "$0 $*"
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
  Fail "Unable to determine CRAM file for '$bamid'"
fi

Started
GetNWDID $bamid

#======================================================================
#   Copy CRAM to Google Cloud
#======================================================================
center=`$topmedcmd show $bamid center`
if [ "$center" = "" ]; then
  Fail "Unable to get center for bamid '$bamid'"
fi
run=`$topmedcmd show $bamid run`
if [ "$run" = "" ]; then
  Fail "Unable to get run for bamid '$bamid'"
fi
datayear=`$topmedcmd show $bamid datayear`
if [ "$datayear" = "3" ]; then
  build=`$topmedcmd show $bamid build`
  if [ "$build" = "38" ]; then
    echo "Datayear $datayear build $build jobs should not be remapped (bamid '$bamid')"
    Successful
    exit 0
  else
    echo "Datayear $datayear build $build needs remapping"
  fi
fi

stime=`date +%s`
echo "Copying CRAM to $incominguri/$center/$run/$nwdid.src.cram"
$gsutil cp $cramfile $incominguri/$center/$run/$nwdid.src.cram
if [ "$?" != "0" ]; then
  Fail "Failed to copy file to GCE"
fi

#   Tell the Google Cloud processes that a new file has shown up
tmperr=/tmp/$$.curlstderr
tmpout=/tmp/$$.curloutput

curl -f -i -u `cat /usr/cluster/topmed/etc/.db_connections/104.198.71.226.cred` --insecure --data $nwdid --stderr $tmperr --output $tmpout  https://104.198.71.226/api/unprocessed-samples
rc=$?
cat $tmpout $tmperr                         # So curl results are in log
a=`grep 'Unprocessable Entity' $tmperr`     # Ignore errors we do not care about
rm $tmperr $tmpout
if [ "$a" = "" ]; then
  if [ "$rc" != "0" ]; then
    Fail "Failed to notify Google Cloud that new file arrived"
  else
    echo "Told Jonathan's system there is a new sample to remap"
  fi
else
  echo "I think an Unprocessable Entity means this file was just being sent by accident"
  echo "If that is not the case, get the sample in GCE marked to be reprocessed and repush"
fi

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of CRAM to Google CLoud completed in $etime seconds"
SetDB $bamid state_gce38bcf_push 0  # Mark bcf as not done yet
SetDB $bamid state_gce38bcf_pull 0
SetDB $bamid state_gce38bcf 0

SetDB $bamid state_gce38pull 0      # Mark B38 remapping as not done yet
SetDB $bamid state_b38 0

SetDB $bamid state_gce38copy 0      # Mark copy files to GCE as not done yet

Successful
Log $etime
exit
