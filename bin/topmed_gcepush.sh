#!/bin/bash
#
#   topmed_gcepush.sh -submit| bamid [cramfile]
#
#	Copy CRAM for a sample to Google Cloud
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
me=gcepush
markverb=$me
incominguri='gs://topmed-incoming'

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='2:00:00'
  SubmitJob $bamid "$PROJECT-gce" '4G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
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

center=`GetDB $bamid center`
if [ "$center" = "" ]; then
  Fail "Unable to get center for bamid '$bamid'"
fi
run=`GetDB $bamid run`
if [ "$run" = "" ]; then
  Fail "Unable to get run for bamid '$bamid'"
fi
datayear=`GetDB $bamid datayear`
if [ "$datayear" = "" ]; then
  Fail "Unable to get datayear for bamid '$bamid'"
fi
nwdid=`GetNWDID $bamid`
stime=`date +%s`

#======================================================================
#   Special hack - year 4 topmed is not remapped, but we pretend it is
#   Note: when G changes his mind, we must delete all symlinks this creates
#   and realize, of course, SOME year 4 BROAD files were already remapped
#
#   We also do not remap wave6 for inpsyght
#======================================================================
remap='Y'
if [ "$PROJECT" = "topmed" -a "$datayear" = "4" -a "$center" = "broad" ]; then
  remap='N'
fi
if [ "$PROJECT" = "inpsyght" -a "$run" = "wave6" ]; then
  remap='N'
fi

if [ "$remap" = "N" ]; then
  echo "Special hack for $PROJECT $center $run datayear $datayear - pretend remapping was done"
  recabfile=`$topmedpath wherefile $bamid b38`
  if [ "$recabfile" = "" ]; then
    Fail "Unable to find path for recab file for $nwdid [$bamid]"
  fi
  if [ -f $recabfile ]; then
    rm -f $recabfile $recabfile.crai
    echo "Recab $recabfile for $nwdid [$bamid] already exists, removed"
  fi
  #   No local recab exists, make fake one
  p=`dirname $recabfile`
  mkdir -p $p
  if [ "$?" != "0" ]; then
    Fail "Unable to create recab directory $nwdid [$bamid]"
  fi
  ln -s $cramfile $recabfile            # Finally make symlink
  if [ "$?" != "0" ]; then
    Fail "Unable to create recab symlink $recabfile for $nwdid [$bamid]"
  fi
  ln -s $cramfile.crai $recabfile.crai
  if [ "$?" != "0" ]; then
    Fail "Unable to create recab symlink $recabfile.crai for $nwdid [$bamid]"
  fi
  x=`GetDB $bamid bamflagstat`
  SetDB $bamid b38flagstat $x
  x=`GetDB $bamid checksum`
  SetDB $bamid b38cramchecksum $x
  x=`CalcMD5 $bamid $recabfile.crai`
  SetDB $bamid b38craichecksum $x
  $topmedcmd setdate $bamid datemapping_b38 $recabfile
  SetDB $bamid state_gce38pull 20       # B38 remapping is done

  #   Done with fake pull, finish up
  SetDB $bamid state_gce38bcf 0
  SetDB $bamid state_b38 20
  SetDB $bamid state_aws38copy 0
  if [ "$PROJECT" = "topmed" ]; then
    SetDB $bamid state_gce38copy 1        # Be sure this gets copied to GCE
  fi

  etime=`date +%s`
  etime=`expr $etime - $stime`
  Successful
  Log $etime
  exit
fi

#======================================================================
#   Copy CRAM to Google Cloud
#======================================================================
echo "Copying CRAM to $incominguri/$center/$run/$nwdid.src.cram"
$gsutil cp $cramfile $incominguri/$center/$run/$nwdid.src.cram
if [ "$?" != "0" ]; then
  Fail "Failed to copy file to GCE"
fi

#   Tell the Google Cloud processes that a new file has shown up
tmperr=/run/shm/$$.curlstderr
tmpout=/run/shm/$$.curloutput

curl -f -i -u `cat /usr/cluster/topmed/etc/.db_connections/104.198.71.226.cred` --insecure --data $nwdid --stderr $tmperr --output $tmpout  https://104.198.71.226/api/unprocessed-samples
rc=$?
cat $tmpout $tmperr                         # So curl results are in log
a=`grep 'Unprocessable Entity' $tmperr`     # Ignore errors we do not care about
rm -f $tmperr $tmpout
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
SetDB $bamid state_gce38bcf 0
SetDB $bamid state_gce38pull 0      # Mark B38 remapping as not done yet
SetDB $bamid state_b38 0
SetDB $bamid state_gce38copy 0      # Mark copy files to GCE as not done yet
SetDB $bamid state_aws38copy 0

Successful
Log $etime
exit
