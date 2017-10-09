#!/bin/bash
#
#   topmed_gcepull.sh -submit| bamid
#
#	Copy remapped CRAM for a sample from GCE
#
. /usr/cluster/topmed/bin/topmed_actions.inc
topmed_check_recab=/usr/cluster/topmed/bin/topmed_check_recab.pl

me=gcepull
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  MayIRun $me  $bamid
  MyRealHost $bamid "b$build"
  SubmitJob $bamid "topmed-$me" '4G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Copy remapped CRAM for a sample from GCE"
  exit 1
fi
bamid=$1

#   Get remapped cram file location
crampath=`$topmedpath wherepath $bamid b$build`
if [ "$crampath" = "" ]; then
  Fail "Unable to determine where remapped CRAM file should go for '$bamid'"
fi
mkdir -p $crampath
if [ "$?" != "0" ]; then
  Fail "Unable to create '$crampath' for remapped CRAM for '$bamid'"
fi

Started
nwdid=`GetNWDID $bamid`
cramfile=$crampath/$nwdid.recab.cram

#======================================================================
#   Copy remapped CRAM from GCE, check flagstat, fix up database
#======================================================================
cramflagstat=`GetDB $bamid cramflagstat`
if [ "$cramflagstat" = "" ]; then
  Fail "Unable to get cramflagstat for bamid '$bamid'"
fi

stime=`date +%s`

#   Remapped cram could be > one place (arrgh!)  Figure out where it is
#   E.G. unuri=gs://topmed-recabs/NWD947950/NWD947950.recab.cram
inuri=''
p="$incominguri/$nwdid/$nwdid.recab.cram"
$gsutil stat "$p"
if [ "$?" = "0" ]; then
  inuri=$p
fi
if [ "$inuri" = "" ]; then
  Fail "Unable to find $nwdid/$nwdid.recab.cram in: $incominguri $bcfuri"
fi

#   Now know where to look for data. Check flagstat
echo "Checking if flagstat is as we expect from $inuri"
$gsutil cp  $inuri.flagstat $crampath
if [ "$?" != "0" ]; then
  Fail "Failed to copy flagstat from GCE: $inuri.flagstat"
fi
#   Get number of interest from flagstat file and check it
n=`grep 'paired in sequencing' $crampath/$nwdid.recab.cram.flagstat | awk '{print $1}'`
if [ "$n" != "$cramflagstat" ]; then
  # Renaming the flagstat file stops pull from happening again
  $gsutil mv $inuri.flagstat $inuri.flagstat.nomatch
  Fail "Flagstat '$n' did not match cramflagstat '$cramflagstat' for bamid '$bamid' nwdid $nwdid  -- URL=$inuri"
fi
echo "Flagstat value is correct: $n"
SetDB $bamid b38flagstat $n

#   See if we have already done this
f=$crampath/$nwdid.recab.cram
if [ -f $f ]; then
  echo "Replacing existing CRAM $f"
fi

echo "Copying remapped CRAM to local file $crampath"
$gsutil cp $inuri $crampath
if [ "$?" != "0" ]; then
  Fail "Failed to copy file from GCE $inuri $crampath"
fi

#   Remapping can still result in a trashed file
set -o pipefail
$samtools view -H $f | grep '^@RG' | $topmed_check_recab
if [ "$?" != "0" ]; then
  Fail "Remapped file '$f' header has multiple ids"
fi
echo "Only one sample found in the header"

#   Post processing needed here
echo "Calculating MD5 for local file ($cramfile)"
md5=`CalcMD5 $bamid $cramfile`
echo "Set checksum and flagstat for b$build file"
SetDB $bamid b${build}cramchecksum $md5
SetDB $bamid b${build}flagstat $cramflagstat

#   Save date of file in database
$topmedcmd setdate $bamid datemapping_b38 $cramfile

#   Clean up data in GCE if data found in incoming.  Move remapped data to bcf bucket
$gsutil mv $inuri.flagstat  $bcfuri/$nwdid
$gsutil mv $inuri         $bcfuri/$nwdid
echo "Moved $inuri files to $bcfuri/$nwdid"
#   Remove any left over cruft in recabs bucket
echo "Removing $incominguri/$nwdid"
$gsutil rm -r $incominguri/$nwdid

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Copy of remapped CRAM from GCE to $crampath completed in $etime seconds"

######### Temporary - sometimes remapping screwed up header, rebuild if needed
year=`GetDB $bamid datayear`
if [ "$year" != "3" ]; then
  CheckRGMap $bamid
  if [ "$?" != "0" ]; then
    echo "#======================================================================"
    echo "#    RGMAP for cramfile is broken, rebuild"
    echo "#======================================================================"
    SetDB $bamid state_fix 2          # Fix action started
    /usr/cluster/topmed/topmed_rgmap.sh $bamid $markverb
    if [ "$?" != "0" ]; then
      Fail "RGMAP correction failed"
    fi
  fi
  SetDB $bamid state_fix 20           # Mark sample as fixed or no rgmap needed
fi
######### Remove this when we believe remapping is correct (2018)

SetDB $bamid state_b${build} 20     # Mark b38 as done
SetDB $bamid state_gce38bcf 1       # We need more reprocessing
SetDB $bamid state_gce38copy 0
Successful
Log $etime
exit
