#!/bin/bash
#
#   topmed_gcepull.sh -submit bamid -retry
#
#	Copy remapped CRAM for a sample from GCE
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
topmed_check_recab=/usr/cluster/$PROJECT/bin/topmed_check_recab.pl

me=gcepull
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  MyRealHost $bamid b$build
  MayIRun $me $bamid $realhost
  timeout='4:00:00'
  SubmitJob $bamid "$PROJECT-gce" '4G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid [-retry]"
  echo ""
  echo "Copy remapped CRAM for a sample from GCE"
  exit 1
fi
bamid=$1
nwdid=`GetNWDID $bamid`
bamid=`GetDB $nwdid bamid`

retry=0
if [ "$2" = "-retry" ]; then        # See if nomatch really was okay
  retry=1
fi

#   Get where our remapped file is to go
crampath=`$topmedpath wherepath $bamid b$build`
if [ "$crampath" = "" ]; then
  Fail "Unable to determine where remapped CRAM file should go for '$bamid'"
fi
mkdir -p $crampath
if [ "$?" != "0" ]; then
  Fail "Unable to create '$crampath' for remapped CRAM for '$bamid'"
fi

Started
cramfile=$crampath/$nwdid.recab.cram

#======================================================================
#   Copy remapped CRAM from GCE, check flagstat, fix up database
#======================================================================
cramflagstat=`GetDB $bamid cramflagstat`
if [ "$cramflagstat" = "" ]; then
  SetDB $bamid state_gce38bcf 0
  SetDB $bamid state_gce38copy 0
  Fail "Cramflagstat is missing from database for bamid '$bamid'"
fi

stime=`date +%s`

inuri="$incominguri/$nwdid/$nwdid.recab.cram"
$gsutil stat "$inuri"
if [ "$?" != "0" ]; then        # Remote file is not there
  Fail "Unable to find $nwdid/$nwdid.recab.cram in: $incominguri"
fi

#   If retrying, try to rename the possible nomatch file. Failure okay
if [ "$retry" = "1" ]; then
  $gsutil mv $inuri.flagstat.nomatch $inuri.flagstat
fi

#   Now know where to look for data. Check flagstat
echo "Checking if flagstat is as we expect from $inuri"
$gsutil cp  $inuri.flagstat $crampath
if [ "$?" != "0" ]; then
  SetDB $bamid state_gce38bcf 0
  SetDB $bamid state_gce38copy 0
  Fail "Failed to copy flagstat from GCE: $inuri.flagstat"
fi
#   Get number of interest from flagstat file and check it
n=`CalcFlagstatFromFile $crampath/$nwdid.recab.cram.flagstat`
if [ "$n" != "$cramflagstat" ]; then
  # Renaming the flagstat file stops pull from happening again
  $gsutil mv $inuri.flagstat $inuri.flagstat.nomatch
  SetDB $bamid state_gce38bcf 0
  SetDB $bamid state_gce38copy 0
  Fail "Flagstat '$n' did not match cramflagstat '$cramflagstat' for bamid '$bamid' nwdid $nwdid  -- URL=$inuri"
fi
echo "Flagstat value is correct: $n"

#   See if we have already done this
f=$crampath/$nwdid.recab.cram
if [ -f $f ]; then
  echo "Replacing existing CRAM $f"
  rm -f $f $f.crai
fi

echo "Copying remapped CRAM to local file $crampath"
$gsutil cp $inuri $f
if [ "$?" != "0" ]; then
  SetDB $bamid state_gce38bcf 0
  SetDB $bamid state_gce38copy 0
  Fail "Failed to copy file from GCE $inuri to $f"
fi

#   Remapping can still result in a trashed file.  Make sure this is a CSG file
set -o pipefail
$samtools view -H $f | $topmed_check_recab -csg
if [ "$?" != "0" ]; then
  SetDB $bamid state_gce38bcf 0
  SetDB $bamid state_gce38copy 0q
  Fail "Remapped file '$f' header has multiple ids"
fi
echo "Only one sample found in the header and is a CSG remapped file"

#   Clean up data in GCE if data found in incoming.  Move remapped data to bcf bucket
$gsutil mv $inuri         $bcfuri/$nwdid/$nwdid.recab.cram
echo "Moved $inuri files to $bcfuri/$nwdid"
#   Remove any left over cruft in recabs bucket
echo "Removing $incominguri/$nwdid"
$gsutil rm -rf $incominguri/$nwdid

#   Post processing needed here
echo "Begin post-processing of $f"
echo "Create index for remapped sample"
CreateIndex $bamid $f
echo "Calculating MD5s for local files"
md5cram=`CalcMD5 $bamid $f`
md5crai=`CalcMD5 $bamid $f.crai`

echo "Set checksums and flagstat for b$build sample"
SetDB $bamid b${build}cramchecksum $md5cram
SetDB $bamid b${build}craichecksum $md5crai
SetDB $bamid b${build}flagstat $cramflagstat

#   Save date of file in database
$topmedcmd setdate $bamid datemapping_b38 $f

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Copy of remapped CRAM from GCE to $crampath completed in $etime seconds"

SetDB $bamid state_b${build} 20     # Mark b38 as done
SetDB $bamid state_gce38bcf 1       # We need more reprocessing
SetDB $bamid state_gce38copy 0
Successful
Log $etime
exit
