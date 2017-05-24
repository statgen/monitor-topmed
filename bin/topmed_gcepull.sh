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
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me  $bamid
  MyRealHost $bamid "b$build"
  SubmitJob $bamid "topmed-$me" '2G' "$0 $*"
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
GetNWDID $bamid
cramfile=$crampath/$nwdid.recab.cram

#======================================================================
#   Copy remapped CRAM from GCE, check flagstat, fix up database
#======================================================================
cramflagstat=`$topmedcmd show $bamid cramflagstat`
if [ "$cramflagstat" = "" ]; then
  Fail "Unable to get cramflagstat for bamid '$bamid'"
fi

stime=`date +%s`

#   Remapped cram could be > one place (arrgh!)  Figure out where it is
inuri=''
for i in $incominguri $bcfuri; do
  p="$i/$nwdid/$nwdid.recab.cram"
  $gsutil stat "$p"
  if [ "$?" = "0" ]; then
    inuri=$i
    break
  fi
done
if [ "$inuri" = "" ]; then
  Fail "Unable to find $nwdid/$nwdid.recab.cram in: $incominguri $bcfuri"
  exit 3
fi

#   Now know where to look for data. Check flagstat
echo "Checking if flagstat is as we expect from $inuri"
$gsutil cp  $inuri/$nwdid/$nwdid.recab.cram.flagstat $crampath
if [ "$?" != "0" ]; then
  Fail "Failed to copy flagstat from GCE: $inuri/$nwdid/$nwdid.recab.cram.flagstat"
fi
#   Get number of interest from flagstat file and check it
n=`grep 'paired in sequencing' $crampath/$nwdid.recab.cram.flagstat | awk '{print $1}'`
if [ "$n" != "$cramflagstat" ]; then
  # Renaming the flagstat file stops pull from happening again
  $gsutil mv $inuri/$nwdid/$nwdid.recab.cram.flagstat $inuri/$nwdid/$nwdid.recab.cram.flagstat.nomatch
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
$gsutil cp $inuri/$nwdid/$nwdid.recab.cram $crampath
if [ "$?" != "0" ]; then
  Fail "Failed to copy file from GCE $inuri/$nwdid/$nwdid.recab.cram $crampath"
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
md5=(`md5sum $cramfile`)
md5=${md5[0]}
if [ "$md5" = "" ]; then
  Fail "Unable to calculate MD5 for remapped '$bamid' [$nwdid] cramfile=$cramfile"
fi
echo "Set checksum and flagstat for b$build file"
SetDB $bamid b${build}cramchecksum $md5
SetDB $bamid b${build}flagstat $cramflagstat

#   Save date of file in database
$topmedcmd setdate $bamid datemapping_b38 $cramfile

#   Clean up data in GCE if data found in incoming.  Move remapped data to bcf bucket
if [ "$inuri" = "$incominguri" ]; then
  $gsutil mv $incominguri/$nwdid/$cramfile.flagstat  $bcfuri/$nwdid/$cramfile.flagstat
  $gsutil mv $incominguri/$nwdid/$cramfile           $bcfuri/$nwdid/$cramfile
  echo "Moved $incominguri/$nwdid to $bcfuri/$nwdid"
  #   Remove any left over cruft in recabs bucket
  echo "Removing $incominguri/$nwdid"
  $gsutil rm -r $incominguri/$nwdid
else
  echo "Data was not found in $incominguri so we leave it where it was"
fi

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Copy of remapped CRAM from GCE to $crampath completed in $etime seconds"

SetDB $bamid state_b${build} 20     # Mark b38 as done
Successful
Log $etime
exit
