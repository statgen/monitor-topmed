#!/bin/bash
#
#   topmed_fix.sh -submit bamid|nwdid
#
#	Fix some sort of problem. This script changes all the time
#
. /usr/cluster/topmed/bin/topmed_actions.inc
me=fix
markverb=$me
export IGNORE_PERMIT=1         # Ignore MayIRun stuff

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='4:00:00'
  SubmitJob $bamid "topmed" '4G' "$0 $*"
  exit
fi

#
#   Rerun topmedcheck on all the centers
#
/usr/cluster/topmed


#
#   Fetch completed recab and put back in localstore as normal
#
n=$1
stime=`date +%s`
b38=`$topmedpath wherefile $n b38`
bcf=`$topmedpath wherefile $n bcf`
b38path=`$topmedpath wherepath $n b38`
bcfpath=`$topmedpath wherepath $n bcf`
d=/tmp/nwd.$n
mkdir -p $d
cd $d || exit 3
gsutil cp -r gs://topmed-bcf/$n/\* . || exit 4
mkdir -p $b38path
mkdir -p $bcfpath
cp -p $n.recab.cram $b38 || exit 4
chmod 0640 $b38
cp -p $n.recab.cram.crai $b38.crai
chmod 0640 $b38.crai
cp -p $n.bcf $bcf || exit 4
chmod 0640 $bcf
cp -p $n.bcf.csi $bcf.csi
chmod 0640 $bcf.csi
rm -rf $d

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Fetched recab for '$n' and updated local files in $etime seconds"
SetDB $bamid state_fix 0        # Reset this state
exit

#
#   Calculate the cram flagstat again
#
#   Get the paired reads count for this file
bamid=$1
echo "Calculate flagstat for $bamid"
bf=`GetDB $bamid bamflagstat`
cf=`GetDB $bamid cramflagstat`
if [ "$bf" != "0" -a "$bf" = "$cf" ]; then
  echo "Flagstat $bf already set for bam and cram"
  exit 0
fi
cramfile=`$topmedpath wherefile $bamid cram`
flagstat=`CalcFlagstat $bamid $cramfile`
if [ "$flagstat" -lt "250000000" ]; then
  Fail "Flagstat $bamflagstat seems too small, try again cause samtools is a flake"
fi
SetDB $bamid bamflagstat $flagstat
SetDB $bamid cramflagstat $flagstat
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Calculated flagstat in $etime seconds"
SetDB $bamid state_fix 0        # Reset this state
exit
