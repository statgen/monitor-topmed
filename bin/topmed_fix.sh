#!/bin/bash
#
#   topmed_fix.sh -submit bamid|nwdid
#   Another possibility is:  topmed_fix.sh  realhost  file-of-bamids
#
#	Fix some sort of problem. This script changes all the time
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
me=fix
me=$1
markverb=$me

#   Special case - run a command against a set of bamids
if [ "$1" = "-submit" ]; then
  shift
  realhost=$1       # Force to run here
  shift
  timeout='24:00:00'
  s=`date +%s`
  me="$me-$s"
  SubmitJob $realhost $PROJECT '8G' "$0 $*"
  echo "Ignore error from topmedcmd about invalid bamid"
  exit
fi

stime=`date +%s`
n=0
e=0
k=0
#   e.g. topmed_fix.sh  [-submit] topmed3 ~/coldline/x00
for bamid in `cat $1`; do
  nwdid=`GetNWDID $bamid`
  studyname=`GetDB $bamid studyname`
  k=`expr $k + 1`
  echo "#========== `date` Running sample $b $nwdid $studyname from $1 ========#"
  st=`date +%s`
  rc=0
  for x in recal.cram recal.cram.crai recal.cram.flagstat recal.cram.md5; do
    y=`echo ${studyname,,} | grep copd`
    if [ "$y" != "" ]; then         # Studyname is not always part of the path  POS
      studyname=copd
    fi
    $gsutil rewrite -s coldline gs://topmed-irc-working/remapping/b37/${studyname,,}/$nwdid.$x
    rc=`expr $rc + $?`
  done
  if [ "$rc" = "0" ]; then
    SetDB $bamid state_fix 20
    n=`expr $n + 1`
  else
    SetDB $bamid state_fix 99
    e=`expr $e + 1`
  fi
  et=`date +%s`
  et=`expr $et - $st`
  echo "[$k] Coldline attempted for $nwdid $bamid in $et seconds"
done

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Complete changing $n samples to coldline storage in $etime seconds - $e failed"
exit


exit
#
#   Standard submit of the job
#
#   export IGNORE_PERMIT=1         # Ignore MayIRun stuff
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
#   Fetch completed recab and put back in localstore as normal
#
n=$1
stime=`date +%s`
b38=`$topmedpath wherefile $n b38`
bcf=`$topmedpath wherefile $n bcf`
b38path=`$topmedpath wherepath $n b38`
bcfpath=`$topmedpath wherepath $n bcf`
d=/run/shm/nwd.$n
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
