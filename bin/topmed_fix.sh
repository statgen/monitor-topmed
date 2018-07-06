#!/bin/bash
#
#   topmed_fix.sh -submit bamid|nwdid
#
#	Fix some sort of problem. This script changes all the time
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
me=fix
markverb=$me

#   Special case - run a command against a set of bamids
#   topmed_fix.sh host cmd file-of-bamids
if [ "$1" = "-submit" ]; then
  shift
  realhost=$1       # Force to run here
  shift
  timeout='24:00:00'
  SubmitJob $realhost $PROJECT '4G' "$0 $*"
  exit
fi

#   Command here is:  topmed_fix.sh  cmd  file-of-bamids
#   e.f.  topmed_fix.sh  [-submit] /usr/cluster/topmed/bin/gcebackup.sh ~/set1.bamids
#
#   ./xsql.sh 'SELECT count(*) from bamfiles where state_gcebackup!=20'
#   ./xsql.sh 'SELECT bamid from bamfiles where state_gcebackup!=20 order by rand()' >setofbams
#   head -9000 setofbams | split -d -l 1000    # Makes xdd files of 1000 each
#   bin/topmed_fix.sh -submit topmed2 /usr/cluster/topmed/bin/topmed_gcebackup.sh ~/x02
#   bin/topmed_fix.sh -submit topmed3 /usr/cluster/topmed/bin/topmed_gcebackup.sh ~/x03
#   bin/topmed_fix.sh -submit topmed4 /usr/cluster/topmed/bin/topmed_gcebackup.sh ~/x04
#   bin/topmed_fix.sh -submit topmed5 /usr/cluster/topmed/bin/topmed_gcebackup.sh ~/x05
#   bin/topmed_fix.sh -submit topmed6 /usr/cluster/topmed/bin/topmed_gcebackup.sh ~/x06
#   bin/topmed_fix.sh -submit topmed7 /usr/cluster/topmed/bin/topmed_gcebackup.sh ~/x07
#   bin/topmed_fix.sh -submit topmed9 /usr/cluster/topmed/bin/topmed_gcebackup.sh ~/x08
#   bin/topmed_fix.sh -submit topmed10 /usr/cluster/topmed/bin/topmed_gcebackup.sh ~/x01
#   bin/topmed_fix.sh -submit topmed  /usr/cluster/topmed/bin/topmed_gcebackup.sh ~/x00
n=`basename $1`
for b in `cat $2`; do
  echo "#========== `date` Running $1 on sample $b ========#"
  $1 $b > $console/$b-$n.out 2>&1
done
exit


exit
#
#   Standard submit of the job
#
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
