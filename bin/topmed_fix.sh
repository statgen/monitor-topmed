#!/bin/bash
#
#   topmed_fix.sh -submit bamid|nwdid
#
#	Fix some sort of problem. This script changes all the time
#
. /usr/cluster/topmed/bin/topmed_actions.inc
me=fix
markverb=fix

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  b38=`GetDB $bamid state_b38`
  if [ "$b38" != "20" ]; then
    SetDB $bamid state_fix 0            # We will do this later
    echo "No B38 cram, do not submit job for $bamid"
    exit
  fi
  #MayIRun $me  $bamid
  #h=(`$topmedpath wherepath $bamid b38 | sed -e 's:/: :g'`)
  #realhost=${h[1]}
  CheckRGMap $bamid
  if [ "$?" = "0" ]; then
    SetDB $bamid state_fix 20           # Mark sample as fixed so no rgmap needed
    echo "No need to submit job for $1"
    exit
  fi
  RandomRealHost $bamid
  qoshost=$realhost
  if [ "$qoshost" = "topmed" ]; then
    qoshost=topmed1
  fi
  echo "Submitting $bamid to run on $realhost"
  SubmitJob $bamid "$qoshost-fix" '8G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid|nwdid"
  echo ""
  echo "Fix some sort of problem."
  exit 1
fi

bamid=$1

b38=`GetDB $bamid state_b38`
if [ "$b38" != "20" ]; then
  SetDB $bamid state_fix 0            # We will do this later
  echo "No B38 cram, do not attempt fix for $bamid"
  exit
fi


Started

echo "# $0 $bamid $nwdid  -- redo header for mis-mapped b38 crams"
fixlog=$console/fix.log
CheckRGMap $bamid
if [ "$?" = "0" ]; then
  SetDB $bamid state_fix 20           # Mark sample as fixed so no rgmap needed
  echo "No need to fix $bamid" >> $fixlog
  echo "No need to fix $bamid"
  Successful
  exit
fi

stime=`date +%s`
rgmapmsg=/tmp/$bamid.fixerr
/usr/cluster/topmed/bin/topmed_rgmap.sh $bamid $markverb 2> $rgmapmsg
if [ "$?" != "0" ]; then
  remap=''
  a=`grep 'Strange order' $rgmapmsg`        # Maybe this must be remapped
  if [ "$a" = "" ]; then
    a=`grep 'must be remapped' $rgmapmsg`
  fi
  if [ "$a" != "" ]; then
    remap='Remap this sample'
    SetDB $bamid state_fix 0
    SetDB $bamid state_b38 0
    SetDB $bamid state_gce38copy 0
    SetDB $bamid state_aws38copy 0
    SetDB $bamid state_gce38bcf 0
    SetDB $bamid state_gce38pull 0
    #SetDB $bamid state_gce38push 0
  fi
  e=`cat $rgmapmsg`
  rm -f $rgmapmsg
  echo "Unable to correct rgmap bamid=$bamid $remap: $e" >> $fixlog
  Fail "$e"
fi
rm -f $rgmapmsg

echo "Corrected rgmap for recab for $bamid" >> $fixlog
Successful
etime=`date +%s`
etime=`expr $etime - $stime`
Log $etime
exit

#================= Old runs =================#
exit

. /usr/cluster/topmed/bin/topmed_actions.inc
me=fix
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  MayIRun $me  $bamid
  RandomRealHost $bamid
  SubmitJob $bamid "topmed-fix" '2G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid|nwdid"
  echo ""
  echo "Fix some sort of problem."
  exit 1
fi

Started
bamid=$1
b=b37
echo "# $0 $bamid $nwdid  -- set ${b}cramchecksum for $b remapped cram"

#   Get file to process
file=`$topmedpath wherefile $bamid $b`
if [ "$file" = "" ]; then
  echo "Unable to find file for '$b'"
  exit 4
fi 

#   Verify the MD5
stime=`date +%s`
md5=`md5sum < $file`
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM  completed in $etime seconds"
if [ "$rc" != "0" ]; then
  echo "MD5sum for $bamid build $b failed" >> $console/fix.log
  Fail "MD5sum failed" 
fi
SetDB $bamid ${b}cramchecksum $md5
echo "MD5sum for $bamid build $b successful" >> $console/fix.log
Successful
Log $etime
exit


nwdid=`GetNWDID $1`
bamid=`GetDB $1 bamid`

#   Get file to process
bcffile=`$topmedpath wherefile $bamid bcf`
if [ "$bcffile" = "" ]; then
  echo "Unable to find file for '$b'"
  exit 4
fi
if [ ! -f $bcffile ]; then
  echo "No BCF file found: $bcffile"
  exit 2
fi

vt=/usr/cluster/software/trusty/topmed-year1-freeze3a/master/vt/vt
vtref=/net/topmed/working/mapping/gotcloud/ref/hg38/hs38DH.fa
bcftools=/usr/cluster/bin/bcftools

$bcftools index $bcffile
if [ "$?" != "0" ]; then
  echo "Unable to run BCFTOOLS for bamid '$bamid' [$nwdid] on $bcffile"
  exit 4
fi
echo "Created index BCF file for $bcffile"
exit
#================= Old runs =================#
Started
bamid=$1
b=b37
echo "# $0 $bamid $nwdid  -- set ${b}cramchecksum for $b remapped cram"

#   Get file to process
file=`$topmedpath wherefile $bamid $b`
if [ "$file" = "" ]; then
  echo "Unable to find file for '$b'"
  exit 4
fi 

#   Verify the MD5
stime=`date +%s`
md5=`md5sum < $file`
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM  completed in $etime seconds"
if [ "$rc" != "0" ]; then
  echo "MD5sum for $bamid build $b failed" >> $console/fix.log
  Fail "MD5sum failed" 
fi

SetDB $bamid ${b}cramchecksum $md5

echo "`date   MD5sum for $bamid build $b successful" >> $console/fix.log
Successful
Log $etime
exit
