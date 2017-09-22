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
  #MayIRun $me  $bamid
  #h=(`$topmedpath wherepath $bamid b38 | sed -e 's:/: :g'`)
  #realhost=${h[1]}
  #SubmitJob $bamid "$realhost-fix" '8G' "$0 $*"
  RandomRealHost $bamid
  SubmitJob $bamid "topmed-fix" '8G' "$0 $*"
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

Started

echo "# $0 $bamid $nwdid  -- redo header for mis-mapped b38 crams"
fixlog=$console/fix.log
stime=`date +%s`
/usr/cluster/topmed/bin/topmed_rgmap.sh $bamid $markverb
if [ "$?" != "0" ]; then
  echo "Unable to correct rgmap bamid=$bamid" >> $fixlog
  Fail "Unable to correct rgmap" 
fi

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
