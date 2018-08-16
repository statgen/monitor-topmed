#!/bin/bash
#
#   topmed_verify.sh -submit bamid
#
#	Verify the MD5 checksum for a BAM file
#   Do not specify a QOS for verify so it runs before QPLOT
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
me=verify
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='4:00:00'
  SubmitJob $bamid $PROJECT '4G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Verify checksum for a bam file and update database"
  exit 1
fi
bamid=$1
checksum=`GetDB $bamid checksum`
bamfile=`$topmedpath wherefile $bamid bam`

Started

#   Verify the MD5 - special case, if md5 is all zeroes, we use what we calculate
stime=`date +%s`
tmpfile=/run/shm/$$.md5
if [ "$checksum" = "00000000000000000000000000000000" ]; then
  sum=`md5sum $bamfile`
  if [ "$sum" = "" ]; then
    Fail "MD5sum calculation failed"
  fi
  SetDB $bamid checksum ${sum[0]}
  echo "Checksum forced to ${sum[0]}"
else  
  echo "$checksum  $bamfile" > $tmpfile
  md5sum -c $tmpfile
  if [ "$?" != "0" ]; then
    Fail "MD5sum failed: `cat $tmpfile`"
  fi
  rm -f $tmpfile
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM  completed in $etime seconds"

#   Set bamsize again to be sure
sz=`ls -L -l $bamfile | awk '{print $5}'`
SetDB $bamid 'bamsize' $sz

#   Get the paired reads count for this file
stime=`date +%s`
echo "Calculate flagstat"
bamflagstat=`CalcFlagstat $bamid $bamfile`
if [ "$bamflagstat" -lt "250000000" ]; then
  Fail "Flagstat $bamflagstat seems too small, try again cause samtools is a flake"
fi
SetDB $bamid bamflagstat $bamflagstat
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Calculated bamflagstat in $etime seconds"

#   Create the index file as necessary
CreateIndex $bamid $bamfile

#   If original file was cram, then some fields are the same for both cram and bam
extension="${bamfile##*.}"
if [ "$extension" = "cram" ]; then
  a=`GetDB $bamid bamflagstat`
  SetDB $bamid 'cramflagstat' $a
  a=`GetDB $bamid checksum`
  SetDB $bamid 'cramchecksum' $a
else
  #   Rename the BAM file
  $topmedrename $bamid $bamfile
  if [ "$?" != "0" ]; then
    Fail "$topmedrename failed"
    exit 1
  fi
fi

rm -f $bamfile.tmp          # Some pgm creates this crap file
Successful
Log $etime
exit
