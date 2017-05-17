#!/bin/bash
#
#   topmed_verify.sh -submit bamid
#
#	Verify the MD5 checksum for a BAM file
#   Do not specify a QOS for verify so it runs before QPLOT
#
. /usr/cluster/topmed/bin/topmed_actions.inc
me=verify
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me  $bamid
  MyRealHost $bamid 'bam'
  SubmitJob $bamid "$realhost-verify" '2G' "$0 $*"
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
checksum=`$topmedcmd show $bamid checksum`
bamfile=`$topmedpath wherefile $bamid bam`

Started

#   Verify the MD5
stime=`date +%s`
tmpfile=/tmp/$$.md5
echo "$checksum  $bamfile" > $tmpfile
md5sum -c $tmpfile
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM  completed in $etime seconds"
rm -f $tmpfile
if [ "$rc" != "0" ]; then
  Fail "MD5sum failed"
fi

#   Set bamsize again to be sure
sz=`ls -L -l $bamfile | awk '{print $5}'`
SetDB $bamid 'bamsize' $sz

#   Get the paired reads count for this file
stime=`date +%s`
echo "Calculate flagstat"
$topmedflagstat $bamfile $bamid bamflagstat
if [ "$?" != "0" ]; then
  Fail "$topmedflagstat $bamfile failed"
  exit 1
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Calculated bamflagstat in $etime seconds"

#   If necessary, create the index file
extension="${bamfile##*.}"
bai=$bamfile.bai
if [ "$extension" = "cram" ]; then
  bai=$bamfile.crai
fi
if [ ! -f $bai ]; then
  echo "Creating index file '$bai'"
  $samtools index $bamfile 2>&1
  if [ "$?" != "0" ]; then
    #   This might be a trashed reference index for samtools, if so remove it
    a=`grep 'cram_ref_load: Assertion' $console/$bamid-$me.out`
    if [ "$a" != "" ]; then
      rm -rf $HOME/.cache/hts-ref/*/*
      Fail "Unable to create index file for '$bamfile' - removed dirty reference cache. Just restart me."
    fi
    Fail "Unable to create index file for '$bamfile'"
  fi
fi

#   If original file was cram, then some fields are the same for both cram and bam
if [ "$extension" = "cram" ]; then
  a=`$topmedcmd -persist show $bamid bamflagstat`
  SetDB $bamid 'cramflagstat' $a
  a=`$topmedcmd -persist show $bamid checksum`
  SetDB $bamid 'cramchecksum' $a
else
  #   Rename the BAM file
  $topmedrename $bamid $bamfile
  if [ "$?" != "0" ]; then
    Fail "$topmedrename failed"
    exit 1
  fi
fi

Successful
Log $etime
exit
