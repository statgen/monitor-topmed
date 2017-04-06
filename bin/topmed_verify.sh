#!/bin/bash
#
#   topmed_verify.sh -submit bamid
#
#	Verify the MD5 checksum for a BAM file
#   Do not specify a QOS for verify so it runs before QPLOT
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=verify
mem=2G
markverb=$me
qos="--qos=topmed-$me"
realhost=''

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedpermit permit test $me $1
  if [ "$?" = "0" ]; then
    echo "$me $1 not permitted" | tee $console/$1-$me.out
    exit 4
  fi 

  # Run this on node where bam lives
  h=`$topmedpath whathost $1 bam`
  if [ "$h" != "" ]; then
    realhost="--nodelist=$h"
    #qos="--qos=$h-$me"
  fi

  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $qos $realhost -J $1-$me --output=$console/$1-$me.out $0 $*`)
  if [ "$?" != "0" ]; then
    $topmedcmd -emsg "Failed to submit command to SLURM - $l" mark $1 $markverb failed
    echo "Failed to submit command to SLURM - $l" > $console/$1-$me.out
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $qos $realhost -J $1-$me --output=$console/$1-$me.out $0 $*" >> $console/$1-$me.out
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` $me ${l[3]} $slurp $slurmqos $mem >> $console/$1.jobids
  fi
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

#   Is this a cram or bam
extension="${bamfile##*.}"

tmpfile=/tmp/$$.md5
echo "$checksum  $bamfile" > $tmpfile

#   Mark this as started
$topmedcmd mark $bamid $markverb started
d=`date +%Y/%m/%d`
s=`hostname`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid checksum=$checksum bamfile=$bamfile ========="

#   Verify the MD5
stime=`date +%s`
md5sum -c $tmpfile
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM  completed in $etime seconds"
rm -f $tmpfile
if [ "$rc" != "0" ]; then
  $topmedcmd -persist -emsg "Md5sum failed" mark $bamid $markverb failed
  exit 1
fi

#   Set bamsize again to be sure
sz=`ls -L -l $bamfile | awk '{print $5}'`
$topmedcmd -persist set $bamid bamsize $sz

#   Get the paired reads count for this file
stime=`date +%s`
echo "Calculate flagstat"
$topmedflagstat $bamfile $bamid bamflagstat
if [ "$?" != "0" ]; then
  $topmedcmd -persist -emsg "$topmedflagstat $bamfile failed" mark $bamid $markverb failed
  exit 1
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Calculated bamflagstat in $etime seconds"

chmod 0444 $bamfile     # This might fail, but that's OK

#   If original file was cram, then some fields are the same for both cram and bam
if [ "$extension" = "cram" ]; then
  a=`$topmedcmd -persist show $bamid bamflagstat`
  $topmedcmd -persist set $bamid cramflagstat $a
  a=`$topmedcmd -persist show $bamid checksum`
  $topmedcmd -persist set $bamid cramchecksum $a
else
  #   Rename the BAM file
  $topmedrename $bamid $bamfile
  if [ "$?" != "0" ]; then
    $topmedcmd -persist -emsg "$topmedrename failed" mark $bamid $markverb failed
    exit 1
  fi
fi

$topmedcmd -persist mark $bamid $markverb completed
echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit
