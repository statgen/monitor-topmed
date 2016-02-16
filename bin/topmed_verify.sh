#!/bin/bash
#
#   topmed_verify.sh -submit bamid checksum bamfile
#
#	Verify the MD5 checksum for a BAM file
#   Do not specify a QOS for verify so it runs before QPLOT
#
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
topmedrename=/usr/cluster/monitor/bin/topmedrename.pl
console=/net/topmed/working/topmed-output
mem=1G
if [ "$TOPMED_MEMORY" != "" ]; then mem=$TOPMED_MEMORY; fi
realhost=topmed
#if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
markverb=md5verified
qos=verify
if [ "$TOPMED_QOS" != "" ]; then qos=$TOPMED_QOS; fi

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test verify $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #   Figure where to submit this to run - should be local
  l=(`$topmedcmd where $1 bam`)         # Get pathofbam and host for bam
  h="${l[1]}"
  if [ "$h" != "" ]; then realhost=$h; fi
  if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
  slurmp="$realhost-incoming"
  slurmqos="$realhost-$qos"

  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos -J $1-verify --output=$console/$1-verify.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos -J $1-verify --output=$console/$1-verify.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` verify ${l[3]} $slurp $slurmqos $mem >> $console/$1.jobids
  fi
  exit
fi

if [ "$3" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid checksum bamfile"
  echo ""
  echo "Verify checksum for a bam file and update database"
  exit 1
fi
bamid=$1
checksum=$2
bamfile=$3

tmpfile=/tmp/$$.md5
echo "$checksum  $bamfile" > $tmpfile

#   Mark this as started
$topmedcmd mark $bamid $markverb started
d=`date +%Y/%m/%d`
s=`hostname`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid checksum=$checksum bamfile=$bamfile ========="
stime=`date +%s`
md5sum -c $tmpfile
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Command completed in $etime seconds"
rm -f $tmpfile
if [ "$rc" != "0" ]; then
  $topmedcmd mark $bamid $markverb failed
  exit 1
fi
$topmedcmd mark $bamid $markverb completed
echo `date` verify $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids

#   Set bamsize again to be sure
sz=`ls -L -l $bamfile | awk '{print $5}'`
$topmedcmd set $bamid bamsize $sz

#   Rename the BAM file and change the MD5 entry
chmod 0444 $bamfile
$topmedrename $bamid $bamfile
if [ "$?" != "0" ]; then
  $topmedcmd mark $bamid $markverb failed
  exit 1
fi
exit
