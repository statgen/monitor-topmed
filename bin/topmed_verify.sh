#!/bin/bash
#
#   topmed_verify.sh -submit bamid checksum bamfile
#
#	Verify the MD5 checksum for a BAM file
#   Do not specify a QOS for verify so it runs before QPLOT
#
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
console=/net/topmed/working/topmed-output
mem=${TOPMED_MEM:-1G}

if [ "$1" = "-submit" ]; then
  shift

  #   Figure where to submit this to run - should be local
  l=(`$topmedcmd where $1`)             # Get bampath backuppath bamname realhost realhostindex
  realhost="${l[3]}"
  realhostindex="${l[4]}"
  slurmp="$realhost-incoming"
  slurmqos="$realhost-verify"

  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos -J $1-verify --output=$console/$1-verify.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos -J $1-verify --output=$console/$1-verify.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 md5verified submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobid
    $topmedcmd set $1 jobidmd5ver ${l[3]}
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
$topmedcmd mark $bamid md5verified started
d=`date +%Y/%m/%d`
s=`hostname`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 checksum=$checksum bamfile=$bamfile ========="
stime=`date +%s`
md5sum -c $tmpfile
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Command completed in $etime seconds"
rm -f $tmpfile
if [ "$rc" != "0" ]; then
  $topmedcmd mark $bamid md5verified failed
  exit 0
fi
$topmedcmd mark $bamid md5verified completed
#   Set bamsize again to be sure
sz=`ls -l $bamfile | awk '{print $5}'`
$topmedcmd set $bamid bamsize $sz
exit $rc
