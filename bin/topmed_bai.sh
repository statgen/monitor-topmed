#!/bin/bash
#
#   topmed_bai.sh -submit bamid
#
#	Create BAI index for a BAM file
#
. /usr/cluster/topmed/bin/topmed_actions.inc

mem=4G
me=bai
markverb="${me}d"
slurmp=topmed-working
qos="--qos=topmed-$me"
constraint="--constraint eth-10g"
realhost=''

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedpermit permit test $me $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  # Keep this off 9 and 10
  n=`expr $1 % 5`
  if [ "$n" = "0" ]; then
    n='5';
  fi
  if [ "$n" = "1" ]; then
    n='';
  fi
  nodelist="--nodelist=topmed$n"

  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $qos $constraint $nodelist  --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $qos $nodelist $constraint --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` $me ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Create index for a cram/bam file and update database"
  exit 1
fi
bamid=$1
bamfile=`$topmedpath wherefile $bamid bam`

#   Is this a cram or bam
extension="${bamfile##*.}"

#   Mark this as started
$topmedcmd mark $bamid $markverb started
d=`date +%Y/%m/%d`
s=`hostname`
echo "#========= $d host=$s $SLURM_JOB_ID $0 bamid=$bamid bamfile=$bamfile ========="
stime=`date +%s`
bai=$bamfile.bai
if [ "$extension" = "cram" ]; then
  bai=$bamfile.crai
fi
if [ -f $bai ]; then
  chmod 0444 $bai
  echo "Using existing index file '$bai'"
  $topmedcmd -persist mark $bamid $markverb completed
  exit 0
fi

echo "Creating index file '$bai'"
$samtools index $bamfile 2>&1
if [ "$?" != "0" ]; then
  echo "Unable to create index file"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Created index '$bai' (at second $etime)"
$topmedcmd -persist mark $bamid $markverb completed
echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids

