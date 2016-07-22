#!/bin/bash
#
#   topmed_bai.sh -submit bamid bamfile
#
#	Create BAI index for a BAM file
#
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
samtools=/net/mario/gotcloud/bin/samtools
mem=4G                  # Artificially high so not too many on small nodes
console=/net/topmed/working/topmed-output
slurmp=topmed
qos=topmed-bai
constraint="--constraint eth-10g"
realhost=''

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test bai $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos $constraint --workdir=$console -J $1-bai --output=$console/$1-bai.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos $constraint --workdir=$console -J $1-bai --output=$console/$1-bai.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 baid submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` bai ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
  exit
fi

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid bamfile"
  echo ""
  echo "Create index for a cram/bam file and update database"
  exit 1
fi
bamid=$1
bamfile=$2

#   Is this a cram or bam
extension="${bamfile##*.}"

#   Mark this as started
$topmedcmd mark $bamid baid started
d=`date +%Y/%m/%d`
s=`hostname`
echo "#========= $d host=$s $SLURM_JOB_ID $0 bamid=$bamid bamfile=$bamfile ========="
stime=`date +%s`
bai=$bamfile.bai
if [ "$extension" = "cram" ]; then
  bai=$bamfile.crai
fi
if [ -f $bai ]; then
  echo "Using existing index file '$bai'"
  $topmedcmd mark $bamid baid completed
  exit 0
fi

echo "Creating index file '$bai'"
$samtools index $bamfile 2>&1
if [ "$?" != "0" ]; then
  echo "Unable to create index file"
  $topmedcmd mark $bamid baid failed
  exit 2
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Created index '$bai' (at second $etime)"
$topmedcmd mark $bamid baid completed
echo `date` bai $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids

