#!/bin/bash
#
#   topmed_bai.sh -submit bamid bamfile
#
#	Create BAI index for a BAM file
#
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
gcbin=/net/mario/gotcloud/bin
mem=${TOPMED_MEM:-1G}
console=/net/topmed/working/topmed-output
slurmp=${TOPMED_PARTITION:-topmed-incoming}
slurmqos=${TOPMED_QOS:-topmed-bai}

if [ "$1" = "-submit" ]; then
  shift
  homehost=`echo $2 | cut -d / -f 3`    # Should be topmed/topmed2
  if [ "$homehost" != "" ]; then
    console=/net/$homehost/working/topmed-output
    slurmp="$homehost-incoming"
  fi

  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-bai --output=$console/$1-bai.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-bai --output=$console/$1-bai.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 baid submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobid
    $topmedcmd set $1 jobidbai ${l[3]}
  fi
  exit
fi

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid bamfile"
  echo ""
  echo "Create BAI for a bam file and update database"
  exit 1
fi
bamid=$1
bamfile=$2

#   Mark this as started
$topmedcmd mark $bamid baid started
d=`date +%Y/%m/%d`
s=`hostname`
echo "#========= $d host=$s $SLURM_JOB_ID $0 bamid=$bamid bamfile=$bamfile ========="
stime=`date +%s`
bai=$bamfile.bai
if [ -f $bai ]; then
  echo "Using existing BAI file '$bai'"
  $topmedcmd mark $bamid baid completed
  exit 0
fi

echo "Creating BAI file '$bai'"
$gcbin/samtools index $bamfile 2>&1
if [ "$?" != "0" ]; then
  echo "Unable to create BAI file"
  $topmedcmd mark $bamid baid failed
  exit 2
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Created BAI '$bai' (at second $etime)"
$topmedcmd mark $bamid baid completed

