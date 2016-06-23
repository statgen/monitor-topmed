#!/bin/bash
#
#	Calculate the flagstat paired in sequencing number for a file
#   and save it in a database column
#   This can handle bam or cram input files
#   If bamid is set, this can handle illumina files
#   If verb is set, this will mark a step for bamid as failed
#   If colname is set, the flagstat value will be saved in this column of the database
samtools=/net/mario/gotcloud/bin/samtools
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
ref=/net/mario/gotcloud/gotcloud.ref/hs37d5.fa
illuminaref=/net/topmed/incoming/study.reference/study.reference/illumina.hg19.fa

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me amfile [[[bamid] verb] colname]"
  echo ""
  echo "Calculate the flagstat paired in sequencing number for a file"
  exit 1
fi
amfile=$1
bamid=$2
verb=$3
colname=$4

outfile=`basename $amfile`
outfile="$outfile.tmp"

if [ "$bamid" != "" ]; then
  center=`$topmedcmd show $bamid center`
  if [ "$center" = "illumina" ]; then    # Special hack for this center
    ref=$illuminaref
  fi
fi

#   Get flagstat values for bam or cram
a=`echo $amfile | grep .cram`
if [ "$a" = "" ]; then
  $samtools flagstat  $amfile > $outfile
  rc=$?
else
  $samtools view  -u -T  $ref $amfile | $samtools flagstat  - >  $outfile
  rc=$?
fi
if [ "$rc" != "0" ]; then
  echo "Command failed: $samtools flagstat $amfile"
  if [ "$verb" != "" ]; then
    $topmedcmd mark $bamid $verb failed
  fi
  exit 2
fi

#   outfile has results of flagstat call, get paired reads
a=(`grep 'paired in sequencing' $outfile`)
if [ "${a[0]}" = "" ]; then
  echo "Unable to get reads of paired in sequencing from '$output' for '$amfile'"
  if [ "$verb" != "" ]; then
    $topmedcmd mark $bamid $verb failed
  fi
  exit 3
fi

#   Maybe put this in the database
if [ "$colname" != "" ]; then
  $topmedcmd  set $bamid $colname ${a[0]}
else
  echo "Flagstat value for '$colname' = ${a[0]}"
fi
rm -f $outfile
