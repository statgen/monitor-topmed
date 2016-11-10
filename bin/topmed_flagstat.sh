#!/bin/bash
#
#   topmed_flagstat.sh amfile [[bamid colname] outfile]
#
#	Calculate the flagstat paired in sequencing number for a file
#   and save it in a database column
#   This can handle bam or cram input files
#   If bamid is set, this can handle illumina files
#   If colname is provided, the flagstat value will be saved in this column of the database
#   If outfile is provided the flagstat results will be copied to this file
samtools=/usr/cluster/bin/samtools
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
ref=/net/mario/gotcloud/gotcloud.ref/hs37d5.fa
illuminaref=/net/topmed/incoming/study.reference/study.reference/illumina.hg19.fa

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me amfile [[bamid colname] outfile]"
  echo ""
  echo "Calculate the flagstat paired in sequencing number for a file"
  exit 1
fi
amfile=$1
bamid=$2
colname=$3
myoutfile=$4

outfile=`basename $amfile`              # Calculate temp file for flagstat results
outfile="/tmp/$outfile.tmp"

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
  echo "$me - Command failed: $samtools flagstat $amfile. Results in $outfile"
  exit 2
fi

#   Outfile has results of flagstat call, get paired reads
a=(`grep 'paired in sequencing' $outfile`)
if [ "${a[0]}" = "" ]; then
  echo "$me - Unable to get reads of paired in sequencing from '$output' for '$amfile'"
  exit 3
fi

#   Maybe put this in the database
if [ "$colname" != "" ]; then
  $topmedcmd  set $bamid $colname ${a[0]}
fi
echo "$me - Flagstat value for '$colname' = ${a[0]}"

if [ "$myoutfile" != "" ]; then
  cp $outfile $myoutfile
fi
rm -f $outfile

