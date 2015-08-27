#!/bin/bash
#
#   topmed_report.sh -submit bamid bamfile
#
#	Collect QPLTOT statistics for a BAM file
#
echo "This is dead code. Irina is doing this herself using her own methods   June 2015"
exit 6


topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
console=/home/topmed/output
topoutdir=/incoming/qc.results
topmedreport=/usr/cluster/monitor/bin/topmed_reportstats.pl

if [ "$1" = "-submit" ]; then
  shift
  sbatch -p topmed-incoming --qos=topmed-incoming-limit -J $1-report --output=$console/$1-report.out $0 $*
  if [ "$rc" = "0" ]; then
    echo "Failed to submit command to SLURM"
    exit 1
  fi
  $topmedcmd mark $1 report submitted
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me   bamid bamfile"
  echo ""
  echo "Collect output from qplot for a report database"
  exit 1
fi
bamid=$1
bamfile=$2

#   Mark this as started
$topmedcmd mark $bamid reported started || exit $?
d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid bamfile=$bamfile ========="
stime=`date +%s`

#   Calculate output directory and CD there
o=`dirname $bamfile`
o=`echo $o | sed -e s:/net/topmed/incoming/topmed/::`
outdir=$topoutdir/$o
mkdir -p $outdir
cd $outdir
if [ "$?" != "0" ]; then
  echo "Unable to CD to qplot output directory for bam '$bamid' - $outdir"
  $topmedcmd mark $bamid reported failed
  exit 2
fi
basebam=`basename $bamfile .bam`
qpstats=$basebam.qp.stats
if [ ! -f "$qpstats" ]; then
  echo "Unable to find qplot stats output for bam '$bamid' - $outdir/$qpstats"
  $topmedcmd mark $bamid reported failed
  exit 3
fi

#   Run program to collect qplot stats
$topmedreport $bamid $qpstats
if [ "$?" != "0" ]; then
  echo "Unable to collect qplot stats for bam '$bamid' - $outdir/$qpstats"
  $topmedcmd mark $bamid reported failed
  exit 2
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "REPORT on '$bamfile' successful (at second $etime)"

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Command completed in $etime seconds. Created files:"
$topmedcmd mark $bamid reported completed
exit 0
