#!/bin/bash
#
#   topmed_qplot.sh [-submit [-mem NG]] bamid bamfile
#
#	Run QPLOT on a BAM file
#
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
console=/net/topmed/working/topmed-output
gcbin=/net/mario/gotcloud/bin
gcref=/net/mario/nodeDataMaster/local/ref/gotcloud.ref
topoutdir=/net/topmed/incoming/qc.results
mem=8G
if [ "$TOPMED_MEMORY" != "" ]; then mem=$TOPMED_MEMORY; fi
qos=qplot
if [ "$TOPMED_QOS" != "" ]; then qos=$TOPMED_QOS; fi

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test qplot $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #   Figure where to submit this to run - should be local
  l=(`$topmedcmd where $1`)     # Get bampath backuppath bamname realhost realhostindex
  realhost="${l[3]}"
  if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
  realhostindex="${l[4]}"
  slurmp="$realhost-incoming"   # Sometimes we think this should be nomosix
  slurmqos="$realhost-$qos"

  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-qplot --output=$console/$1-qplot.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-qplot --output=$console/$1-qplot.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 qploted submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobid
    echo `date` qplot ${l[3]} >> $console/$1.jobids
  fi
  exit
fi

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit [-mem nG]] bamid bamfile"
  echo ""
  echo "Run qplot on a bam file and update database"
  exit 1
fi
bamid=$1
bamfile=$2

homehost=`echo $bamfile | cut -d / -f 3`    # Should be topmed or topmned2
if [ "$homehost" != "" ]; then
  console=/net/$homehost/working/topmed-output
fi

#   Mark this as started
$topmedcmd mark $bamid qploted started
d=`date +%Y/%m/%d`
stime=`date +%s`
s=`hostname`
echo "#========= $d host=$s $SLURM_JOB_ID $0 bamid=$bamid bamfile=$bamfile ========="
bai=$bamfile.bai
if [ ! -f $bai ]; then
  echo "BAI '$bai' does not exist"
  $topmedcmd mark $bamid qploted failed
  exit 3
  echo "Creating BAI file '$bai'"
  $gcbin/samtools index $bamfile 2>&1
  if [ "$?" != "0" ]; then
    echo "Unable to create BAI file"
    $topmedcmd mark $bamid qploted failed
    exit 2
  fi
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "Created BAI '$bai' (at second $etime)"
fi

#   Create output directory and CD there
o=`dirname $bamfile`
o=`echo $o | sed -e s:/net/topmed/incoming/topmed/::`
outdir=$topoutdir/$o
mkdir -p $outdir
cd $outdir
if [ "$?" != "0" ]; then
  echo "Unable to CD to qplot output directory for '$bamid' - $outdir"
  $topmedcmd mark $bamid qploted failed
  exit 3
fi

#   Run qplot, output written to current working directory
basebam=`basename $bamfile .bam`
$gcbin/qplot --reference  $gcref/hs37d5.fa --dbsnp $gcref/dbsnp_142.b37.vcf.gz \
  --label $basebam --stats $basebam.qp.stats --Rcode $basebam.qp.R $bamfile 2>&1
if [ "$?" != "0" ]; then
  echo "QPLOT failed for '$bamfile'"
  $topmedcmd mark $bamid qploted failed
  rm -f $basebam.*
  exit 4
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "QPLOT on '$bamfile' successful (at second $etime)"

#   Run verifybamid, output written to current working directory
$gcbin/verifyBamID --bam  $bamfile --vcf $gcref/hapmap_3.3.b37.sites.vcf.gz	\
  --site --free-full --chip-none --ignoreRG --precise --maxDepth 80 \
  --grid 0.02	--out  $basebam.vb 2>&1
if [ "$?" != "0" ]; then
  echo "VerifyBAMID failed for '$bamid'"
  $topmedcmd mark $bamid qploted failed
  rm -f $basebam.*
  exit 5
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "VerifyBAMID on '$bamfile' successful (at second $etime)"
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Command completed in $etime seconds. Created files:"
ls -la $basebam.*
$topmedcmd mark $bamid qploted completed
exit 0
