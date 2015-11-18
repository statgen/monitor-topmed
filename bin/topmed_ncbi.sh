#!/bin/bash
#
#   topmed_ncbi.sh -submit bamid bamfile
#
#	Send the proper set of files to NCBI
#
bindir=/usr/cluster/bin
samtools=/net/mario/gotcloud/bin/samtools
ref=/net/mario/gotcloud/gotcloud.ref/hs37d5.fa
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
ncbidir=/net/topmed/working/cp2ncbi
mem=${TOPMED_MEM:-8G}
console=/net/topmed/working/topmed-output

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test ncbi $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #   Figure where to submit this to run - should be local
  l=(`$topmedcmd where $1`)     # Get bampath backuppath bamname realhost realhostindex
  realhost="${l[3]}"
  realhostindex="${l[4]}"
  slurmp="$realhost-incoming"
  slurmqos="$realhost-ncbi"

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-ncbi --output=$console/$1-ncbi.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-ncbi --output=$console/$1-ncbi.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 cped2ncbi submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobid
    $topmedcmd set $1 jobidcp2ncbi ${l[3]}
  fi
  exit
fi

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid bamfile"
  echo ""
  echo "Copy the BAM and CRAM files to NCBI"
  exit 1
fi
bamid=$1
bamfile=$2

#   This is not implemented for the moment
$topmedcmd mark $bamid cped2ncbi completed
$topmedcmd mark $bamid cped2ncbi delivered
#$topmedcmd mark $bamid cped2ncbi failed
exit 2

#========================== dead code ==========================
#   Now calc destination directory
l=(`$topmedcmd where $bamid`)           # Get bampath backuppath bamname
bampath="${l[0]}"
backupdir="${l[1]}"
bamname="${l[2]}"

if [ "$bampath/$bamname" != "$bamfile" ]; then
  echo ""
  echo ""
  echo "This sure looks wrong - the input is not where I think it should be. Continuing anyway"
  echo "  input    bamfile=$bamfile"
  echo "  expected bamfile=$bampath/$bamname"
  echo ""
  echo ""
fi

d=`date +%Y/%m/%d`
mkdir -p $backupdir
cd $backupdir
if [ "$?" != "0" ]; then
  echo "Unable to CD $backupdir"
  $topmedcmd mark $bamid cped2ncbi failed
  exit 2
fi
s=`hostname`
echo "#========= '$d' host=$s $SLURM_JOB_ID squeezed=$squeezed $0 bamid=$bamid bamfile=$bamfile backupdir=$backupdir ========="

#   Mark this as started
$topmedcmd mark $bamid cped2ncbi started
stime=`date +%s`
#   Try to force everything as readonly
chmod 555 $bamfile 2> /dev/null

#   Run 'bam squeeze' on .bam file with result written as .cram
#   Checking with 'samtools flagstat' is quick, zero means success.
#   Running 'samtools index' on a .cram file does NOT require 
#   an explicit genome reference sequence, but 'samtools view' 
#   does.  (And, 'samtools flagstat' does accept .sam input.)
#
#   All this code used to be just a simple
#       cp -p $bamfile $backupdir
#
chkname=`basename $bamfile .bam`
nwdid=`$samtools view -H $bamfile | grep '^@RG' | grep -o 'SM:\S*' | sort -u | cut -d \: -f 2`
if [ "$nwdid" = "" ]; then
  echo "Unable to extract NWDID from header of '$bamfile'"
  $topmedcmd mark $bamid cped2ncbi failed
  exit 2
fi
newname="$nwdid.src.cram"

now=`date +%s`
$samtools flagstat  $bamfile  >  ${chkname}.init.stat
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools flagstat  $bamfile"
  $topmedcmd mark $bamid cped2ncbi failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "samtools flagstat completed in $s seconds"

now=`date +%s`
if [ "$squeezed" = "n" ]; then
  $bindir/bam squeeze  --in $bamfile  --out - --rmTags "BI:Z;BD:Z;PG:Z"  --keepDups --binMid  --binQualS  2,3,10,20,25,30,35,40,50 | $samtools view  -C -T  $ref  -  >  $newname
else 
  $samtools view  -C -T $ref  $bamfile  >  $newname
fi
if [ "$?" != "0" ]; then
  echo "Command failed: $bindir/bam squeeze  --in $bamfile ..."
  $topmedcmd mark $bamid cped2ncbi failed
  exit 1
fi
s=`date +%s`; s=`expr $s - $now`; echo "squeeze completed in $s seconds"

now=`date +%s`
$samtools index $newname
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools index $newname"
  $topmedcmd mark $bamid cped2ncbi failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "samtools index completed in $s seconds"

now=`date +%s`
$samtools view  -u -T  $ref $newname | $samtools flagstat  - >  ${chkname}.cram.stat
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools view  -h -T  $ref $newname ..."
  $topmedcmd mark $bamid cped2ncbi failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "samtools view completed in $s seconds"

diff=`diff  ${chkname}.init.stat  ${chkname}.cram.stat | wc -l`
if [ "$diff" != "0" ]; then
  echo "Stat for backup CRAM file differs from that for original BAM"
  diff  ${chkname}.init.stat  ${chkname}.cram.stat
  $topmedcmd mark $bamid cped2ncbi failed
  exit 2
fi
echo "Stat for CRAM file matches that of original"
rm -f ${chkname}.init.stat  ${chkname}.cram.stat ${chkname}.bam ${chkname}.bam.md5  ${chkname}.bai ${chkname}.bai.md5

echo "Calculating new MD5"
now=`date +%s`
md5sum $newname
rc=$?
if [ "$rc" != "0" ]; then
  echo "Command failed: md5sum $newname"
  $topmedcmd mark $bamid cramed failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "md5 calculated in $s seconds"

etime=`date +%s`
etime=`expr $etime - $stime`
here=`pwd`
echo "BAM to CRAM backup completed in $etime seconds, created $here/$newname"
$topmedcmd set $bamid expt_sampleid $nwdid
if [ "$?" != "0" ]; then
  echo "Command failed: $topmedcmd set $bamid expt_sampleid $nwdid"
  $topmedcmd mark $bamid cramed failed
  exit 3
fi

#   All was good
$topmedcmd mark $bamid cped2ncbi delivered