#!/bin/bash
#
#   topmed_cram.sh -submit bamid bamfile
#
#	Backup the original BAM as a CRAM file
#
bindir=/usr/cluster/bin
samtools=/net/mario/gotcloud/bin/samtools
ref=/net/mario/gotcloud/gotcloud.ref/hs37d5.fa
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
backupdir=/net/topmed/working/backups
mem=${TOPMED_MEM:-8G}
console=/net/topmed/working/topmed-output
slurmp=${TOPMED_PARTITION:-nomosix}
slurmqos=topmed-cram

if [ "$1" = "-submit" ]; then
  shift  
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-cram --output=$console/$1-cram.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-cram --output=$console/$1-cram.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 cramed submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobid
    $topmedcmd set $1 jobidcram ${l[3]}
  fi
  exit
fi

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid bamfile"
  echo ""
  echo "Backup a bam file and update database"
  exit 1
fi
bamid=$1
bamfile=$2

#   This is a bit of a hack to insure backups of data on topmed
#   go to topmed2 and vice versa
homehost=`echo $2 | cut -d / -f 3`      # Should be topmed/topmed2
if [ "$homehost" != "" ]; then
  if [ "$homehost" = "topmed" ]; then
    homehost=topmed2
  fi
  if [ "$homehost" = "topmed2" ]; then
    homehost=topmed
  fi
  backupdir=/net/$homehost/working/backups
fi

#   Mark this as started  Calc dest directory in a poor manner
d=`echo $bamfile | sed -e s:/net/topmed/::`
d=`echo $d | sed -e s:/net/topmed2/::`
d=`echo $d | sed -e s:/net/topmed3/::`
dest=`dirname $d`
d=`date +%Y/%m/%d`
cd $backupdir
if [ "$?" != "0" ]; then
  echo "Unable to CD $backupdir"
  $topmedcmd mark $bamid cramed failed
  exit 2
fi
mkdir -p $dest
cd $dest
if [ "$?" != "0" ]; then
  echo "Unable to CD $backupdir/$dest"
  $topmedcmd mark $bamid cramed failed
  exit 2
fi
s=`hostname`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid bamfile=$bamfile dest=$backupdir/$dest ========="

#   Mark this as started
$topmedcmd mark $bamid cramed started
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
#       cp -p $bamfile $backupdir/$dest
#
chkname=`basename $bamfile .bam`
nwdid=`$samtools view -H $bamfile | grep '^@RG' | grep -o 'SM:\S*' | sort -u | cut -d \: -f 2`
if [ "$nwdid" = "" ]; then
  echo "Unable to extract NWDID from header of '$bamfile'"
  $topmedcmd mark $bamid cramed failed
  exit 2
fi
newname="$nwdid.src.cram"

now=`date +%s`
$samtools flagstat  $bamfile  >  ${chkname}.init.stat
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools flagstat  $bamfile"
  $topmedcmd mark $bamid cramed failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "samtools flagstat completed in $s seconds"

now=`date +%s`
$bindir/bam squeeze  --in $bamfile  --out - --rmTags "BI:Z;BD:Z;PG:Z"  --keepDups --binMid  --binQualS  2,3,10,20,25,30,35,40,50 | $samtools view  -C -T  $ref  -  >  $newname
if [ "$?" != "0" ]; then
  echo "Command failed: $bindir/bam squeeze  --in $bamfile ..."
  $topmedcmd mark $bamid cramed failed
  exit 1
fi
s=`date +%s`; s=`expr $s - $now`; echo "squeeze completed in $s seconds"

now=`date +%s`
$samtools index $newname
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools index $newname"
  $topmedcmd mark $bamid cramed failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "samtools index completed in $s seconds"

now=`date +%s`
$samtools view  -h -T  $ref $newname | $samtools flagstat  - >  ${chkname}.cram.stat
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools view  -h -T  $ref $newname ..."
  $topmedcmd mark $bamid cramed failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "samtools view completed in $s seconds"

diff=`diff  ${chkname}.init.stat  ${chkname}.cram.stat | wc -l`
if [ "$diff" != "0" ]; then
  echo "Stat for backup CRAM file differs from that for original BAM"
  diff  ${chkname}.init.stat  ${chkname}.cram.stat
  $topmedcmd mark $bamid cramed failed
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
$topmedcmd set $bamid nwdid $nwdid
if [ "$?" != "0" ]; then
  echo "Command failed: $topmedcmd set $bamid nwdid $nwdid"
  $topmedcmd mark $bamid cramed failed
  exit 3
fi

#   All was good
$topmedcmd mark $bamid cramed completed
