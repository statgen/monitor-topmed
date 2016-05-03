#!/bin/bash
#
#   topmed_cram.sh -submit|-squeezed bamid bamfile
#
#	Backup the original BAM as a CRAM file
#
bindir=/usr/cluster/bin
samtools=/net/mario/gotcloud/bin/samtools
ref=/net/mario/gotcloud/gotcloud.ref/hs37d5.fa
illuminaref=/net/topmed/incoming/study.reference/study.reference/illumina.hg19.fa
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
backupdir=/net/topmed/working/backups
medir=`dirname $0`
calcmd5=$medir/topmed_calcmd5.sh
mem=8G
if [ "$TOPMED_MEMORY" != "" ]; then mem=$TOPMED_MEMORY; fi
realhost=topmed
#if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
console=/net/topmed/working/topmed-output
markverb=cramed
squeezed=n
qos=cram
if [ "$TOPMED_QOS" != "" ]; then qos=$TOPMED_QOS; fi

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test cram $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #   Figure where to submit this to run - should be local
  l=(`$topmedcmd where $1 backup`)     # Get backupdir and backupfile and host
  h="${l[2]}"
  if [ "$h" != "" ]; then realhost=$h; fi
  if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
  slurmp="$realhost-incoming"
  slurmqos="$realhost-$qos"

  #  Is this squeezed or not?  For now this is only files from the Broad usually,
  #  however, not quite always.  $qual is the number of distinct base call quality 
  #  score values plus one if the bam is not squeezed.  (from Tom)
  sq=''
  qual=`$samtools view $2 2>/dev/null | head -10000 | cut -f 11 | sed 's/./&\n/g' | sort | uniq -c | wc -l`
  if [ $qual -le 11 ]; then
    sq='-squeezed'
  fi

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-cram --output=$console/$1-cram.out $0 $sq $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-cram --output=$console/$1-cram.out $0 $sq $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` cram ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
  exit
fi

#   Watch for additional options
if [ "$1" = "-squeezed" ]; then
  squeezed=y
  shift
fi

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit|-squeezed] bamid bamfile"
  echo ""
  echo "Backup the original BAM as a CRAM file and update database"
  exit 1
fi
bamid=$1
bamfile=$2

#   Get destination directory for backup files
l=(`$topmedcmd where $1 backup`)        # Get backupdir and backupfile and host
backupdir="${l[0]}"

d=`date +%Y/%m/%d`
mkdir -p $backupdir
cd $backupdir
if [ "$?" != "0" ]; then
  echo "Unable to CD $backupdir"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
s=`hostname`
echo "#========= '$d' host=$s $SLURM_JOB_ID squeezed=$squeezed $0 bamid=$bamid bamfile=$bamfile backupdir=$backupdir ========="

#   Mark this as started
$topmedcmd mark $bamid $markverb started
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
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
newname="$nwdid.src.cram"

now=`date +%s`
$samtools flagstat  $bamfile >  ${chkname}.init.stat
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools flagstat  $bamfile"
  $topmedcmd mark $bamid $markverb failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "samtools flagstat completed in $s seconds"

now=`date +%s`

#   This was added so we could remake the cram files for strange illumina data
center=`$topmedcmd show $bamid center`
if [ "$center" = "illumina" ]; then
  ref=$illuminaref
  echo "BAM to CRAM for '$center' requires a different fasta file '$ref'"
fi

#   Create the CRAM file
if [ "$squeezed" = "n" ]; then
  $bindir/bam squeeze  --in $bamfile  --out - --rmTags "BI:Z;BD:Z;PG:Z"  --keepDups --binMid  --binQualS  2,3,10,20,25,30,35,40,50 | $samtools view  -C -T  $ref  -  >  $newname 
else 
  $samtools view  -C -T $ref  $bamfile  >  $newname
fi
if [ "$?" != "0" ]; then
  echo "Command failed: $bindir/bam squeeze  --in $bamfile ..."
  $topmedcmd mark $bamid $markverb failed
  exit 1
fi
s=`date +%s`; s=`expr $s - $now`; echo "squeeze completed in $s seconds"

now=`date +%s`
$samtools index $newname
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools index $newname"
  $topmedcmd mark $bamid $markverb failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "samtools index completed in $s seconds"

now=`date +%s`
$samtools view  -u -T  $ref $newname | $samtools flagstat  - >  ${chkname}.cram.stat
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools view  -h -T  $ref $newname ..."
  $topmedcmd mark $bamid $markverb failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "samtools view completed in $s seconds"

diff=`diff  ${chkname}.init.stat  ${chkname}.cram.stat | wc -l`
if [ "$diff" != "0" ]; then
  echo "Stat for backup CRAM file differs from that for original BAM"
  diff  ${chkname}.init.stat  ${chkname}.cram.stat
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
echo "Stat for CRAM file matches that of original"
rm -f ${chkname}.init.stat  ${chkname}.cram.stat ${chkname}.bam ${chkname}.bam.md5  ${chkname}.bai ${chkname}.bai.md5

echo "Calculating new MD5"
now=`date +%s`
#   Calculate the MD5 for the cram
md5=`$calcmd5.sh $newname | awk '{print $1}'`
if [ "$md5" = "" ]; then
  echo "Command failed: md5sum $newname"
  $topmedcmd mark $bamid $markverb failed
  exit 3
fi
$topmedcmd set $bamid cramchecksum $md5
s=`date +%s`; s=`expr $s - $now`; echo "md5 calculated in $s seconds"

etime=`date +%s`
etime=`expr $etime - $stime`

#   Be sure that NWDID is set in database
here=`pwd`
echo "BAM to CRAM backup completed in $etime seconds, created $here/$newname"
$topmedcmd set $bamid expt_sampleid $nwdid
if [ "$?" != "0" ]; then
  echo "Command failed: $topmedcmd set $bamid expt_sampleid $nwdid"
  $topmedcmd mark $bamid $markverb failed
  exit 3
fi

#   All was good
$topmedcmd mark $bamid $markverb completed
echo `date` cram $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
