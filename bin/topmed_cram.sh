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
topmedpath=/usr/cluster/monitor/bin/topmedpath.pl
topmedflagstat=/usr/cluster/monitor/bin/topmed_flagstat.sh
backupdir=/net/topmed/working/backups
medir=`dirname $0`
calcmd5=/usr/cluster/monitor/bin/topmed_calcmd5.sh
mem=8G
console=/net/topmed/working/topmed-output
markverb=cramed
squeezed=n
constraint="--constraint eth-10g"
qos=''
slurmp=topmed
realhost=''

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test cram $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  # Run this on node where bam lives
  h=`$topmedpath whathost $1 bam`
  if [ "$h" != "" ]; then
    realhost="--nodelist=$h"
    qos="--qos=$h-cram"
  fi

  #   Get destination directory for backup files
  backupdir=`$topmedpath wherepath $1 backup`
  if [ "$backupdir" = "" ]; then
    echo "Unable to determine backup directory for '$bamid'. No job submitted."
    exit 2
  fi

  #  Is this squeezed or not?  For now this is only files from the Broad usually,
  #  however, not quite always.  $qual is the number of distinct base call quality 
  #  score values plus one if the bam is not squeezed.  (from Tom)
  sq=''
  qual=`$samtools view $2 2>/dev/null | head -10000 | cut -f 11 | sed 's/./&\n/g' | sort | uniq -c | wc -l`
  if [ $qual -le 11 ]; then
    sq='-squeezed'
  fi

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-cram --output=$console/$1-cram.out $0 $sq $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-cram --output=$console/$1-cram.out $0 $sq $*"
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
  echo "Convert original BAM to a CRAM file and update database"
  exit 1
fi
bamid=$1
bamfile=$2

#   Is this a cram or bam
extension="${bamfile##*.}"

#   Maybe this is already running or completed?
s=`$topmedcmd show $bamid state_cram`
if [ "$s" = "20" -o "$s" = "3" ]; then
  echo "bamid=$bamid bamfile=$bamfile already running or completed"
  exit
fi

#   Get destination directory for backup files
backupdir=`$topmedpath wherepath $bamid backup`
if [ "$backupdir" = "" ]; then
  echo "Unable to determine backup directory for '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
s=`hostname`
p=`pwd`
echo "#========= '$d' host=$s $SLURM_JOB_ID squeezed=$squeezed $0 bamid=$bamid bamfile=$bamfile backupdir=$backupdir pwd=$p========="

#   Mark this as started
$topmedcmd mark $bamid $markverb started
stime=`date +%s`
#   Try to force everything as readonly
chmod 555 $bamfile 2> /dev/null

#   Be sure that NWDID is set in database
nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "NWDID not set for $bamid $bamfile"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
newname="$nwdid.src.cram"

#======================================================================
#   Original input file was a cram
#
#   Some centers deliver us a cram. In this case we do not need to
#   create a cram, but just copy it elsewhere.
#======================================================================
if [ "$extension" = "cram" ]; then
  d=`$topmedcmd -persist show $bamid run`
  if [ "$d" = "" ]; then
    echo "Unable to determine name of run for '$bamid'"
    $topmedcmd -persist mark $bamid $markverb failed
    exit 3
  fi
  offsite=`$topmedcmd -persist show $d offsite`
  if [ "$offsite" = "N" ]; then
    now=`date +%s`
    cp -p $bamfile.crai $backupdir/$newname.crai
    if [ "$?" != "0" ]; then
      echo "Copy command failed: cp -p $bamfile.crai $backupdir/$newname.crai"
      $topmedcmd -persist mark $bamid $markverb failed
      exit 3
    fi

    cp -p $bamfile $backupdir/$newname
    if [ "$?" != "0" ]; then
      echo "Copy command failed: cp -p $bamfile $backupdir/$newname"
      $topmedcmd -persist mark $bamid $markverb failed
      exit 3
    fi
    chmod 0444 $backupdir/$newname
    s=`date +%s`; s=`expr $s - $now`; echo "Copy completed in $s seconds"
  else
    echo "Create symlink to original file since backup is offsite"
    mkdir -p $backupdir             # Safe to make backup dir cause it is small
    cd $backupdir
    if [ "$?" != "0" ]; then
        echo "Unable to CD $backupdir. This directory must be created first."
        $topmedcmd -persist mark $bamid $markverb failed
        exit 2
    fi
    ln -sf $bamfile $newname        # Create backup file as symlink
    ln -sf $bamfile.crai $newname.crai
    if [ ! -f $newname ]; then
        echo "Unable to create backup file '$newname' in '$backupdir'"
        $topmedcmd -persist mark $bamid $markverb failed
        exit 2
    fi
  fi

  # The md5 for the backup is the same as that for the input
  md5=`$topmedcmd show $bamid checksum`
  $topmedcmd -persist set $bamid cramchecksum $md5

  #   Paired reads count for this file is unchanged
  reads=`$topmedcmd -persist show $bamid bamflagstat`
  $topmedcmd -persist set $bamid cramflagstat $reads

  #   All was good
  $topmedcmd -persist mark $bamid $markverb completed
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo `date` cram $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
  echo "CRAM backup completed, created $backupdir/$newname"
  exit
fi

#======================================================================
#   Original input file was a bam
#======================================================================
#mkdir -p $backupdir            # We run out of space if we make it ourselves
cd $backupdir
if [ "$?" != "0" ]; then
  echo "Unable to CD $backupdir. This directory must be created first."
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

chkname=`basename $bamfile .bam`

#   Get flagstat for original input file
now=`date +%s`
$topmedflagstat $bamfile $bamid bamflagstat /tmp/${chkname}.init.stat
if [ "$?" != "0" ]; then
  echo "Command failed: $flagstat $bamfile"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "$flagstat completed in $s seconds"

#   Illumina cram flles require a different fasta
center=`$topmedcmd show $bamid center`
if [ "$center" = "illumina" ]; then
  ref=$illuminaref
  echo "BAM to CRAM for '$center' requires a different fasta file '$ref'"
fi

now=`date +%s`
#   Create the CRAM file
#   Run 'bam squeeze' on .bam file with result written as .cram
#   Checking with 'samtools flagstat' is quick, zero means success.
#   Running 'samtools index' on a .cram file does NOT require 
#   an explicit genome reference sequence, but 'samtools view' 
#   does.  (And, 'samtools flagstat' does accept .sam input.)
if [ "$squeezed" = "n" ]; then
  $bindir/bam squeeze  --in $bamfile  --out - --rmTags "BI:Z;BD:Z;PG:Z"  --keepDups --binMid  --binQualS  2,3,10,20,25,30,35,40,50 | $samtools view  -C -T  $ref  -  >  $newname 
else 
  $samtools view  -C -T $ref  $bamfile  >  $newname
fi
if [ "$?" != "0" ]; then
  echo "Command failed: $bindir/bam squeeze  --in $bamfile ..."
  $topmedcmd -persist mark $bamid $markverb failed
  exit 1
fi
s=`date +%s`; s=`expr $s - $now`; echo "Cram created in $s seconds"

#   Build index for the new file
now=`date +%s`
$samtools index $newname
if [ "$?" != "0" ]; then
  echo "Command failed: $samtools index $newname"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "Cram index created in $s seconds"

#   Get flagstat for this cram
now=`date +%s`
$topmedflagstat $newname $bamid cramflagstat /tmp/${chkname}.cram.stat
if [ "$?" != "0" ]; then
  echo "Command failed: $flagstat $bamfile"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
s=`date +%s`; s=`expr $s - $now`; echo "$flagstat for cram completed in $s seconds"

#   Did flagstat output for cram match that for the input bam?
diff=`diff  /tmp/${chkname}.init.stat  /tmp/${chkname}.cram.stat | wc -l`
if [ "$diff" != "0" ]; then
  echo "Stat for backup CRAM file differs from that for original BAM"
  diff  /tmp/${chkname}.init.stat  /tmp/${chkname}.cram.stat
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
echo "Stat for CRAM file matches that of original"
rm -f /tmp/${chkname}.init.stat  /tmp/${chkname}.cram.stat

echo "Calculate MD5 for cram"
now=`date +%s`
#   Calculate the MD5 for the cram
md5=`$calcmd5 $newname | awk '{print $1}'`
if [ "$md5" = "" ]; then
  echo "Command failed: md5sum $newname"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
$topmedcmd set $bamid cramchecksum $md5
s=`date +%s`; s=`expr $s - $now`; echo "MD5 calculated in $s seconds"

here=`pwd`
etime=`date +%s`
etime=`expr $etime - $stime`
echo "BAM to CRAM backup completed in $etime seconds, created $here/$newname"
$topmedcmd -persist mark $bamid $markverb completed
echo `date` cram $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
