#!/bin/bash
#
#   topmed_cram.sh -submit|-squeezed bamid
#
#	Backup the original BAM as a CRAM file
#
. /usr/cluster/topmed/bin/topmed_actions.inc

bindir=/usr/cluster/bin
ref=/net/mario/gotcloud/gotcloud.ref/hs37d5.fa
illuminaref=/net/topmed/incoming/study.reference/study.reference/illumina.hg19.fa
backupdir=/net/topmed/working/backups

me=cram
markverb=$me
squeezed=n

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me $bamid
  MyRealHost $bamid 'bam'

  #   Get destination directory for backup files
  backupdir=`$topmedpath wherepath $1 backup`
  if [ "$backupdir" = "" ]; then
    Fail "Unable to determine backup directory for '$bamid'. No job submitted."
  fi
  #  Is this squeezed or not?  For now this is only files from the Broad usually,
  #  however, not quite always.  $qual is the number of distinct base call quality 
  #  score values plus one if the bam is not squeezed.  (from Tom)
  sq=''
  qual=`$samtools view $2 2>/dev/null | head -10000 | cut -f 11 | sed 's/./&\n/g' | sort | uniq -c | wc -l`
  if [ $qual -le 11 ]; then
    sq='-squeezed'
  fi
  SubmitJob $bamid "$realhost-$me" '8G' "$0 $sq $*"
  exit
fi

#   Watch for additional options
if [ "$1" = "-squeezed" ]; then
  squeezed=y
  shift
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit|-squeezed] bamid"
  echo ""
  echo "Convert original BAM to a CRAM file and update database"
  exit 1
fi
bamid=$1
bamfile=`$topmedpath wherefile $bamid bam`

#   Get destination directory for backup files
backupdir=`$topmedpath wherepath $bamid backup`
if [ "$backupdir" = "" ]; then
  Fail "Unable to determine backup directory for '$bamid'"
fi

Started
GetNWDID $bamid
stime=`date +%s`
newname="$nwdid.src.cram"

#   Try to force everything as readonly
chmod 555 $bamfile 2> /dev/null

#======================================================================
#   Original input file was a cram
#
#   Some centers deliver us a cram. In this case we do not need to
#   create a cram, but just copy it elsewhere.
#======================================================================
extension="${bamfile##*.}"
if [ "$extension" = "cram" ]; then
  d=`$topmedcmd -persist show $bamid run`
  if [ "$d" = "" ]; then
    Fail "Unable to determine name of run for '$bamid'"
  fi
  offsite=`$topmedcmd -persist show $d offsite`
  if [ "$offsite" = "N" ]; then
    now=`date +%s`
    cp -p $bamfile.crai $backupdir/$newname.crai
    if [ "$?" != "0" ]; then
      Fail "Copy command failed: cp -p $bamfile.crai $backupdir/$newname.crai"
    fi

    cp -p $bamfile $backupdir/$newname
    if [ "$?" != "0" ]; then
      Fail "Copy command failed: cp -p $bamfile $backupdir/$newname"
    fi
    chmod 0444 $backupdir/$newname
    s=`date +%s`; s=`expr $s - $now`; echo "Copy completed in $s seconds"
  else
    echo "Create symlink to original file since backup is offsite"
    mkdir -p $backupdir             # Safe to make backup dir cause it is small
    cd $backupdir
    if [ "$?" != "0" ]; then
        Fail "Unable to CD $backupdir. This directory must be created first."
    fi
    ln -sf $bamfile $newname        # Create backup file as symlink
    ln -sf $bamfile.crai $newname.crai
    if [ ! -f $newname ]; then
        Fail "Unable to create backup file '$newname' in '$backupdir'"
    fi
  fi

  # The md5 for the backup is the same as that for the input
  md5=`$topmedcmd show $bamid checksum`
  SetDB $bamid 'cramchecksum' $md5

  #   Paired reads count for this file is unchanged
  reads=`$topmedcmd -persist show $bamid bamflagstat`
  SetDB $bamid 'cramflagstat' $reads

  #   All was good
  Successful
  etime=`date +%s`
  etime=`expr $etime - $stime`
  Log $etime
  echo "CRAM backup completed, created $backupdir/$newname"
  exit
fi

#======================================================================
#   Original input file was a bam
#======================================================================
#mkdir -p $backupdir            # We run out of space if we make it ourselves
cd $backupdir
if [ "$?" != "0" ]; then
  Fail "Unable to CD $backupdir. This directory must be created first." 
fi

chkname=`basename $bamfile .bam`

#   Get flagstat for original input file
now=`date +%s`
$topmedflagstat $bamfile $bamid bamflagstat /tmp/${chkname}.init.stat
if [ "$?" != "0" ]; then
  Fail "Command failed: $flagstat $bamfile"
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
echo "Creating cram from $bamfile"
if [ "$squeezed" = "n" ]; then
  $bindir/bam squeeze  --in $bamfile  --out - --rmTags "BI:Z;BD:Z;PG:Z"  --keepDups --binMid  --binQualS  2,3,10,20,25,30,35,40,50 | $samtools view  -C -T  $ref  -  >  $newname 
else 
  $samtools view  -C -T $ref  $bamfile  >  $newname
fi
if [ "$?" != "0" ]; then
  Fail "Command failed: $bindir/bam squeeze  --in $bamfile ..." 
fi
s=`date +%s`; s=`expr $s - $now`; echo "Cram created in $s seconds"

#   Build index for the new file
now=`date +%s`
$samtools index $newname
if [ "$?" != "0" ]; then
  Fail "Command failed: $samtools index $newname"
fi
s=`date +%s`; s=`expr $s - $now`; echo "Cram index created in $s seconds"

#   Get flagstat for this cram
now=`date +%s`
$topmedflagstat $newname $bamid cramflagstat /tmp/${chkname}.cram.stat
if [ "$?" != "0" ]; then
  Fail "Command failed: $flagstat $bamfile"
fi
s=`date +%s`; s=`expr $s - $now`; echo "$flagstat for cram completed in $s seconds"

#   Did flagstat output for cram match that for the input bam?
diff=`diff  /tmp/${chkname}.init.stat  /tmp/${chkname}.cram.stat | wc -l`
if [ "$diff" != "0" ]; then
  diff  /tmp/${chkname}.init.stat  /tmp/${chkname}.cram.stat
  Fail "Stat for backup CRAM file differs from that for original BAM"
fi
echo "Stat for CRAM file matches that of original"
rm -f /tmp/${chkname}.init.stat  /tmp/${chkname}.cram.stat

echo "Calculate MD5 for cram"
now=`date +%s`
#   Calculate the MD5 for the cram
md5=`$calcmd5 $newname | awk '{print $1}'`
if [ "$md5" = "" ]; then
  Fail "Command failed: md5sum $newname"
fi
$topmedcmd set $bamid cramchecksum $md5
s=`date +%s`; s=`expr $s - $now`; echo "MD5 calculated in $s seconds"

here=`pwd`
etime=`date +%s`
etime=`expr $etime - $stime`
echo "BAM to CRAM backup completed in $etime seconds, created $here/$newname"
Successful
Log $etime
