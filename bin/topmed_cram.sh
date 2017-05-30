#!/bin/bash
#
#   topmed_cram.sh -submit|-squeezed bamid
#
#	Convert the original BAM to a CRAM
#
. /usr/cluster/topmed/bin/topmed_actions.inc

bindir=/usr/cluster/bin
ref=/net/mario/gotcloud/gotcloud.ref/hs37d5.fa
illuminaref=/net/topmed/incoming/study.reference/study.reference/illumina.hg19.fa

me=cram
markverb=$me
squeezed=n

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  MayIRun $me $bamid
  MyRealHost $bamid 'bam'

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

Started
GetNWDID $bamid
stime=`date +%s`

extension="${bamfile##*.}"
#======================================================================
#   If original file was a CRAM, nothing to do, just bookkeeping
#======================================================================
if [ "$extension" = "cram" ]; then  
  # Check if cram directory exists for this bamid
  # If not, attempt to create ot (e.g. link cram direct to original input run)
  cramdir=`$topmedpath wherepath $bamid cram`
  if [ ! -d "$cramdir" ]; then
    echo "No CRAM path known, we will attempt to create a symlink to the original run"
    center=`GetDB $bamid center`
    run=`GetDB $bamid run`
    #   Go to where backups directory for this run should be
    d=/net/topmed/working/backups/incoming/topmed/$center
    cd $d
    if [ "$?" != "0" ];  then
      Fail "Unable to CD to CRAM directory $d" 
    fi
    #   Get home of original source files  (e.g. /net/topmed/incoming/topmed/...)
    #   Now make the /net be relative so it could work in FLUX someday
    bdir=`$topmedpath wherepath $bamid bam | sed -e 's:/net:../../../../../..:'`
    ln -s $bdir .
    #   Check if this worked
    cramdir=`$topmedpath wherepath $bamid cram`
    if [ "$cramdir" = "" ];  then
      Fail "Unable to create CRAM directory for '$bamid'" 
    else
      if [ ! -d $cramdir ]; then
        Fail "Cannot find directory to backup cram (e.g. original source dir)"
      fi
    fi
    echo "Created the CRAM directory for $bamid: $cramdir"
  fi
  #  Be sure the cram file exists or make a symlink to it
  f="$nwdid.src.cram"
  if [ ! -f "$cramdir/$f" ]; then
    cd $cramdir
    if [ "$?" != "0" ]; then
      Fail "Unable to CD to cramdir '$cramdir'"
    fi
    bfile=`GetDB $bamid bamname`
    if [ "$?" != "0" ]; then
      Fail "Unable to find bamname for '$bamid'"
    fi
    ln -s $bfile $f
    if [ "$?" != "0" ]; then
      Fail "Unable to create symlink to '$f' for '$bamid'"
    fi
    echo "Created symlink for $f"
    if [ ! -f "$f.crai" ]; then         # Create link to crai as needed too
      ln -s $bfile.crai $f.crai
      if [ "$?" != "0" ]; then
        Fail "Unable to create symlink to '$f.crai' for '$bamid'"
      fi
      echo "Created symlink for $f.crai"
    fi
  fi

  echo "MD5 for cram, same as original file"
  v=`GetDB $bamid checksum`
  SetDB $bamid 'cramchecksum' $v

  v=`GetDB $bamid checksum`
  SetDB $bamid 'cramchecksum' $v

  flagstat=`GetDB $bamid cramflagstat`
  if [ "$flagstat" = "" ]; then
    echo "New CRAM, calculate flagstat for $cramfile"
    now=`date +%s`
    $topmedflagstat $cramfile $bamid cramflagstat
    if [ "$?" != "0" ]; then
      Fail "Command failed: $flagstat $bamfile"
    fi
    s=`date +%s`; s=`expr $s - $now`;
    echo "$flagstat for CRAM completed in $s seconds"
  else 
    echo "Flagstat for cram, same as original file"
    v=`GetDB $bamid bamflagstat`
    SetDB $bamid 'cramflagstat' $v
    echo "No CRAM actually created, just did database bookkeeping"
  fi
  Successful
  Log 1
  exit
fi

#======================================================================
#   Original input file was a BAM, convert to a CRAM
#======================================================================
cramdir=`$topmedpath wherepath $bamid cram`
if [ "$cramdir" = "" ]; then
  Fail "Unable to CD $cramdir. This directory must be created first." 
fi
cd $cramdir
chkname=`basename $bamfile .bam`

#   Get flagstat for original input file
now=`date +%s`
$topmedflagstat $bamfile $bamid bamflagstat /tmp/${chkname}.init.stat
if [ "$?" != "0" ]; then
  Fail "Command failed: $flagstat $bamfile"
fi
s=`date +%s`; s=`expr $s - $now`; echo "$flagstat completed in $s seconds"

#   Illumina cram files require a different fasta
center=`GetDB $bamid center`
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
newname=`$topmedpath wherefile $bamid cram`
newname=`basename $newname`
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
s=`date +%s`; s=`expr $s - $now`; echo "CRAM index created in $s seconds"

#   Get flagstat for this cram
now=`date +%s`
$topmedflagstat $newname $bamid cramflagstat /tmp/${chkname}.cram.stat
if [ "$?" != "0" ]; then
  Fail "Command failed: $flagstat $bamfile"
fi
s=`date +%s`; s=`expr $s - $now`; echo "$flagstat for CRAM completed in $s seconds"

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
SetDB $bamid 'cramchecksum' $md5
s=`date +%s`; s=`expr $s - $now`; echo "MD5 calculated in $s seconds"

here=`pwd`
etime=`date +%s`
etime=`expr $etime - $stime`
echo "BAM to CRAM backup completed in $etime seconds, created $here/$newname"
Successful
Log $etime
