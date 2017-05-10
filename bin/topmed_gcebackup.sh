#!/bin/bash
#
#   topmed_gcebackup.sh -submit| bamid
#
#	Backup CRAM of a sample to GCE
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcebackup
markverb=$me
cores="--cpus-per-task=2"           # Cores here should be same as gsutil

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me  $bamid
  RandomRealHost $bamid
  SubmitJob $bamid "topmed-$me" '2G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Backup CRAM of a sample to GCE"
  exit 1
fi
bamid=$1

file=`$topmedpath wherefile $bamid bam`
if [ "$file" = "" ]; then
  Fail "Unable to determine original BAM or CRAM file for '$bamid'"
fi

Started
GetNWDID $bamid
stime=`date +%s`

datayear=`$topmedcmd show $bamid datayear`
build=`$topmedcmd show $bamid build`
extension="${file##*.}"
#======================================================================
#   Backup CRAM to GCE for only datayear=3 and build=38,
#======================================================================
if [ "$extension" = "cram" -a "$datayear" = "3" -a "$build" = "38" ]; then  
  backupuri=`$topmedpath wherepath $nwdid backup`

  echo "Backup of CRAM to $backupuri/$nwdid.src.cram to GCE"
  f=`basename $cramfile`
  $gsutil cp $cramfile $backupuri/$f
  if [ "$?" != "0" ]; then
    Fail "Failed to copy file to GCE: cp $cramfile $backupuri/$f"
  fi
  echo "Backup of CRAM to Google CLoud completed in $etime seconds"
  $topmedcmd set $bamid 'offsite' D         # Mark this as backed up offsite

  etime=`date +%s`
  etime=`expr $etime - $stime`
  
  Successful
  Log $etime
  exit
fi

#======================================================================
#   Original file was a bam
#======================================================================
if [ "$extension" = "bam" ]; then
  cramfile=`$topmedpath wherefile $bamid cram`
  if [ "$cramfile" = "" ]; then
    Fail "Unable to determine CRAM file for '$bamid'"
  fi
  if [ ! -f "$cramfile" ]; then
    Fail "CRAM file '$cramfile' for '$bamid' was never created"
  fi
  echo "Original BAM was converted by topmed_cram to a cram - nothing to do"
  Successful
  Log 1
  exit
fi

#======================================================================
#   Original file was a cram - should have been backed up offsite
#======================================================================
if [ "$extension" = "cram" ]; then
  offsite=`$topmedcmd show $bamid offsite`
  if [ "$offsite" != "D" ]; then
    Fail "Original file was a CRAM ($bamid) but was not backed up offsite ($offsite)"
  fi
  echo "Original CRAM was backed up offsite"

  cramdir=`$topmedpath wherepath $bamid cram`
  if [ "$cramdir" = "" ]; then
    Fail "Unable to determine backup directory for '$bamid'"
  fi
  echo "Create symlink to original file since backup is offsite"
  mkdir -p $cramdir             # Safe to make backup dir cause it is small
  cd $cramdir
  if [ "$?" != "0" ]; then
      Fail "Unable to CD $cramdir. This directory must be created first."
  fi
  newname=`basename $cramfile`
  ln -sf $bamfile $newname        # Create backup file as symlink in case of screwy names
  ln -sf $bamfile.crai $newname.crai
  if [ ! -f $newname ]; then
    Fail "Unable to create backup file '$newname' in '$cramdir'"
  fi

  Successful
  Log 1
  exit
fi

#   We should never get here
Fail "Somehow $0 $bamid failed to figure out what to backup"
exit 3
