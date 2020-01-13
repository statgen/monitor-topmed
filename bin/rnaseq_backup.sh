#!/bin/bash
#
#   rnaseq_backup.sh -submit| sampleid
#
#	Backup of a sample original file to local storage
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
topmedcmd="$topmedcmd -datatype rnaseq"
topmedpath="$topmedpath -datatype rnaseq"
me=backup
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 txseqid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='2:00:00'
  SubmitJob $bamid "$PROJECT-$me" '4G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] sampleid"
  echo ""
  echo "Backup data of a sample to GCE"
  exit 1
fi

#	Must set bamid, legacy decision from genome
bamid=$1
sampleid=$1
nwdid=`GetNWDID $sampleid`
sampleid=`GetDB $nwdid txseqid`
center=`$topmedcmd show $sampleid center`

Started
stime=`date +%s`

#======================================================================
#   Just once we want to backup the toplevel directory
#======================================================================
backupdir=`$topmedpath wherefile $sampleid localbackup`
origrundir=`$topmedpath wherepath $sampleid rundir`
releasefiles=`$topmedpath wherepath $sampleid releasefiles`
releasefiles=`basename $releasefiles`
if [ "releasefiles" = '' ]; then
  Fail "Failed to get releasefiles directory for sample '$sampleid'"
fi
if [ "$center" = "broad" ]; then	# Broad does it differently :-(
  releasefiles="releaseFiles/$releasefiles"
fi

cd $origrundir
if [ "$?" != "0" ]; then
  Fail "Unable to find rundir '$origrundir' for sample '$sampleid'"
fi

#	Set up target directories for this sample. Copy meta data for all samples too
if [ ! -f "$backupdir/Manifest.txt" ]; then
  mkdir -p $backupdir || chmod 0770 $backupdir	# May already exist
  cp -p *.* $backupdir				# Copy top level files, error on directories
  mkdir -p $backupdir/$releasefiles	# Samples copied here
  # Broad uses a different structure   These commands may fail for some projects
  if [ "$center" = "broad" ]; then	# Broad does it differently :-(
    cp -rp releaseFiles/data  $backupdir/releaseFiles
    cp -rp releaseFiles/hound  $backupdir/releaseFiles
  fi
  echo "Backup of toplevel files created in $backupdir"
else
  echo "Not necessary to backup toplevel files"
fi

#======================================================================
#   Backup original source files for sample to local storage
#======================================================================
echo "Backup of '$sampleid' [$nwdid] data to $backupdir"
if [ "$center" != "broad" ]; then	# Other centers do it differently
  fileprefix=`$topmedcmd show $sampleid fileprefix`
  if [ "$?" != "0" ]; then
    Fail "Failed to get fileprefix for sample '$sampleid'"
  fi
  cp --preserve=timestamps $releasefiles/$fileprefix.* $backupdir/$releasefiles
  if [ "$?" != "0" ]; then
    Fail "Failed to backup sample '$sampleid' [$fileprefix] to $backupdir/$releasefiles"
  fi
fi
if [ "$center" = "broad" ]; then	# Broad does it this way
  cp --preserve=timestamps $releasefiles/$nwdid.* $backupdir/$releasefiles
  if [ "$?" != "0" ]; then
    Fail "Failed to backup sample '$sampleid' [$nwdid] to $backupdir/$releasefiles"
  fi
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Backup $releasefiles/$nwdid to $backupdir/$releasefiles completed in $etime seconds"

Successful
Log $etime
exit
