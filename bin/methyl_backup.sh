#!/bin/bash
#
#   methyl_backup.sh -submit| batchid
#
#	Backup of a batch of files to local storage
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
topmedcmd="$topmedcmd -datatype methyl"
topmedpath="$topmedpath -datatype methyl"
me=backup
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=$1
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='2:00:00'
  SubmitJob $bamid "$PROJECT-$me" '4G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] batchid"
  echo ""
  echo "Backup data of a sample to GCE"
  exit 1
fi

#	Must set bamid, legacy decision from genome
bamid=$1
batchid=$1

Started
stime=`date +%s`

#======================================================================
#   Just once we want to backup the toplevel directory
#======================================================================
backupdir=`$topmedpath wherefile $batchid localbackup`
origrundir=`$topmedpath wherepath $batchid rundir`
batchname=`$topmedcmd show $batchid batchname`

cd $origrundir
if [ "$?" != "0" ]; then
  Fail "Unable to find rundir '$origrundir' for samples '$batchid'"
fi
if [ ! -f "$backupdir/Manifest.txt" ]; then		# Backup toplevel directory?
  mkdir -p $backupdir || chmod 0770 $backupdir	# May already exist
  cp -p `ls|grep -v zip` $backupdir				# Copy non-zip files
  echo "Backup of $batchname toplevel files created in $backupdir"
else
  echo "Not necessary to backup $batchname toplevel files"
fi

#======================================================================
#   Backup original source files for sample to local storage
#======================================================================
n=0
for f in `$topmedcmd list batch $batchid`; do
  echo "Backup of $batchname data to $backupdir"
  cp --preserve=timestamps $f $backupdir
  if [ "$?" != "0" ]; then
    Fail "Failed to backup batch '$batchname' [$f] to $backupdir"
  fi
  n=`expr $n + 1`
done
if [ "$n" != "4" ]; then
  Fail "Did not copy enough zip $batchname batch files for $batchid"
fi

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Backup $n files for $batchname to $backupdir completed in $etime seconds"

Successful
Log $etime
exit
