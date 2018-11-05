#!/bin/bash
#
#   topmed_backup.sh -submit| bamid
#
#	Backup of a sample original file to local storage
#   Note, we no longer copy this file to GCE (July 2018)
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
me=backup
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='2:00:00'
  SubmitJob $bamid "$PROJECT-$me" '4G' "$0 $*"
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
nwdid=`GetNWDID $bamid`
bamid=`GetDB $nwdid bamid`

Started
stime=`date +%s`

#======================================================================
#   Backup original source file to local storage
#======================================================================
origfile=`$topmedpath wherefile $bamid bam`
backupfile=`$topmedpath wherefile $bamid localbackup`
if [ ! -f $origfile ]; then     # If orig does not exist, try cram
   origfile=`$topmedpath wherefile $bamid cram`
   f=`basename $origfile`
   backupfile=`$topmedpath wherefile $bamid localbackup`/$f 
fi
echo "Backup of $origfile to $backupfile"
d=`dirname $backupfile`
mkdir -p $d                         # Just in case center dir needs making
chmod 0770 $d
cp --preserve=timestamps $origfile $backupfile
if [ "$?" != "0" ]; then
  Fail "Failed to backup $origfile to $backupfile"
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Backup of original source file to local storage completed in $etime seconds"

Successful
Log $etime
exit
