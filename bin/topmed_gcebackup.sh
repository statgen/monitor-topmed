#!/bin/bash
#
#   topmed_gcebackup.sh -submit| bamid
#
#	Backup of a sample original file to local storage
#   Note, we no longer copy this file to GCE (July 2018)
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
me=gcebackup
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  RandomRealHost $bamid
  MyRealHost $bamid cram
  MayIRun $me $bamid $realhost
  timeout='2:00:00'
  SubmitJob $bamid $PROJECT '4G' "$0 $*"
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

file=`GetDB $bamid bamname`
if [ "$file" = "" ]; then
  Fail "Unable to determine original source file BAMNAME for '$bamid'"
fi
extension="${file##*.}"

#======================================================================
#   Backup original source file to local storage
#======================================================================
backupfile=`$topmedpath wherefile $bamid localbackup`
echo "Backup of $file to $backupfile"
cp -p $file $backupfile
if [ "$?" != "0" ]; then
  Fail "Failed to backup $file to $backupfile"
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Backup of original source file to local storage completed in $etime seconds"

Successful
Log $etime
exit
