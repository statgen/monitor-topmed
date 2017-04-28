#!/bin/bash
#
#   topmedforcew.sh dir
#
#	Force in incoming directory to be writable so the monitor can deal with it

#   We only play with files in /incoming/topmed
topdir=/incoming/topmed
topuser=topmed
a=`echo $1 | grep $topdir`
if [ "$a" = "" ]; then
  echo "$0 $topdir/center/run"
  echo ""
  echo "Force directories in $topdir to be writable for user $topuser"
  exit 1
fi

cd $1
if [ "$?" != "0" ]; then
  echo "$0 - Unable to CD to '$1'"
  exit 4
fi
# echo "Forcing '$1' to be writable for the monitor"
chgrp topmed .
chmod 770 .
chmod 660 Manifest.txt *.md5 *.cram *.crai *.bam *.bai 2> /dev/null   # Failure is OK
