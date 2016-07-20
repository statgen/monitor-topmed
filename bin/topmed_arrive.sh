#!/bin/bash
#
#   topmed_arrive.sh -submit bamid bamfile
#
#	Note a BAM file has arrived. Extract header information
#   This command almost certainly runs locally, however
#   the batch -submit has been left in
#
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
topmednwd=/usr/cluster/monitor/bin/topmed_nwdid.pl
topmedrename=/usr/cluster/monitor/bin/topmedrename.pl
console=/net/topmed/working/topmed-output

#   We do not expect -submit to be used here
if [ "$1" = "-submit" ]; then
  shift
  echo "Ignoring option -submit"
fi

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me   bamid bamfile"
  echo ""
  echo "Mark a bam as arrived and extract details from the header"
  exit 1
fi
bamid=$1
bamfile=$2

#   Mark this as started
$topmedcmd mark $bamid arrived started || exit $?

#   Set NWDID and other database fields
$topmednwd -bamid $bamid $bamfile
rc=$?
if [ "$rc" != "0" ]; then
  $topmedcmd mark $bamid arrived failed
  exit $rc
fi

#   Rename the bamfile to NWD and fix the database
$topmedrename $bamid $bamfile
rc=$?
if [ "$rc" != "0" ]; then
  $topmedcmd mark $bamid arrived failed
  exit $rc
fi
$topmedcmd mark $bamid arrived completed
exit 0
