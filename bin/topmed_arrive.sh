#!/bin/bash
#
#   topmed_arrive.sh -submit bamid
#
#	Note a BAM file has arrived. Extract header information
#   This command almost certainly runs locally, however
#   the batch -submit has been left in
#
. /usr/cluster/topmed/bin/topmed_actions.inc
topmednwd=/usr/cluster/topmed/bin/topmed_nwdid.pl

me=arrive
markverb=$me

#   We do not expect -submit to be used here
if [ "$1" = "-submit" ]; then
  shift
  echo "Ignoring option -submit"
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me   bamid"
  echo ""
  echo "Mark a bam as arrived and extract details from the header"
  exit 1
fi
bamid=$1                            # This MUST be bamid, may not use NWD yet
#nwdid=`GetNWDID $bamid`            # Only case where we do not get NWDID right away

bamfile=`$topmedpath wherefile $bamid bam`
Started

#   Determine build used to generate this bam/cram
c=`$topmedcmd show $bamid center`
if [ "$c" = "illumina" ]; then      # Hardcode build for this center
  build=37
else
  build=`GetBuild $bamid $bamfile`
fi
echo "File '$bamfile' [$bamid] is from build $build"
SetDB $bamid 'build' $build

#   Set NWDID and other database fields
$topmednwd -bamid $bamid $bamfile
rc=$?
if [ "$rc" != "0" ]; then
  Fail "$topmednwd failed"
fi

#   Rename the bamfile to NWD and fix the database
$topmedrename $bamid $bamfile
rc=$?
if [ "$rc" != "0" ]; then
  Fail "$topmedrename failed"
fi

Successful
exit 0
