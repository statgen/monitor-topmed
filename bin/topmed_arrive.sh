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
bamid=$1
bamfile=`$topmedpath wherefile $bamid bam`

#   Mark this as started
$topmedcmd -persist mark $bamid $markverb started || exit $?

#   Determine build used to generate this bam/cram
build=''
l=(`$samtools view -H $bamfile | grep '@SQ' | cut -f 1-3 2>/dev/null | head -2`)
if [ "${l[2]}" = "LN:249250621" -a "${l[5]}" = "LN:243199373" ]; then
  build=37
else
  if [ "${l[2]}" = "LN:248956422" -a "${l[5]}" = "LN:242193529" ]; then
    build=38
  else
    a=`echo $bamfile | grep illumina`
    if [ "$a" != "" ]; then       # Path has illumina in it
      build=37
    fi
  fi
fi
if [ "$build" = "" ]; then
  $topmedcmd -emsg "Unable to determine build for '$bamfile' (L2=${l[2]} L5=${l[5]}" mark $bamid $markverb failed
  exit 4
fi
echo "File '$bamfile' [$bamid] is from build $build"
$topmedcmd -persist set $bamid build $build

#   Set NWDID and other database fields
$topmednwd -bamid $bamid $bamfile
rc=$?
if [ "$rc" != "0" ]; then
  $topmedcmd -persist -emsg "$topmednwd failed" mark $bamid $markverb failed
  exit $rc
fi

#   Rename the bamfile to NWD and fix the database
$topmedrename $bamid $bamfile
rc=$?
if [ "$rc" != "0" ]; then
  $topmedcmd -persist -emsg "$topmedrename failed" mark $bamid $markverb failed
  exit $rc
fi
$topmedcmd -persist mark $bamid $markverb completed
exit 0
