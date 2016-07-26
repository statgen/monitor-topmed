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
samtools=/net/mario/gotcloud/bin/samtools

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

#   Determine build used to generate this bam/cram
build=''
l=(`$samtools view -H $bamfile | grep '@SQ' | cut -f 1-3 | head -2`)
if [ "${l[2]}" = "LN:249250621" -a "${l[5]}" = "LN:243199373" ]; then
  build=37
else
  if [ "${l[2]}" = "LN:248956422" -a "${l[5]}" = "LN:242193529" ]; then
    build=38
  else
    a=`echo $1 | grep illumina`
    if [ "$a" != "" ]; then       # Path has illumina in it
      build=37
    fi
  fi
fi
if [ "$build" = "" ]; then
  echo "Unable to determine build for '$bamfile' (L2=${l[2]} L5=${l[5]}"
  $topmedcmd mark $bamid arrived failed
  exit 4

fi

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
  $topmedcmd -persist mark $bamid arrived failed
  exit $rc
fi
$topmedcmd -persist mark $bamid arrived completed
exit 0
