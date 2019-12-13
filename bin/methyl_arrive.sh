#!/bin/bash
#
#   methyl_arrive.sh -submit batchid
#
#	Note that methylation batch files for a sample have arrived
#   This command runs locally
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
topmedcmd="$topmedcmd -datatype methyl"
topmedpath="$topmedpath -datatype methyl"
me=arrive
markverb=$me

#   We do not expect -submit to be used here
if [ "$1" = "-submit" ]; then
  shift
  echo "Ignoring option -submit"
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me   sampleid"
  echo ""
  echo "Mark a sample as arrived and extract any details as needed"
  exit 1
fi
#	Must set bamid, legacy decision from genome
bamid=$1                       	# This MUST be sampleid, not TOR 
batchid=$bamid

Started

Successful
exit 0
