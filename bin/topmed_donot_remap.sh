#!/bin/bash
#
#   topmed_donot_remap.sh [runname | runid] build
#
#	Some centers provide CRAMs based a run that we like how it was done
#   Use this script to mark all samples for that run so that it is not remapped
#
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
col=donot_remap

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me  [runname | runid] build"
  echo ""
  echo "Mark all samples for that run so that it is not remapped"
  echo "E.g.  $me 398 38  - marks all samples in runid 398 so build 38 will not be run on them"
  echo "E.g.  $me 2016.1122.whi.01 38  - same as above"
  exit 1
fi
runid=$1
build=$2

#   For now, build can only be 38 until it changes
if [ "$build" != "38" ]; then
  echo "Invalid build '$build'. For now, build can only be 38 until it changes"
  exit 3
fi

#   Get list of samples for this run
tmp=/tmp/donot_remap
$topmedcmd -with-bamid list samples $runid > $tmp || exit $?
if [ "$?" != "0" ]; then
  echo "Unable to get list of samples for run '$runid'"
  exit $rc
fi

#   Mark each sample
while read line; do
  a=($line)
  bamid=${a[1]}
  $topmedcmd set $bamid $col $build || exit 3
done < $tmp
n=`wc -l $tmp | awk '{print $1}'`
echo "Set $n samples '$col=$build' run=$runid"
rm -f $tmp

