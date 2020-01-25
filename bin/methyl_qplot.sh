#!/bin/bash
#
#   methyl_qplot.sh -submit| sampleid
#
#	Run quality controls on methylation data
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
topmedcmd="$topmedcmd -datatype methyl"
topmedpath="$topmedpath -datatype methyl"
me=qplot
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='2:00:00'
  SubmitJob $bamid "$PROJECT-$me" '8G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Run quality controls on methylation data"
  exit 1
fi
#	Must set bamid, legacy decision from genome

bamid=$1
sampleid=$1
nwdid=`GetNWDID $sampleid`

#Started
stime=`date +%s`

echo "I don't know how to find the source BAM data for this"
echo "This script is not ready yet"; exit 2

#	Figure out file to run QC on
d=`$topmedpath wherepath $sampleid releasefiles`	# Where files to check are
cd $d || exit 4

bam=`ls -1 | grep $nwdid | grep -v .bai`
if [ "$bam" = "" ]; then
  Fail "$me Unable to find BAM for $sampleid"
fi

#	Run quality control here. Check all return codes and set fail=y or use Fail
fail=y
echo "Run QC on $d/$bam"

if [ "$fail" = "y" ]; then
  Fail "$me Quality control calculation failed"
fi

etime=`date +%s`
etime=`expr $etime - $stime`
echo "$me  completed in $etime seconds"

Successful
Log $etime
exit

