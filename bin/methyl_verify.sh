#!/bin/bash
#
#   methyl_verify.sh -submit batchid
#
#	Verify the MD5 checksum for a set of batch files
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
topmedcmd="$topmedcmd -datatype methyl"
topmedpath="$topmedpath -datatype methyl"
me=verify
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=$1
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='4:00:00'
  SubmitJob $bamid "$PROJECT-$me" '4G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] batchid"
  echo ""
  echo "Verify checksum for a file and update database"
  exit 1
fi
#	Must set bamid, legacy decision from genome
bamid=$1
batchid=$1

Started
stime=`date +%s`

d=`$topmedpath wherepath $batchid rundir`	# Where batch data lives
cd $d || exit 4

#	Get list of files and verify checksum for each. Do them one at a time for clarity
tmpfile=/run/shm/$$.md5
fail=n
n=0
$topmedcmd list batch $batchid | while read fileN; do
  if [ "$fileN" != "" ]; then
    md5=`$topmedcmd show $batchid file${n}checksum`
    echo "$md5 $fileN" > $tmpfile
    md5sum -c $tmpfile
    if [ "$?" != "0" ]; then
      echo "MD5sum failed: `cat $tmpfile`"
      fail=y
    fi
    n=`expr $n + 1`
  fi
done
rm -f $tmpfile
if [ "$fail" = "y" ]; then
  Fail "MD5sum calculation failed"
fi

etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM  completed in $etime seconds"

Successful
Log $etime
exit
