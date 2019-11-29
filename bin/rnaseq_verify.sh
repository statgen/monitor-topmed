#!/bin/bash
#
#   rnaseq_verify.sh -submit sampleid
#
#	Verify the MD5 checksum for a file file
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
topmedcmd="$topmedcmd -datatype rnaseq"
topmedpath="$topmedpath -datatype rnaseq"
me=verify
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 txseqid`      #	Must set bamid, legacy decision from genome
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='4:00:00'
  SubmitJob $bamid "$PROJECT-$me" '4G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] sampleid"
  echo ""
  echo "Verify checksum for a file and update database"
  exit 1
fi
#	Must set bamid, legacy decision from genome
bamid=$1
sampleid=$1

Started
stime=`date +%s`

d=`$topmedpath wherepath $sampleid releasefiles`	# Where files to check are
cd $d || exit 4

#	Get list of files and verify checksum for each. Do them one at a time for clarity
tmpfile=/run/shm/$$.md5
fail=n
$topmedcmd list files $sampleid | while read l; do
  ll=($l)
  if [ "${ll[1]}" != "" ]; then
    echo "${ll[1]}  ${ll[0]}" > $tmpfile
    md5sum -c $tmpfile
    if [ "$?" != "0" ]; then
      echo "MD5sum failed: `cat $tmpfile`"
      fail=y
    fi
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
