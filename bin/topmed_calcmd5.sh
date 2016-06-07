#!/bin/bash
#
#   topmed_calcmd5.sh - calculate the md5 for a bam or cram
#
#   You might find a truncated cram for testing in
#       /net/topmed/incoming/study.reference/study.reference/NWD850567.trunc.test.cram
#
#   If the file has not completed (no end of file marker, yet)
#   we wait a bit and try again, failing after a few minutes of this.
#
#	prints md5sum checksum to stdout
if [ "$1" = "" ]; then
  echo "$0  cram"
  echo "Use this in other scripts to calculate the md5 for a cram"
  echo "It will attempt to wait until the cram is complete before trying"
  exit 1
fi

samtools=/usr/cluster/bin/samtools
awkcalcend=/usr/cluster/monitor/bin/topmed_calc_cram_end.awk

for n in 1 2 3 4 5 6 7 8 9; do
  offset=`samtools view -H $1 | awk -f $awkcalcend`
  l=`$samtools view $1 $offset 2>&1 | grep truncate`
  if [ "$l" = "" ]; then
    md5sum $1               # Let caller trap results
    exit $?
  fi
  echo "Try $n : File '$1' is not complete, waiting a little bit to try again" >&2
  sleep 60
done
echo "File '$1' was always truncated, giving up" >&2
exit 3
