#!/bin/bash
#
#   topmed_calcmd5.sh - calculate the md5 for a bam or cram
#
#   If the file has not completed (no end of file marker, yet)
#   we wait a bit and try again, failing after a few minutes of this.
#
#	prints md5sum checksum to stdout
samtools=/usr/cluster/bin/samtools

for n in 1 2 3 4 5 6 7 8 9; do
  l=`$samtools view $1 1:3,000,000-3,010,000 2>&1 | grep truncate`
  if [ "$l" = "" ]; then
    md5sum $1               # Let caller trap results
    exit $?
  fi
  echo "Try $n : File '$1' is not complete, waiting a little bit to try again" >&2
  sleep 60
done
echo "File '$1' was always truncated, giving up" >&2
exit 3
