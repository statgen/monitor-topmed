#!/bin/bash
#
#   topmedforcew.sh dir
#
#	Force in incoming directory to be writable so the monitor can deal with it
#   This depends on a file set by topmednewrun.sh   /tmp/run.forcew
#
#   08 2,9,14,19 * * * /usr/cluster/topmed/bin/topmedforcew.sh
topuser=topmed
topgrp=topmed
topdir=/incoming/topmed

for f in `ls /tmp/*.forcew 2> /dev/null`; do
  path=`cat $f`
  cd $topdir/$path
  if [ "$?" = "0" ]; then
    echo "Setting ownership and permissions for $path"
    chown -R topmed .
    chgrp topmed .
    chmod 770 .
    chmod 660 Manifest.txt *.md5 *.cram *.crai *.bam *.bai 2> /dev/null   # Failure is OK
  fi
  rm $f
done
