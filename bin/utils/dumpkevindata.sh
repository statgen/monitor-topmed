#!/bin/bash
#
#  dumpdata.sh  [outcsvfile]
#	Dump the qc_results table
#
table=qc_results
ofile=/tmp/$table.csv
if [ "$1" != "" ]; then
  ofile=$1
fi
/usr/cluster/boehnke/bin/sql2csv.pl -re topmed --sqlcmd "select * from $table" -hdr  $ofile

exit $?

