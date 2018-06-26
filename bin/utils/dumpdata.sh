#!/bin/bash
#
#  dumpdata.sh  [outcsvfile]
#	Dump the bamfiles table
#
table=states
table2=bamfiles
ofile=/tmp/$table2.csv
if [ "$1" != "" ]; then
  ofile=$1
fi
/usr/cluster/boehnke/bin/sql2csv.pl -re topmed --sqlcmd "select * from $table order by id" -hdr /tmp/states.csv
/usr/cluster/boehnke/bin/sql2csv.pl -re topmed --dropcols emsg --sqlcmd "select * from $table2" -hdr  $ofile

exit $?

