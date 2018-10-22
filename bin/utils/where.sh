#!/bin/bash
#
#  Usage:  where.sh NWDnnnnn|bamid bam|backup|qcresults|console|b37|b38
#
bamid=$1
fcn=bam
pth=/usr/cluster/topmed/bin/topmedpath.pl
shift;

#  better handling of arguments 2,3,...  tb 9/7/2018 :

if [ $# -ge 1 ]; then x=$*;
else x="bam cram b38";
fi
for t in $x; do
  printf '%-12s' $t; $pth wherefile $bamid $t
done
