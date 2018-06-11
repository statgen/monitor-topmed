#!/bin/bash
#
#  Usage:  where.sh NWDnnnnn|bamid bam|backup|qcresults|console|b37|b38
#
bamid=$1
fcn=bam
if [ "$2" != "" ]; then
  /usr/cluster/topmed/bin/topmedpath.pl wherefile $bamid $2
  exit $? 
fi
#echo -n "FCN=$fcn  "
echo -n 'bam   '; /usr/cluster/topmed/bin/topmedpath.pl wherefile $bamid bam
echo -n 'cram  '; /usr/cluster/topmed/bin/topmedpath.pl wherefile $bamid cram
echo -n 'b38   '; /usr/cluster/topmed/bin/topmedpath.pl wherefile $bamid b38
echo -n 'b37   '; /usr/cluster/topmed/bin/topmedpath.pl wherefile $bamid b37
echo ''
