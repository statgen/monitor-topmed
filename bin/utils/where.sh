#!/bin/bash
#
#  Usage:  where.sh NWDnnnnn|bamid bam|backup|qcresults|console|b37|b38
#
bamid=$1
fcn=bam
pth=/usr/cluster/topmed/bin/topmedpath.pl
if [ "$2" != "" ]; then
  $pth wherefile $bamid $2
  exit $? 
fi
#echo -n "FCN=$fcn  "
echo -n 'bam   '; $pth wherefile $bamid bam
echo -n 'cram  '; $pth wherefile $bamid cram
echo -n 'b38   '; $pth wherefile $bamid b38
echo -n 'b37   '; $pth wherefile $bamid b37
echo ''
