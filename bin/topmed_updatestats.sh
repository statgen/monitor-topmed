#!/bin/bash
#
#   topmed_updatestats.sh [yyyy mm dd]
#
#	Update stats once a day
#   Because jobs might take a long time to complete, we must redo
#   the previous day or two, just to be sure we get all the stats
#
yyyy=`date +%Y`
mm=`date +%m`
dd=`date +%d`
if [ "$3" != "" ]; then
  yyyy=$1
  mm=$2
  dd=$3
fi

d="$yyyy/$mm/$dd"
range=`date +%Y/%m/%d`
x0=`date +%Y/%m/%d --date "$d"`
x1=`date +%Y/%m/%d --date "$d -1 day"`
x2=`date +%Y/%m/%d --date "$d -2 day"`
x3=`date +%Y/%m/%d --date "$d -3 day"`
range="$x3 $x2 $x1 $x0"

for yyyymmdd in $range; do
  /usr/cluster/monitor/bin/topmedstats.pl jobid   $yyyymmdd
  /usr/cluster/monitor/bin/topmedstats.pl summary $yyyymmdd
done
