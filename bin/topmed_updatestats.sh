#!/bin/bash
#
#   topmed_updatestats.sh
#
#	Update stats once a day
#   Because jobs might take a long time to complete, we must redo
#   the previous day or two, just to be sure we get all the stats
#
yyyy=`date +%Y`
mm=`date +%m`
dd=`date +%d`
range=""
if [ "$dd" = "01" -o "$dd" = "02" ]; then
  yyy=$yyyy
  m=`expr $mm - 1`
  if [ "$m" = "0" ]; then
    yyy=`expr $yyyy - 1`
    m="12"
  fi
  if [ "$m" -lt "10" ]; then
    m="0$m"
  fi
  range="$yyy/$m/28 $yyy/$m/29 $yyy/$m/30 $yyy/$m/31 $yyyy/$mm/01 $yyyy/$mm/02"
fi
if [ "$dd" = "10" -o "$dd" = "11" ]; then
  range="$yyyy/$mm/09 $yyyy/$mm/10 $yyyy/$mm/11"
fi
if [ "$range" = "" ]; then
  d=`expr $dd - 2`
  for n in `seq -w $d $dd`; do
    range="$range $yyyy/$mm/$n"
  done
fi
for yyyymmdd in $range; do
  /usr/cluster/monitor/bin/topmedstats.pl jobid   $yyyymmdd
  /usr/cluster/monitor/bin/topmedstats.pl summary $yyyymmdd
done
