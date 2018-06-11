#!/bin/bash
#  Save information on SLURM queues. Used by topmedthrottle and topmedpermit
#f=/run/shm/Slurm.summary
#/usr/cluster/topmed/bin/topmedcluster.pl summary > $f

#  Check the state of SLURM hosts, remember which are down
f=/run/shm/RandomRealHost.down
rm -f $f
for n in 2 3 4 5 6 7 9 10; do
  h=topmed$n
  l=`/usr/cluster/bin/scontrol show node $h | grep State= | sed -e 's/=/ /'`
  ll=($l)
  if [ "${ll[1]}" = "DOWN" -o "${ll[1]}" = "DOWN*" ]; then
    n=`expr $n - 1`
    echo "$n $h is down" >> $f
    #now=`date '+%Y/%m/%d %H:%M'`
    #echo "$now $h is down"
  fi
done
