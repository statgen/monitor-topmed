#!/bin/bash
#
#   topmed_drained.sh hostname hostname2 ...
#
#	SLURM gets drained too frequently so we check it's state
#   every now and then and resume it
#   This script must be run by a user with SLURM admin rights
#

for host in $*; do
  state=`/usr/cluster/bin/scontrol show node $host | grep State= | grep DRAIN`
  if [ "$state" != "" ]; then
    #echo "Resuming DRAINED SLURM daemon at '$host'"
    scontrol update nodename=$host state=RESUME 2>/dev/null
  fi
done
