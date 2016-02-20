#!/bin/bash
#
#   topmed_drained.sh hostname
#
#	SLURM gets drained too frequently so we check it's state
#   every now and then and resume it
#   This script must be run by a user with SLURM admin rights
#
hosts=topmed
if [ "$1" != "" ]; then
  hosts="$*"
fi
for host in $hosts; do
  state=`/usr/cluster/bin/sinfo -p ${host}-incoming show node | grep drain`
  if [ "$state" != "" ]; then
    echo "Resuming DRAINED SLURM daemon at '$host' (state=$state)"
    scontrol update nodename=$host State=RESUME
    /usr/cluster/bin/sinfo -p ${host}-incoming show node
  #else 
    #echo "Node '$host' does not need to be resumed"
    #/usr/cluster/bin/sinfo -p ${host}-incoming show node
  fi
done
