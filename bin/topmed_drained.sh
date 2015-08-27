#!/bin/bash
#
#   topmed_drained.sh
#
#	SLURM gets drained too frequently so we check it's state
#   every now and then and resume it
#   This script must be run by a user with SLURM admin rights
#
host=topmed
state=`/usr/cluster/bin/sstate | grep $host | grep DRAINED`
if [ "$state" = "DRAINED" ]; then
  echo "Resuming DRAINED SLURM daemon"
  scontrol update NodeName=$host State=RESUME
fi
