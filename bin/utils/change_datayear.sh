#!/bin/bash
#
#  Usage:  change_datayear.sh runid newval
#
dosql=/usr/cluster/boehnke/bin/dosql.pl 
topmedcmd=/usr/cluster/topmed/bin/topmedcmd.pl

if [ "$2" = "" ]; then
  echo "$0 runid newval"
  echo ""
  echo "Change the datayear for a run"
  echo "Warning - no checking is done here"
  exit 1
fi

$dosql -cmd "UPDATE runs SET datayear=$2 WHERE runid=$1"
if [ "$?" != "0" ]; then
  echo "Unable to set datayear for run '$1'"
  exit 2
fi

$dosql -cmd "UPDATE bamfiles SET datayear=$2 WHERE runid=$1"
if [ "$?" != "0" ]; then
  echo "Unable to set datayear for all samples in run '$1'"
  exit 2
fi
