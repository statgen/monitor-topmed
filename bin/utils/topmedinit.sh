#!/bin/bash
#
#	topmedinit.sh  center
#
topmedcmd=/usr/cluster/topmed/bin/topmedcmd.pl

if [ "$1" = "" ]; then
  echo "$0  center    Force monitor to look for newly arrived runs"
  echo ""
  echo "This can be very slow -- as it will access ALL incoming"
  echo "directories on all topmed nodes"
  echo "In each it looks for center/* runs which are not known"
  echo "It will then ask you if this run should be linked into topmed:/incoming/center"
  echo "and then the monitor will attempt to find samples in each"
  exit 1
fi
center=$1

declare -A runs
for r in `$topmedcmd list runs $center`; do
  runs[$r]=1
done
if [ "${#runs[@]}" = "0" ]; then
  echo "$now Center '$center' has no runs - really?"
  exit 2
fi

#   Look for a run we do not know about yet
#   ddn added at front because the same run will show up in topmed9, we just want ddn run
n=0
for h in ddn topmed topmed2 topmed3 topmed4 topmed5 topmed6 topmed7 topmed9 topmed10; do
  cd /net/$h/incoming/topmed/$center 2> /dev/null
  if [ "$?" != "0" ]; then
    continue
  fi
  for r in `ls .`; do
    if [ ! -d $r ]; then            # If we cannot read directory, it is not local
      #echo "$r is not a local directory"
      continue;
    fi
    #   Check if this is a known run for this senter
    if [ "${#runs[$r]}" = "1" ]; then
      #echo "$r is known, ignoring"
      continue;
    fi
    #   This is not something I know about
    echo -n "Center=$center  run='$r'  Shall we set this new run up? (y/n): "
    read a
    if [ "$a" = "y" ]; then
      echo "Setting up $center/$r"
      n=`expr $n + 1`
      newdir=$h/incoming/topmed/$center/$r
      cd /net/topmed/incoming/topmed/$center || exit 4
      ln -s ../../../../$newdir . || exit 3
      cd /net/$h/incoming/topmed/$center
      runs[$r]=1           # Now this center is known
    fi
  done
done

if [ "$n" = "0" ]; then
  echo -n "I found no new runs, you still want to run topmed_init.pl ? (y/n): "
  read a
  if [ "$a" = "y" ]; then
    n=1
  fi
fi

if [ "$n" != "0" ]; then
  cd $HOME
  /usr/cluster/topmed/bin/topmed_init.pl -center $center updatedb
fi
exit $?

