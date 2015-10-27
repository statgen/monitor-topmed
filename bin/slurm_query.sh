#!/bin/bash
#
#   slurm_query.sh -squeue partition
#
#   Run squeue -p partition
#	This is needed because too many machines in our cluster cannot
#   run squeue
#
url=http://topmed:48109
tmpfile=/tmp/showslurm.out
lastn=25

if [ "$2" = "" ]; then
  echo "Usage: $0 -squeue partition"
  echo "Usage: $0 -df     path-is-ignored"
  echo "Usage: $0 -logs   path-is-ignored"
  echo "Usage: $0 -errorcheck  path-is-ignored"
  echo "Usage: $0 -restart password"
  echo "Usage: $0 -stop    password"
  exit 1
fi

if [ "$1" = "-df" ]; then
  /usr/bin/wget -o /dev/null -O $tmpfile $url/df
  if [ "$?" != "0" ]; then
    echo "<b>Query for disk usage failed. Perhaps daemon at '$url' is not running</b>"
    exit 2
  fi
  echo "<pre>"
  cat $tmpfile
  echo "</pre>"
  rm -f $tmpfile
  exit
fi

if [ "$1" = "-errorcheck" ]; then
  /usr/bin/wget -o /dev/null -O $tmpfile $url/errorcheck
  if [ "$?" != "0" ]; then
    echo "<b>Unable to run topmed_monitor.pl check. Perhaps daemon at '$url' is not running</b>"
    exit 2
  fi
  echo "<pre>"
  cat $tmpfile
  echo "</pre>"
  rm -f $tmpfile
  exit
fi

if [ "$1" = "-logs" ]; then
  /usr/bin/wget -o /dev/null -O $tmpfile $url/logs
  if [ "$?" != "0" ]; then
    echo "<b>Query for log information failed. Perhaps daemon at '$url' is not running</b>"
    exit 2
  fi
  echo "<pre>"
  cat $tmpfile
  echo "</pre>"
  rm -f $tmpfile
  exit
fi

if [ "$1" = "-restart" ]; then
  /usr/bin/wget -o /dev/null -O $tmpfile $url/restart/$2
  if [ "$?" != "0" ]; then
    echo "<b>Restart of daemon failed. Perhaps daemon at '$url' is not running</b>"
    exit 2
  fi
  echo "<pre>"
  cat $tmpfile
  echo "</pre>"
  rm -f $tmpfile
  exit
fi

if [ "$1" = "-stop" ]; then
  /usr/bin/wget -o /dev/null -O $tmpfile $url/stop/$2
  if [ "$?" != "0" ]; then
    echo "<b>Stop of daemon failed. Perhaps daemon at '$url' is not running</b>"
    exit 2
  fi
  echo "<pre>"
  cat $tmpfile
  echo "</pre>"
  rm -f $tmpfile
  exit
fi

if [ "$1" = "-squeue" ]; then
  /usr/bin/wget -o /dev/null -O $tmpfile $url/squeue/$2
  if [ "$?" != "0" ]; then
    echo "<b>Query for partition '$2' failed. Perhaps daemon at '$url' is not running</b>"
    exit 2
  fi
  n=`cat $tmpfile | wc -l`
  n=`expr $n - 1`
  if [ "$n" -lt "0" ]; then
    n=0
    echo "<b>Unable to get details for partition '$2'</b>"
    exit
  fi
  echo "<p>Partition <b>$2</b> has $n jobs queued"
  if [ "$2" = "nomosix" ]; then
    grep topmed $tmpfile > $tmpfile.tmp         # Only look at topmed jobs
    mv $tmpfile.tmp $tmpfile
    for t in ver bac bai qpl cra ncb; do
      s=`grep $t $tmpfile|wc -l`
      r=`grep $t $tmpfile|grep ' R '|wc -l`
      if [ "$s" != "0" ]; then
        echo "<br/>$t jobs: $s   ($r running)"
      fi
    done
  fi
  if [ "$2" = "topmed-incoming" -o "$2" = "topmed2-incoming" ]; then
    for t in ver bac bai qpl cra ncb; do
      s=`grep $t $tmpfile|wc -l`
      r=`grep $t $tmpfile|grep ' R '|wc -l`
      if [ "$s" != "0" ]; then
        echo "<br/>$t jobs: $s   ($r running)"
      fi
    done
  fi
  echo "<br/>Last $lastn queued jobs are:</p><pre>"
  grep backup $tmpfile | tail -$lastn $tmpfile
  echo "</pre>"
  rm -f $tmpfile
  exit
fi

echo "Unknown option '$1'"
exit 3
