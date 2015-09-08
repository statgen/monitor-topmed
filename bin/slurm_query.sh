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
    s=`grep qpl $tmpfile|wc -l`
    r=`grep align-topmed $tmpfile|grep ' R '|wc -l`
    echo "<br/>align jobs: $s   ($r running)"
    s=`grep cram $tmpfile|wc -l`
    r=`grep cram $tmpfile|grep ' R '|wc -l`
    echo "<br/>cram jobs: $s   ($r running)"
    s=`grep backup $tmpfile|wc -l`
    r=`grep backup $tmpfile|grep ' R '|wc -l`
    echo "<br/>backup jobs: $s   ($r running)"
  fi
  if [ "$2" = "topmed" ]; then
    s=`grep bac $tmpfile|wc -l`
    r=`grep bac $tmpfile|grep ' R '|wc -l`
    echo "<br/>backup jobs: $s   ($r running)"
    s=`grep qpl $tmpfile|wc -l`
    r=`grep qpl $tmpfile|grep ' R '|wc -l`
    echo "<br/>qplot jobs: $s   ($r running)"
  fi
  if [ "$2" = "topmed-incoming" -o "$2" = "topmed2-incoming" ]; then
    s=`grep ver $tmpfile|wc -l`
    r=`grep ver $tmpfile|grep ' R '|wc -l`
    echo "<br/>verify jobs: $s   ($r running)"
    s=`grep bac $tmpfile|wc -l`
    r=`grep bac $tmpfile|grep ' R '|wc -l`
    echo "<br/>backup jobs: $s   ($r running)"
    s=`grep bai $tmpfile|wc -l`
    r=`grep bai $tmpfile|grep ' R '|wc -l`
    echo "<br/>bai jobs: $s   ($r running)"
    s=`grep cram $tmpfile|wc -l`
    r=`grep cram $tmpfile|grep ' R '|wc -l`
    echo "<br/>cram jobs: $s   ($r running)"
  fi
  echo "<br/>Last $lastn queued jobs are:</p><pre>"
  grep backup $tmpfile | tail -$lastn
  grep cram $tmpfile | tail -$lastn
  echo "</pre>"
  rm -f $tmpfile
  exit
fi

echo "Unknown option '$1'"
exit 3
