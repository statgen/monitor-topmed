#!/bin/bash
#
#   topmedforcew.sh 
#
#   Look for new runs on a host and force in incoming directory writable
#   This scripts runs as root, but executes some commands as user topmed
#
#   Making this work requires two users - root to change ownership and
#   the topmed user to be able to make changes to the NFS mounted data
#   This runs as root, but uses commands like 'su topmed -c somecommand'
#   for the few cases where topmed authorizatio is needed.
#
topuser=topmed
topgrp=topmed
topdir=/incoming/topmed
topmedcmd=/usr/cluster/topmed/bin/topmedcmd.pl
export PROJECT=topmed
myhost=`hostname`               # Where I am running

#------------------------------------------------------------------
# Subroutine:
#   ProcessRuns(center)
#
#   Look for new runs for this center on this hose
#------------------------------------------------------------------
function ProcessRuns {
  local center=$1
  now=`date '+%Y/%m/%d %H:%M'`
  #   Get list of all known runs for this center
  declare -A runs
  for r in `$topmedcmd list runs $center`; do
  runs[$r]=1
  done
  if [ "${#runs[@]}" = "0" ]; then
    echo "$now Center '$center' has no runs - really?"
    return 2
  fi

  #   Look for a run we do not know about yet
  if [ ! -d $topdir/$center ]; then
    return 3
  fi
  cd $topdir/$center || exit 2

  n=0
  for r in `ls .`; do
    if [ ! -d $r ]; then            # If we cannot read directory, it is not local
      continue;
    fi
    #   Manifest must exist
    owner=$(stat -c '%U' $r/Manifest.txt 2>/dev/null)
    if [ "$owner" = "" ]; then
      echo "$now No Manifest.txt exists for $center/$r yet"
      continue;
    fi

    #   Run must exist, not owned by $topuser
    if [ "$owner" != "$topuser" ]; then
      newrun=/net/topmed/incoming/topmed/$center/$r
      echo -n "$now Monitor knows about ${#runs[@]} runs for $center  "
      # Determine if this is year 2 (bam) or 3 (cram) - determines where the backup dir is
      #year=0
      #if [ -f $r/Manifest.txt ]; then
      #  s=`head -1 $r/Manifest.txt grep cram 2> /dev/null`
      #  if [ "$s" != "" ]; then year=3; fi
      #  s=`head -1 $r/Manifest.txt grep bam 2> /dev/null`
      #  if [ "$s" != "" ]; then year=2; fi
      #fi
      #  With luck we have determined if this year 2 or 3
      #if [ "$year" = "0" ]; then
      #  echo "$now Cannot figure out the year for '$r'"
      #  ls -l $r/Manifest.txt
      #  continue
      #fi
      #   We have a new run, set it up
      date
      echo "$now New run '$r' year $year for '$center' found"
      if [ "$myhost" != "topmed" ]; then
        su $topuser -c "ln -s ../../../../$myhost/incoming/topmed/$center/$r $newrun"
        echo "$now Set symlink to $newrun"
        backuphost=$myhost          # Year 3 uses the incoming directory for 'backup'
        masterbackuprun="../../../../../../$backuphost/incoming/topmed/$center/$r"
      fi
      #     Sort out backup directory for run
      #if [ "$year" = "2" ]; then
      #  backuphost=topmed7          # Year 2 we make a directory elsewhere for backups
      #  su $topuser -c "mkdir -p /net/$backuphost/working/backups/incoming/topmed/$center/$r"
      #  masterbackuprun="../../../../../../$backuphost/working/backups/incoming/topmed/$center/$r"
      #fi
      #if [ "$year" = "3" ]; then
      #  backuphost=$myhost          # Year 3 uses the incoming directory for 'backup'
      #  masterbackuprun="../../../../../../$backuphost/incoming/topmed/$center/$r"
      #fi

      newbackuprun=/net/topmed/working/backups/incoming/topmed/$center/$r
      su $topuser -c "ln -s $masterbackuprun $newbackuprun"
      echo "$now BACKUP Symlink set to host '$backuphost' for $newbackuprun ($masterbackuprun)"
      n=`expr $n + 1`
      #     Force ownership of run and files since we are root
      chown -R $topuser $r
      chgrp -R $topgrp $r
      chmod 770 $r
      chmod 660 $r/Manifest.txt $r/*.cram $r/*.bam 2> /dev/null
      echo "$now Ownership and permissions set for $center/$r"
    fi
  done

  #  Because ownership of files still might be wrong, we must find all
  #  files in all directories that are not owned by topmed and force ownership
  for r in `ls .`; do
    if [ ! -d $r ]; then            # If we cannot read directory, it is not local
      continue;
    fi
    #   Manifest must exist
    stat $r/Manifest.txt 2>/dev/null >/dev/null
    if [ "$?" != "0" ]; then
      echo "$now No Manifest.txt exists for $center/$r yet"
      continue;
    fi
    #   Any files owned by not-topmed?  If so, force ownership again
    n=`ls -l  $r/ |grep -v 'topmed topmed'|grep -v total|grep -v /:|wc -l`
    if [ "$n" != "0" ]; then
      chown -R $topuser $r 
      chgrp -R $topgrp $r
      chmod 770 $r
      chmod 660 $r/Manifest.txt $r/*.cram $r/*.bam 2> /dev/null
      echo "$now Forced ownership of '$r'"
    fi
  done
  return
}

#=========================================================
#   Main line code
#=========================================================
wait=$((1 + RANDOM % 20))
sleep $wait                     # Wait a little bit to avoid messing up logs

#   Capture output in local file and if there is something append to useful log
locallog=/tmp/topmedforcew.log.$$
ProcessRuns baylor > $locallog
ProcessRuns washu  >> $locallog
if [ -s $locallog ]; then
  su $topuser -c "(cat $locallog >> /net/topmed/working/topmed-output/topmed_init.log)"
  chown topmed $locallog          # So I can remove these
else
  rm $locallog
fi
exit
