#!/bin/bash
#
#   topmednewrun.sh
#
#	Look for new runs on a host.  This runs as user topmed
#
#   07 2,9,14,19 * * * /usr/cluster/topmed/bin/topmednewrun.sh >> /net/topmed/working/topmed-output/topmed_init.log
topmedcmd=/usr/cluster/topmed/bin/topmedcmd.pl

HOST=`hostname`

#------------------------------------------------------------------
# Subroutine:
#   ProcessRuns(center)
#
#   Look for new runs for this center on this hose
#------------------------------------------------------------------
function ProcessRuns {
  local center=$1

  #   Get list of all known runs for this center
  declare -A runs
  for r in `$topmedcmd list runs $center`; do
  runs[$r]=1
  done
  if [ "${#runs[@]}" = "0" ]; then
    echo "Center '$center' has no runs - really?"
    return 2
  fi

  #   Look for a run we do not know about yet
  cd /net/$HOST || exit 2
  if [ ! -d incoming/topmed/$center ]; then
    return
  fi
  cd incoming/topmed/$center || exit 2

  n=0
  for r in `ls .`; do
    if [ "${runs[$r]}" = "" ]; then       # Run not known by monitor yet
      newrun=/net/topmed/incoming/topmed/$center/$r
      if [ ! -d $newrun ]; then
        echo -n "Monitor knows about ${#runs[@]} runs for $center      "
        date
        echo "New run '$r' for '$center' detected"
        ln -s ../../../../$HOST/incoming/topmed/$center/$r $newrun
        echo $center/$r > /tmp/$r.forcew  # Used by topmedforcew.sh to set ownership
        # Determine if this is year 2 (bam) or 3 (cram) - determines where the backup dir is
        year=0
        if [ -f $newrun/Manifest.txt ]; then
          s=`grep cram $newrun/Manifest.txt 2> /dev/null`
          if [ "$s" != "" ]; then year=3; fi
          s=`grep bam $newrun/Manifest.txt 2> /dev/null`
          if [ "$s" != "" ]; then year=2; fi
        fi
        if [ "$year" = "0" ]; then        # No Manifest yet, look for file
          s=`ls $newrun/*.cram 2> /dev/null`
          if [ "$s" != "" ]; then year=3; fi
          s=`ls $newrun/*.bam 2> /dev/null`
          if [ "$s" != "" ]; then year=2; fi
        fi
        #  With any luck we have determined if this year 2 or 3
        if [ "$year" = "0" ]; then
          echo "Cannot set the backup directory because I cannot figure out the year for '$r'"
        else
          backuphost=$HOST                # Year 3 uses the incoming directory for 'backup'
          if [ "$year" = "2" ]; then      # For now, force all backups for year 2 here
            backuphost=topmed7
            mkdir -p /net/$backuphost/working/backups/incoming/topmed/$center/$r
          fi
          newbackuprun=/net/topmed/working/backups/incoming/topmed/$center/$r
          ln -s ../../../../../../$backuphost/incoming/topmed/$center/$r $newbackuprun
          echo "Symlinks set up for '$r' for monitor processes"
          n=`expr $n + 1`
        fi
      fi
    fi
  done
#  if [ "$n" = "0" ]; then
#    echo "No new runs detected for $HOST $center"
#    return 2
#  fi
  return
}

#=========================================================
#   Main line code
#=========================================================
if [ ! -d /net/$HOST/incoming/topmed/$center ]; then exit; fi   # Center not on this host

ProcessRuns baylor
ProcessRuns washu
exit


