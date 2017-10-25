#!/bin/bash
#
#   topmed_redo.sh bamid error message
#     or
#   topmed_redo.sh bamid action args-for-action (based on error message)
#
#   This is designed to redo bits of the monitor processes
#   Ideally it is driven by error messages from topmedcheck.pl -fix topmed_redo.sh
#   Each unique error message will require particular processing to set the arguments properly
. /usr/cluster/topmed/bin/topmed_actions.inc
me=redo
markverb=''

#   If we are passed lots of arguments, this is full message, submit as job
if [ "$#" -gt "5" ]; then
  bamid=`GetDB $1 bamid`
  shift
  RandomRealHost $bamid 'bam'         # Sets random $realhost
  myargs=''
  # Figure out exactly what action to do, only pass relevent arguments to batch job
  if [ "$4 $5" = "must match" ]; then
    myargs="$bamid calcflagstat $2 $6 $7"   # e.g. redo.sh bamid calcflagstat b38flagstat bamflagstat 924109044
    if [ "$7" = "0" ]; then
      echo "No point in calculating flagstats for $2 since $6 = $7"
      exit
    fi
  fi
  if [ "$5 $6 $7" = "Index is older:" ]; then
    myargs="$bamid makeindex $8"            # e.g. redo.sh bamid makindex .../NWD219226-HG5J7ALXX-3.hgv.bam
  fi
  if [ "$2 $3" = "Index missing" ]; then
    myargs="$bamid makeindex $7"
  fi
  if [ "$1 $2 $3 $4" = "GCE recab.cram.flagstat file missing:" ]; then
    myargs="$bamid gceflagstat $5"
  fi
  if [ "$1 $2 $3 $4" = "GCE recab.cram.md5 file missing:" ]; then
    myargs="$bamid gcemd5 $5"
  fi
  if [ "$1 $2 $3 $4" = "GCE recab.cram.crai file missing:" ]; then
    myargs="$bamid gcecrai $5"
  fi
  if [ "$1 $2 $3 $4" = "GCE bcf file missing:" ]; then      # ??? remove
    myargs="$bamid gcebcf $5"
  fi
  if [ "$1 $2 $3 $4" = "GCE bcf.csi file missing:" ]; then
    myargs="$bamid gcecsi $5"
  fi

  if [ "$myargs" != "" ]; then
    echo "Submitting job to topmed-redo($realhost): $0 $myargs"
    SubmitJob $bamid "topmed-redo" '4G' "$0 $myargs"
    exit
  else
    echo "$me - unable to determine exactly what action to redo: $*"
  fi
  exit 2
fi

#-----------------------------------------------------------------------
#   Special argument order here:  bamid action other-args
#-----------------------------------------------------------------------
bamid=$1
action=$2

#   /tmp/redo.sh 91714 Column b38flagstat x must match bamflagstat 0         -> 91714 calcflagstat b38flagstat bamflagstat 0
#   /tmp/redo.sh 51998 Column b38flagstat x must match bamflagstat 755115788 -> 51998 calcflagstat b38flagstat bamflagstat 755115788
#   /tmp/redo.sh 44887 Column b37flagstat x must match bamflagstat 887078342 -> 44887 calcflagstat b37flagstat bamflagstat 887078342
if [ "$action" = "calcflagstat" ]; then
  if [ "$3" = "b38flagstat" ]; then       # Figure out which file we want
    f=`$topmedpath wherefile $bamid b38`
  fi
  if [ "$3" = "b37flagstat" ]; then       # Figure out which file we want
    f=`$topmedpath wherefile $bamid b37`
  fi
  if [ "$f" = "" ]; then
    Fail "$action - Could not determine file: $*"
  fi
  b=`CalcFlagstat $bamid $f`
  if [ "$b" = "" ]; then
    Fail "$action - Unable to calculate flagstat for $bamid f=$f"
  fi
  #  If expected value is 0, pretend all is good, else compare
  if [ "$5" = "0" ]; then
    echo "Setting $3 and $4 to $b for $bamid"
    SetDB $bamid $3 $b
    SetDB $bamid $4 $b
    SetDB $bamid state_gce38copy 1
    exit
  fi
  if [ "$b" = "$5" ]; then
    SetDB $bamid $3 $b
    SetDB $bamid state_gce38copy 1
    echo "Successfully set $3 $b for $bamid"
    exit
  fi
  Fail "$action - flagstat for $f = $b, did not match expected value $4"
fi

#   /tmp/redo.sh 4779  BAM Index is older than path/NWD912733.hg19.bam: path/NWD912733.hg19.bam -> 4779 makeindex path/NWD912733.hg19.bam
#   /tmp/redo.sh 94350 B38 Index missing for a  path/NWD108619.recab.cram: path/NWD108619.recab.cram    (similiar)
if [ "$action" = "makeindex" ]; then
  f=$3
  echo "Creating index for $bamid $f"
  ext="${f##*.}"
  if [ "$ext" = "bam" ]; then
    rm -f $f.bai
  fi
  if [ "$ext" = "cram" ]; then
    rm -f $f.crai
  fi
  CreateIndex $bamid $f
  SetDB $bamid state_gce38copy 1
  SetDB $bamid state_aws38copy 1
  echo "Successfully created index $f for $bamid"
  exit
fi

# /tmp/redo.sh 5805 GCE recab.cram.flagstat file missing: NWD537762.recab.cram.flagstat
if [ "$action" = "gceflagstat" ]; then
  flagstat=`GetDB $bamid b38flagstat`
  nwdid=`GetNWDID $bamid`
  echo "$flagstat + 0 paired in sequencing" > /tmp/$$
  gsutil cp /tmp/$$ gs://topmed-bcf/$nwdid/$3
  rm -f /tmp/$$
  echo "Successfully set flagstat in GCE for $nwdid / $bamid"
  exit
fi

# /tmp/redo.sh 48459 GCE recab.cram.md5 file missing: NWD948036.recab.cram.md5
if [ "$action" = "gcemd5" ]; then
  md5=`GetDB $bamid b38cramchecksum`
  nwdid=`GetNWDID $bamid`
  echo "$md5 $nwdid.recab.cram" > /tmp/$$
  gsutil cp /tmp/$$ gs://topmed-bcf/$nwdid/$3
  rm -f /tmp/$$
  echo "Successfully set md5 in GCE for $nwdid / $bamid"
  exit
fi

# /tmp/redo.sh 4908 GCE recab.cram.crai file missing: NWD498844.recab.cram.crai
if [ "$action" = "gcecrai" ]; then
  cramfile=`$topmedpath wherefile $bamid b38`
  nwdid=`GetNWDID $bamid`
  gsutil cp $cramfile.crai gs://topmed-bcf/$nwdid/$3
  echo "Successfully copied CRAI to GCE for $nwdid / $bamid"
  exit
fi

# /tmp/redo.sh 4908 GCE bcf.csi file missing: NWD498844.bcf.csi
if [ "$action" = "gcecsi" ]; then
  bcffile=`$topmedpath wherefile $bamid bcf`
  nwdid=`GetNWDID $bamid`
  gsutil cp $bcffile.csi gs://topmed-bcf/$nwdid/$3
  echo "Successfully copied CSI to GCE for $nwdid / $bamid"
  exit
fi

# /tmp/redo.sh 4908 GCE bcf file missing: NWD498844.bcf
if [ "$action" = "gcebcf" ]; then
  bcffile=`$topmedpath wherefile $bamid bcf`
  nwdid=`GetNWDID $bamid`
  gsutil cp $bcffile gs://topmed-bcf/$nwdid/$3
  echo "Successfully copied CSI to GCE for $nwdid / $bamid"
  exit
fi

echo "Did not know what to do with action: $action"
exit 3


