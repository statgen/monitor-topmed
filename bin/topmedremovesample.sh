#!/bin/bash
#
#  Usage:  topmedremovesample.sh [-noprompt] bamid|nwdid
#
#   Use this to completely remove a sample from the database
#   and local storage. Be careful!
dosql=/usr/cluster/boehnke/bin/dosql.pl
topmedcmd=/usr/cluster/topmed/bin/topmedcmd.pl
topmedpath=/usr/cluster/topmed/bin/topmedpath.pl
prompt=y
rmremote=n

#------------------------------------------------------------------
# Subroutine:
#   GetFile bamid type    Prints path to a type of file that exists or null
#------------------------------------------------------------------
function GetFile {
  local bamid=$1
  local type=$2
  f=`$topmedpath wherefile $bamid $type`
  if [ -f "$f" ]; then
    echo $f
  else
    echo ""
  fi
}

#------------------------------------------------------------------
# Subroutine:
#   xsql cmd    Execute dosql.pl and clean up the output
#------------------------------------------------------------------
function xsql {
  q="'"
  eval $dosql -batch -cmd ${q}$*${q} | awk '{if(NR>5)print}' | grep -v 'rows in set' | grep -v Bye
}

if [ "$1" = "-noprompt" ]; then     # Prompt or just blast through
  prompt=n
  shift
fi
if [ "$1" = "-remote" ]; then       # Remove remote storage?
  rmremote=y
  shift
fi

if [ "$1" = "" ]; then
  echo "$0 [-noprompt] [-remote] bamid|nwdid"
  echo ""
  echo "Completely remove a sample from the database and local storage"
  echo "Warning - once started, there is no recovery"
  exit 1
fi

#   Show details of this sample
echo "=================== $1 ==================="
l=`$topmedcmd whatnwdid $1`
if [ "$l" = "" ]; then
  echo "Sample '$bamid' is unknown"
  exit 1
fi
echo $l
l=($l)
nwdid=${l[0]}
bamid=${l[1]}
run=${l[7]}
piname=${l[9]}
center=${l[12]}
year=${l[14]}

#   Show details of run for this sample
l=`$topmedcmd whatrun $bamid`
echo $l
l=($l)
dirname=${l[6]}
count=`expr ${l[10]} - 1`

#   Get the list of files to remove
bam=`GetFile $bamid bam`
if [ "$bam" != "" ]; then
  extension=${bam##*.}
  if [ "$extension" = "bam" ]; then
    bam="$bam $bam.bai"
  else
    bam=''
  fi
fi
bamdir=`$topmedpath wherepath $bamid bam`
manifest=''
if [ -f "$bamdir/Manifest.txt" ]; then
  manifest=$bamdir/Manifest.txt
fi  
bam=`GetFile $bamid bam`        # Actually original incoming file (no longer a bam)
if [ "$bam" != "" ]; then
  bam="$bam $bam.crai $bam.md5"
fi
cram=`GetFile $bamid cram`
if [ "$cram" != "" ]; then
  cram="$cram $cram.crai $cram.md5"
fi
qcresults=`GetFile $bamid qcresults`
if [ "$qcresults" != "" ]; then
  f=`dirname $qcresults` 
  qcresults="$f/$nwdid.*" 
fi
bcf=`GetFile $bamid bcf`
if [ "$bcf" != "" ]; then
  bcf="$bcf $bcf.csi"
fi
state=`$topmedcmd show $bamid state_b37`
b37=''
if [ "$state" = "20" ]; then
  b37=`GetFile $bamid b37`
  if [ "$b37" != "" ]; then
    f=`dirname $b37`
    b37=`dirname $f`
  fi
fi
state=`$topmedcmd show $bamid state_b38`
b38=''
if [ "$state" = "20" ]; then
  b38=`GetFile $bamid b38`
  if [ "$b38" != "" ]; then
    b38=`dirname $b38`
  fi
fi

#   Prompt to see if we should really do this
echo ""
echo "Delete data for sample $bamid $nwdid:"
if [ "$bam"       != "" ]; then echo "  BAM: $bam"; fi
if [ "$manifest"  != "" ]; then echo "  MANIFEST: $manifest"; fi
if [ "$cram"      != "" ]; then echo "  CRAM: $cram"; fi
if [ "$qcresults" != "" ]; then echo "  QC RESULTS: $qcresults"; fi
if [ "$bcf"       != "" ]; then echo "  BCF: $bcf"; fi
if [ "$b37"       != "" ]; then echo "  B37 CRAM: $b37"; fi
if [ "$b38"       != "" ]; then echo "  B38 CRAM: $b38"; fi
echo "Count for run $run will be set to '$count'"
echo ""
if [ "$prompt" = "y" ]; then
  echo -n "Shall we delete data for sample $bamid $nwdid? (y/n): "
  read a
  if [ "$a" != "y" ]; then
    exit 1
  fi
fi

#   Actually remove files
echo "Removing files ..."
rm -f $bam $cram $bcf || exit 2
rm -rf $b37 $b38 || exit 3
if [ "$manifest" != "" ]; then
  echo "Fix Manifest ..."
  if [ ! -f $manifest.orig ]; then
    cp -p $manifest $manifest.orig
  fi
  grep -v $nwdid $manifest > $manifest.xx
  mv $manifest.xx $manifest
fi
echo "Update database ..."
xsql "delete from qc_results where bam_id=$bamid" || exit 4
xsql "delete from bamfiles where bamid=$bamid" || exit 5
xsql "update runs set count=$count where dirname=\"$run\"" || exit 6

echo "Successfully removed sample $bamid/$nwdid"
exit

