#!/bin/bash
#
#   topmed_gcecleanup.sh -submit bamid
#
#	Remove all backups of the original cram file
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc

#------------------------------------------------------------------
# Subroutine:
#   ipath=`B37ValidateIndex(bamid,file)`
#
#   file must exist
#
#	print indexfile_path valid_or_not
#------------------------------------------------------------------
function B37ValidateIndex {
  local bamid=$1
  local file=$2         # BAM or CRAM

  #  To verify an index, we search for a string - depending on the extension
  ext="${file##*.}"
  buildstr='notset'
  indexfile='notset'
  if [ "$ext" = "bam" ]; then
    indexfile=$file.bai
    buildstr="3:148,100,000-148,110,000"
  fi
  if [ "$ext" = "cram" ]; then
    indexfile=$file.crai
    build=`GetBuild $bamid $file`
    if [ "$build" = "37" ]; then
      buildstr="hs37d5:35,450,000-35,477,943"
      if [ "$project" = "inpsyght" ]; then  # Broad use differ 37 ref for this project
        buildstr="GL000192.1:544,000-547,490"
      fi
    fi
    if [ "$build" = "38" ]; then
      buildstr="chr19_KI270938v1_alt"
    fi
  fi
  if [ "buildstr" = 'notset' ]; then
    Fail "ValidateIndex: Unable to determine string to search for in header of $file"
  fi

  #   Index must exist and not be null
  if [ ! -f $indexfile ]; then
    echo "$indexfile invalid"
    return
  fi
  if [ ! -s $indexfile ]; then
    echo "$indexfile invalid"
    return
  fi

  #   Index exists, let's see if it is any good
  n=`$samtools view $file $buildstr | wc -l`
  if [ "$n" = "" -o "$n" -lt "100" ]; then
    #   This might be a trashed reference index for samtools, if so remove it
    if [ "$me" != "" ]; then
      consfile="$console/$bamid-$me.out"
      a=`grep 'cram_ref_load: Assertion' $consfile`     # Could fail if $me not defined
      if [ "$a" != "" ]; then
        rm -rf $HOME/.cache/hts-ref/*/*
         >&2 echo "Samtools cram_ref_load error, removed dirty reference cache and retrying"
      fi
    fi 
    echo "$indexfile invalid"
    return
  fi
  echo "$indexfile valid"
  return
}

me=gcecleanup
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  SubmitJob $bamid $PROJECT '1G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Remove all backups of the original cram file"
  exit 1
fi
bamid=$1
cramfile=`$topmedpath wherefile $bamid b38`
nwdid=`GetNWDID $bamid`
bamid=`GetDB $nwdid bamid`
build=`GetDB $nwdid build`

Started
stime=`date +%s`

if [ "$build" = "37" ]; then
  x=(`B37ValidateIndex $bamid $cramfile`)
else
  x=(`ValidateIndex $bamid $cramfile`)
fi
if [ "${x[1]}" != "valid" ]; then
  Fail "Unable to validate $cramfile - $x"
fi
echo "File '$cramfile' is valid"

#   Now figure out if there is something to remove in a bucket
run=`$topmedcmd show $bamid run`
center=`$topmedcmd show $bamid center`
if [ "$run" = "" -o "$center" = "" ]; then
  Fail "Unable to figure out run or center for $bamid/$nwdid"
fi
uri=''
echo ""
echo "Searching for a backup of sample '$nwdid' in run '$run' center '$center'"
$gsutil ls gs://topmed-backups/$center/$run/$nwdid.src.cram 2>/dev/null 
if [ "$?" = "0" ]; then
  uri=gs://topmed-backups/$center/$run
  echo "Removing files from $uri"
  $gsutil rm $uri/$nwdid.src.cram $uri/$nwdid.src.cram.crai
  echo ""
fi

$gsutil ls gs://topmed-archives/$center/$run/$nwdid,src,cram 2>/dev/null
if [ "$?" = "0" ]; then
  uri=gs://topmed-archives/$center/$run
  echo "Removing files from $uri"
  $gsutil rm $uri/$nwdid.src.cram $uri/$nwdid.src.cram.crai
  echo ""
fi

$gsutil ls gs://topmed-irc-working/archives/$center/$run/$nwdid.src.cram 2>/dev/null
if [ "$?" = "0" ]; then
  uri=gs://topmed-irc-working/archives/$center/$run
  echo "Removing files from $uri"
  $gsutil rm $uri/$nwdid.src.cram $uri/$nwdid.src.cram.crai
  echo ""
fi

if [ "$uri" = "" ]; then
  echo "No backup 
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Clean backups for '$cramfile' [$nwdid] finished in $etime seconds"
Successful
Log $etime
