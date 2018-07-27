#!/bin/bash
#
#   topmed_gcecleanup.sh -submit bamid
#
#	Remove all backups of the original cram file
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc

#------------------------------------------------------------------
# Subroutine:
#   uri=`WhereIsGCEBackup nwdid`
#
#   Requires cache files to be created first:
#       gsutil ls -l -r gs://topmed-backups/ > topmed-backups.cache.txt
#       gsutil ls -l -r gs://topmed-archives/ > topmed-archives.cache.txt
#       gsutil ls -l -r gs://topmed-irc-working/archives > topmed-irc-working.cache.txt
#
#   prints string of:
#       size of cram/bam  uri to cram/bam   list of uris to delete
#------------------------------------------------------------------
function WhereIsGCEBackup {
  local nwd=$1
  tmp=/tmp/WhereIsGCEBackup.$$
  grep $nwdid /net/topmed/working/home/topmed/*.cache.txt | sed -e 's/.txt:/ /' > $tmp
  if [ "$?" != "0" ]; then
    return
  fi
  #    Peel out interesting bits from cache file
  srcfile=`grep -e 'bam$' $tmp`
  ext=bam
  if [ "$srcfile" = "" ]; then
    srcfile=`grep -e 'cram$' $tmp`
    ext=cram
  fi
  if [ "$srcfile" = "" ]; then
    echo none
    return
  fi
  srcsize=`echo $srcfile | awk '{print $2}'`
  srcfile=`echo $srcfile | awk '{print $4}'`
  rmlist=''
  for f in `awk '{print $4}' < $tmp`; do rmlist="$rmlist $f"; done
  echo "$ext $srcfile $srcsize $rmlist"
  rm -f $tmp
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
nwdid=`GetNWDID $bamid`
bamid=`GetDB $nwdid bamid`

Started
stime=`date +%s`

#   Insure we really have local store cram
cramfile=`$topmedpath wherefile $bamid cram`
if [ ! -f "$cramfile" ]; then
  Fail "Until to find cram for $nwdid [ $bamid] - $cramfile"
fi

#   Get backup files from GCE
x=(`WhereIsGCEBackup $nwdid`)
gceext=${x[0]}
gcefile=${x[1]}
gcesize=${x[2]}
unset x[0]
unset x[1]
unset x[2]
if [ "$gceext" = "none" ]; then
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "Nothing removed from GCE backups [$nwdid] in $etime seconds"
  Successful
  Log $etime
  exit
fi

#   We have something in GCE backups. See if local file size matches
file=`$topmedpath wherefile $bamid $gceext`
if [ "$file" = "" ]; then
  Fail "Unable to find local path to $gceext for $nwdid [ $bamid ]"
fi
sz=`stat --dereference --format %s $file 2>/dev/null`
if [ "$?" != "0" ]; then            # If local file was removed, cannot compare
    if [ "$gceext" = "cram" ]; then
      Fail "Local cram does not exist - how can that be?"
    fi
    echo "No local file for $gceext, just removing backup files"
    $gsutil rm ${x[*]}         # Remove all GCE files
    etime=`date +%s`
    etime=`expr $etime - $stime`
    echo "No local $gceext file, removed GCE backups [$nwdid] in $etime seconds"
    Successful
    Log $etime
    exit
fi

#   We have a local file, see if it's size matches that in backups
if [ "$sz" != "$gcesize" ]; then
  Fail "Size mismatch: $gcefile / $gcesize != $file / $sz"
fi
echo "Local $file size matches that in backups"
$gsutil rm ${x[*]}         # Remove all GCE files

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Removed GCE backups [$nwdid] in $etime seconds"
Successful
Log $etime
exit


#------------------------------------------------------------------
# Subroutine:
#   ipath=`B37ValidateIndex(bamid,file)`
#
#   file must exist
#
#	print indexfile_path valid_or_not
#
########### This function is apparently not working, trying something else ##########
#------------------------------------------------------------------
function B37ValidateIndex {
  local bamid=$1
  local file=$2         # BAM or CRAM

#samtools view -H /net/topmed5/working/backups/incoming/topmed/nygc/nhlbi-umich-20150904_1/NWD979839.src.cram | awk -f /net/topmed/incoming/study.reference/current.6.2016/nhlbi.1530.cram.end.region.awk

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
    build=37
    if [ "$build" = "37" ]; then
      buildstr=`samtools view -H $file | awk -f /net/topmed/incoming/study.reference/current.6.2016/nhlbi.1530.cram.end.region.awk`
      echo $buildstr > /tmp/j
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
    echo "$indexfile invalid $buildstr"
    return
  fi
  if [ ! -s $indexfile ]; then
    echo "$indexfile invalid $buildstr"
    return
  fi

  #   Index exists, let's see if it is any good
  n=`$samtools view $file $buildstr | wc -l`
  if [ "$n" = "" -o "$n" -lt "200" ]; then
    echo "$indexfile invalid $buildstr $n"
    return
  fi
  echo "$indexfile valid $buildstr $n"
  return
}

#------------------------------------------------------------------
# Subroutine:
#   uri=`xxWhereIsGCEBackup nwdid`
#
#   prints uri sizeoffile or nothing
#------------------------------------------------------------------
function xxWhereIsGCEBackup {
  local nwd=$1
  run=`$topmedcmd show $nwd run`
  center=`$topmedcmd show $nwd center`
  if [ "$run" = "" -o "$center" = "" ]; then
    Fail "Unable to figure out run or center for $nwd"
  fi
  >&2 echo "Searching for a backup of sample '$nwdid' in run '$run' center '$center'"
  x=(`$gsutil ls -l gs://topmed-backups/$center/$run/$nwdid.src.cram 2>/dev/null | grep $nwd`)
  if [ "${x[0]}" != "" ]; then
    print ${x[0]} gs://topmed-backups/$center/$run
    return
  fi
  x=(`$gsutil ls -l gs://topmed-archives/\*/$nwdid 2>/dev/null | grep $nwd`)
  if [ "${x[0]}" != "" ]; then
    print ${x[0]} gs://topmed-archives/$center/$run
    return
  fi
  x=(`$gsutil ls -l gs://topmed-irc-working/archives/$center/$run/$nwdid.src.cram 2>/dev/null | grep $nwd`)
  if [ "${x[0]}" != "" ]; then
    print ${x[0]} gs://topmed-irc-working/archives/$center/$run
    return
  fi
  return
}
