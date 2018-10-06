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
  for f in `grep $ext $tmp | awk '{print $4}'`; do rmlist="$rmlist $f"; done
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
  SubmitJob $bamid $PROJECT '4G' "$0 $*"
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
echo "gceext=$gceext gcefile=$gcefile gcesize=$gcesize x=${x[*]}"
exit 99

#   We have something in GCE backups. See if local file size matches
filepath=`$topmedpath wherepath $bamid $gceext`
if [ "$filepath" = "" ]; then
  Fail "Unable to find local path to $gceext for $nwdid [ $bamid ]"
fi
if [ "$gceext" = "bam" ]; then
  f=`$topmedcmd show $bamid bamname_orig`
else
  f=`$topmedcmd show $bamid cramname`
fi
file=$filepath/$f
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
