#!/bin/bash
#
#   topmed_gcecpbcf.sh -submit| bamid
#
#	Copy bcf files from local storage to a GCE bucket
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
me=gcecpbcf
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='2:00:00'
  SubmitJob $bamid "$PROJECT-gce" '3G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Copy BCF files from local storage to a GCE bucket"
  exit 1
fi
bamid=$1
nwdid=`GetNWDID $bamid`
bamid=`GetDB $nwdid bamid`

Started
stime=`date +%s`

#   Get bcf file for remapped cram file
recabbcf=`$topmedpath wherefile $bamid bcf`
if [ "$recabbcf" = "" ]; then
  Fail "Unable to determine BCF file for '$bamid'"
fi
if [ ! -f $recabbcf ]; then
  Fail "BCF file for '$bamid' does not exist: $recabbcf"
fi

CreateBCFIndex $bamid $recabbcf
recabcsi=$recabbcf.csi

#======================================================================
#   All the files of interest exist locally
#   Copy is very fast to GCE, so don't try to be too clever
#======================================================================
copyuri=`$topmedpath wherepath $nwdid gcebcfupload`
lockdir=/run/shm/$$.$nwdid
mkdir -p $lockdir
gsutil="$gsutil -o GSUtil:state_dir=$lockdir"
$gsutil cp $recabbcf $recabcsi $copyuri 
if [ "$?" != "0" ]; then
  Fail "Failed to copy files to GCE as $gcefile"
fi

#   Clean out files we no longer want in BCF
$gsutil rm $copyuri/$nwdid.recab.cram.md5 $copyuri/$nwdid.recab.cram.flagstat $copyuri/$nwdid.metadata.txt

echo "GCE BCF files for $bamid $nwdid are up to date"
$gsutil ls -l $copyuri

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of BCF files to GCE completed in $etime seconds"
Successful
rm -rf $lockdir
Log $etime
exit
