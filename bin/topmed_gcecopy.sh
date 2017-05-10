#!/bin/bash
#
#   topmed_gcecopy.sh -submit| bamid
#
#	Copy files from local storage to a GCE bucket for GCE processing
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcecopy
markverb=$me
cores="--cpus-per-task=2"           # Cores here should be same as gsutil
stat="stat --printf=%s"
bcftools=/usr/cluster/bin/bcftools

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me  $bamid
  RandomRealHost $bamid
  SubmitJob $bamid "topmed-$me" '2G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Copy files from local storage to a GCE bucket"
  exit 1
fi
bamid=$1

Started
GetNWDID $bamid
stime=`date +%s`

#   Get remapped cram file
recabcram=`$topmedpath wherefile $bamid b38`
if [ "$recabcram" = "" ]; then
  Fail "Unable to determine CRAM file for '$bamid'"
fi
if [ ! -f $recabcram ]; then
  Fail "CRAM file for '$bamid' does not exist: $recabcram"
fi
baserecabcram=`basename $recabcram`
sizerecabcram=`$stat $recabcram`

recabcrai=$recabcram.crai
if [ ! -f $recabcrai ]; then
  echo "CRAI missing, trying to create it"
  $samtools index $recabcram
  if [ "$?" != "0" ]; then
    Fail "Unable to create CRAI for $bamid' for: $recabcram"
  fi
fi
if [ ! -f $recabcrai ]; then
  Fail "CRAI file for '$bamid' does not exist: $recabcrai"
fi
baserecabcrai=`basename $recabcrai`
sizerecabcrai=`$stat $recabcrai`

#   Get bcf file for remapped cram file
recabbcf=`$topmedpath wherefile $bamid bcf`
if [ "$recabbcf" = "" ]; then
  Fail "Unable to determine BCF file for '$bamid'"
fi
if [ ! -f $recabbcf ]; then
  Fail "BCF file for '$bamid' does not exist: $recabbcf"
fi
baserecabbcf=`basename $recabbcf`
sizerecabbcf=`$stat $recabbcf`

recabcsi=$recabbcf.csi
if [ ! -f $recabcsi ]; then
  echo "CSI missing, trying to create it"
  $bcftools index $recabbcf
  if [ "$?" != "0" ]; then
    Fail "Unable to create CSI for $bamid' for: $recabbcf"
  fi
fi
if [ ! -f $recabcsi ]; then
  Fail "CSI file for '$bamid' does not exist: $recabcsi"
fi
baserecabcsi=`basename $recabcsi`
sizerecabcsi=`$stat $recabcsi`

#======================================================================
#   All the files of interest exist locally
#   Copy the file to the remote destination only if the size differs
#======================================================================
copylist=''
copyuri=`$topmedpath wherepath $nwdid upload`

l=(`$gsutil stat $copyuri/$baserecabcram | grep Content-Length:`)
gcesize=${l[1]}
#echo "$baserecabcram sizes= $sizerecabcram $gcesize"
if [ $sizerecabcram != $gcesize ]; then
  copylist="$copylist $recabcram"
fi

l=(`$gsutil stat $copyuri/$baserecabcrai | grep Content-Length:`)
gcesize=${l[1]}
#echo "$baserecabcrai sizes= $sizerecabcrai $gcesize"
if [ $sizerecabcrai != $gcesize ]; then
  copylist="$copylist $recabcrai"
fi

l=(`$gsutil stat $copyuri/$baserecabbcf | grep Content-Length:`)
gcesize=${l[1]}
#echo "$baserecabbcf sizes= $sizerecabbcf $gcesize"
if [ $sizerecabbcf != $gcesize ]; then
  copylist="$copylist $recabbcf"
fi

l=(`$gsutil stat $copyuri/$baserecabcsi | grep Content-Length:`)
gcesize=${l[1]}
#echo "$baserecabcsi sizes= $sizerecabcsi $gcesize"
if [ $sizerecabcsi != $gcesize ]; then
  copylist="$copylist $recabcsi"
fi

if [ "$copylist" != "" ]; then
  echo "====> Copying these files to $copyuri:  $copylist"
  for f in $copylist; do
    echo "Copying $f"
    echo $gsutil cp $f $copyuri
    if [ "$?" != "0" ]; then
      Fail "Failed to copy file to GCE"
    fi
  done
else
  echo "No files to copy, everything for $bamid $nwdid is up to date"
fi

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of files to GCE completed in $etime seconds"
Successful
Log $etime
exit
