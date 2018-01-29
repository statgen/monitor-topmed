#!/bin/bash
#
#   topmed_gcecopy.sh -submit| bamid
#
#	Copy files from local storage to a GCE SHARE bucket for GCE processing
#   This is not the same as copying to gs://topmed-bcf or other of OUR buckets
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcecopy
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  RandomRealHost $bamid
  #MyRealHost $bamid b$build
  MayIRun $me $bamid $realhost
  timeout='2:00:00'
  SubmitJob $bamid "topmed" '3G' "$0 $*"
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
nwdid=`GetNWDID $bamid`
bamid=`GetNWDID $bamid`
stime=`date +%s`

#   Get remapped cram file
recabcram=`$topmedpath wherefile $bamid b38`
if [ "$recabcram" = "" ]; then
  Fail "Unable to determine CRAM file for '$bamid'"
fi
if [ ! -f $recabcram ]; then
  Fail "CRAM file for '$bamid' does not exist: $recabcram"
fi

#   Create the index file as necessary
CreateIndex $bamid $recabcram
recabcrai=$recabcram.crai

#   Get checksum for b38 file or calculate it
b38cramchecksum=`GetDB $bamid b38cramchecksum`
if [ "$b38cramchecksum" = "" ]; then
  echo "Calculating MD5 for $recabcram"
  stime=`date +%s`
  b38cramchecksum=`CalcMD5 $bamid $recabcram`
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "MD5SUM for recab completed in $etime seconds"
  SetDB $bamid b38cramchecksum $b38cramchecksum
fi

#======================================================================
#   All the files of interest exist locally
#   Copy is very fast to GCE, so don't try to be too clever
#======================================================================
copyuri=`$topmedpath wherepath $nwdid gceupload`
gcecramname=$nwdid.b${build}.irc.v1.cram

$gsutil $gsshareoption cp $recabcram $copyuri/$gcecramname
if [ "$?" != "0" ]; then
  Fail "Failed to copy $recabcram to GCE as $gcecramname"
fi
echo "Copied CRAM to $copyuri"
gsutil $gsshareoption cp $recabcrai $copyuri/$gcecramname.crai
if [ "$?" != "0" ]; then
  Fail "Failed to copy $recabcrai to GCE as $gcecramname.crai"
fi
echo "Copied CRAI to $copyuri"

echo "GCE SHARE files for $bamid $nwdid"
gsutil $gsshareoption ls -l $copyuri/${nwdid}\*

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of CRAM files to GCE completed in $etime seconds"
Successful
Log $etime
exit
