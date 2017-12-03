#!/bin/bash
#
#   topmed_gcecopy.sh -submit| bamid
#
#	Copy files from local storage to a GCE bucket for GCE processing
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcecopy
markverb=$me
stat="stat --printf=%s --dereference"

#------------------------------------------------------------------
# Subroutine:
#   Copy2GCE(f, gcefile)
#   Copy file to GCE as gcefile
#------------------------------------------------------------------
function Copy2GCE {
  local f=$1
  local gcefile=$2
  if [ "${f##*.}" = "cram" ]; then
    gs=$gsutilbig
  else
    gs=$gsutil
  fi
  echo "====> Copying to $gcefile"
  $gs cp $f $copyuri/$gcefile       # Someday we will upgrade and -q will work
  if [ "$?" != "0" ]; then
    Fail "Failed to copy $f to GCE as $gcefile"
  fi
}

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
stime=`date +%s`

#   Get remapped cram file
recabcram=`$topmedpath wherefile $bamid b38`
if [ "$recabcram" = "" ]; then
  Fail "Unable to determine CRAM file for '$bamid'"
fi
if [ ! -f $recabcram ]; then
  Fail "CRAM file for '$bamid' does not exist: $recabcram"
fi
sizerecabcram=`$stat $recabcram`

#   Create the index file as necessary
CreateIndex $bamid $recabcram
recabcrai=$recabcram.crai
sizerecabcrai=`$stat $recabcrai`

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

#   Get bcf file for remapped cram file
#recabbcf=`$topmedpath wherefile $bamid bcf`
#if [ "$recabbcf" = "" ]; then
#  Fail "Unable to determine BCF file for '$bamid'"
#fi
#if [ ! -f $recabbcf ]; then
#  Fail "BCF file for '$bamid' does not exist: $recabbcf"
#fi
#sizerecabbcf=`$stat $recabbcf`

#recabcsi=$recabbcf.csi
#if [ ! -f $recabcsi ]; then
#  echo "CSI missing, trying to create it"
#  $bcftools index $recabbcf
#  if [ "$?" != "0" ]; then
#    Fail "Unable to create CSI for $bamid' for: $recabbcf"
#  fi
#fi
#if [ ! -f $recabcsi ]; then
#  Fail "CSI file for '$bamid' does not exist: $recabcsi"
#fi
#sizerecabcsi=`$stat $recabcsi`

#======================================================================
#   All the files of interest exist locally
#   Copy the file to the remote destination only if the size differs
#======================================================================
copyuri=`$topmedpath wherepath $nwdid gceupload`

#   Now get the sizes of files of interest in AWS
baserecabcram=`basename $recabcram`
l=(`$gsutil stat $copyuri/$baserecabcram | grep Content-Length:`)
gcesizerecabcram=${l[1]}

baserecabcrai=$baserecabcram.crai
l=(`$gsutil stat $copyuri/$baserecabcrai | grep Content-Length:`)
gcesizerecabcrai=${l[1]}

#baserecabbcf=`basename $recabbcf`
#l=(`$gsutil stat $copyuri/$baserecabbcf | grep Content-Length:`)
#gcesizerecabbcf=${l[1]}

#baserecabcsi=$baserecabbcf.csi
#l=(`$gsutil stat $copyuri/$baserecabcsi | grep Content-Length:`)
#gcesizerecabcsi=${l[1]}

#   Force everthing to be sent, usually for testing
#gcesizerecabcram=`expr $gcesizerecabcram - 1`
#gcesizerecabcrai=`expr $gcesizerecabcrai - 1`
#gcesizerecabbcf=`expr $gcesizerecabbcf - 1`
#gcesizerecabcsi=`expr $gcesizerecabcsi - 1`

#   Send files that have changed size
if [ "$sizerecabcram" != "$gcesizerecabcram" ]; then
  tmpf=/run/shm/$$
  Copy2GCE $recabcram $baserecabcram
  echo "$b38cramchecksum $baserecabcram" > $tmpf
  Copy2GCE $tmpf $baserecabcram.md5
  flagstat=`GetDB $bamid b38flagstat`
  echo "$flagstat + 0 paired in sequencing" > $tmpf
  Copy2GCE $tmpf $baserecabcram.flagstat
  rm -f $tmpf
else
  echo "No need to send unchanged $baserecabcram"
fi
if [ "$sizerecabcrai" != "$gcesizerecabcrai" ]; then
  Copy2GCE $recabcrai $baserecabcrai
else
  echo "No need to send unchanged $baserecabcrai"
fi
#if [ "$sizerecabbcf" != "$gcesizerecabbcf" ]; then
#  Copy2GCE $recabbcf $baserecabbcf
#else
#  echo "No need to send unchanged $baserecabbcf"
#fi
#if [ "$sizerecabcsi" != "$gcesizerecabcsi" ]; then
#  Copy2GCE $recabcsi $baserecabcsi
#else
#  echo "No need to send unchanged $baserecabcsi"
#fi

echo "GCE files for $bamid $nwdid are up to date"
$gsutil ls -l $copyuri

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of files to GCE completed in $etime seconds"
Successful
Log $etime
exit
