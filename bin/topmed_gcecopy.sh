#!/bin/bash
#
#   topmed_gcecopy.sh -submit| bamid
#
#	Copy files from local storage to a GCE bucket for GCE processing
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcecopy
markverb=$me

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

#   Create metadata, same place where other remapped data is
p=`$topmedpath wherepath $bamid b38`
metafile=$p/$nwdid.metadata.txt
crammd5=`GetDB $bamid b${build}cramchecksum`
craimd5=`GetDB $bamid b${build}craichecksum`
flagstat=`GetDB $bamid b${build}flagstat`
studyname=`GetDB $bamid studyname`

echo "$nwdid.b$build.irc.cram-md5 $crammd5 $nwdid.b$build.irc.crai-md5 $craimd5 $nwdid.b$build.irc.nonredundant_sequence_reads $flagstat studyname $studyname" > $metafile

#======================================================================
#   All the files of interest exist locally
#   Copy is very fast to GCE, so don't try to be too clever
#   When a file does not change name, copy as a group
#======================================================================
copyuri=`$topmedpath wherepath $nwdid gceupload`
Copy2GCE $recabcram `basename $recabcram`
$gsutil cp $recabcrai $recabbcf $recabcsi $metafile $copyuri 
if [ "$?" != "0" ]; then
  Fail "Failed to copy files to GCE as $gcefile"
fi

#   For a few months, clean out files we no longer want
$gsutil rm $copyuri/\*.md5 $copyuri/\*.flagstat
rm -f $p/*.md5 $p/*.flagstat

echo "GCE files for $bamid $nwdid are up to date"
$gsutil ls -l $copyuri

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of files to GCE completed in $etime seconds"
Successful
Log $etime
exit
