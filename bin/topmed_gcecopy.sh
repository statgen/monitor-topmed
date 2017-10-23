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

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me  $bamid
  RandomRealHost $bamid
  SubmitJob $bamid "topmed-$me" '3G' "$0 $*"
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
baserecabcram=`basename $recabcram`
sizerecabcram=`$stat $recabcram`
recabcrai=$recabcram.crai

#   Create the index file as necessary
CreateIndex $bamid $recabcram

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
f=`basename $recabcram`
echo "$b38cramchecksum $f" > $recabcram.md5

#   Pass along flagstat too
flagstat=`GetDB $bamid cramflagstat`
echo "$flagstat + 0 paired in sequencing" > $recabcram.flagstat

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
#   Just to be sure if the cram or bcf is missing, copy the index too
#======================================================================
replacelist=''
rlist=''
cmissinglist=''
cmissing=''
copyuri=`$topmedpath wherepath $nwdid gceupload`

l=(`$gsutil stat $copyuri/$baserecabcram | grep Content-Length:`)
gcesize=${l[1]}
if [ "$sizerecabcram" != "$gcesize" ]; then
  echo "$baserecabcram sizes= $sizerecabcram / $gcesize"
  if [ "$gcesize" = "" ]; then
    cmissing="$cmissing $baserecabcram $baserecabcrai $recabcram.md5 $recabcram.flagstat"
    cmissinglist="$cmissinglist $recabcram $recabcrai $recabcram.md5 $recabcram.flagstat"
  else
    replacelist="$replacelist $recabcram $recabcrai $recabcram.md5 $recabcram.flagstat"
    rlist="$rlist $baserecabcram $baserecabcrai  $recabcram.md5 $recabcram.flagstat"
  fi
else
  l=(`$gsutil stat $copyuri/$baserecabcrai | grep Content-Length:`)
gcesize=${l[1]}
  if [ "$sizerecabcrai" != "$gcesize" ]; then
    echo "$baserecabcrai sizes= $sizerecabcrai / $gcesize"
    if [ "$gcesize" = "" ]; then
      cmissing="$cmissing $baserecabcrai $recabcram.md5 $recabcram.flagstat"
      cmissinglist="$cmissinglist $recabcrai $recabcram.md5 $recabcram.flagstat"
    else
      replacelist="$replacelist $recabcrai $recabcram.md5 $recabcram.flagstat"
      rlist="$rlist $baserecabcrai $recabcram.md5 $recabcram.flagstat"
    fi
  fi
fi

l=(`$gsutil stat $copyuri/$baserecabbcf | grep Content-Length:`)
gcesize=${l[1]}
if [ "$sizerecabbcf" != "$gcesize" ]; then
  echo ""
  echo "$baserecabbcf sizes= $sizerecabbcf / $gcesize"
  if [ "$gcesize" = "" ]; then
    cmissing="$cmissing $baserecabbcf $baserecabcsi"
    cmissinglist="$cmissinglist $recabbcf $recabcsi"
  else
    replacelist="$replacelist $recabbcf $recabcsi"
    rlist="$rlist $baserecabbcf $baserecabcsi"
  fi
else
  l=(`$gsutil stat $copyuri/$baserecabcsi | grep Content-Length:`)
gcesize=${l[1]}
  if [ "$sizerecabcsi" != "$gcesize" ]; then
    echo "$baserecabcsi sizes= $sizerecabcsi / $gcesize"
    if [ "$gcesize" = "" ]; then
      cmissing="$cmissing $baserecabcs"
      cmissinglist="$cmissinglist $recabcsi"
    else
      replacelist="$replacelist $recabcsi"
      rlist="$rlist $baserecabcsi"
    fi
  fi
fi

#   First copy just the files that are missing
if [ "$cmissing" != "" ]; then
  echo ""
  echo ""
  echo "====> Copying MISSING files to $copyuri:  $cmissing"
  for f in $cmissinglist; do
    basef=`basename $f`
    echo "Copying $basef"
    $gsutilbig cp $f $copyuri/$basef
    if [ "$?" != "0" ]; then
      Fail "Failed to copy $basef to GCE"
    fi
  done
fi

#   Now, copy all the files that need replacing
if [ "$replacelist" != "" ]; then
  echo ""
  echo "====> Replace in $copyuri:  $rlist"
#Fail "Replace in $copyuri:  $rlist"
  for f in $replacelist; do
    basef=`basename $f`
    echo "Copying $basef"
    $gsutilbig cp $f $copyuri/$basef
    if [ "$?" != "0" ]; then
      Fail "Failed to copy $basef to GCE"
    fi
  done
fi

echo "GCE files for $bamid $nwdid are up to date"
$gsutil ls -l $copyuri

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of files to GCE completed in $etime seconds"
Successful
Log $etime
exit
