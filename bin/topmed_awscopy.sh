#!/bin/bash
#
#   topmed_awscopy.sh -submit| bamid
#
#	Copy files from local storage to AWS for release to the world
#

. /usr/cluster/topmed/bin/topmed_actions.inc
aws="/usr/local/bin/aws --profile nhlbi-data-commons s3"
me=awscopy
markverb=$me
stat="stat --printf=%s --dereference"

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me  $bamid
  RandomRealHost $bamid
  SubmitJob $bamid "topmed" '2G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Copy files from local storage to AWS for release to the world"
  exit 1
fi
bamid=$1

Started
nwdid=`GetNWDID $bamid`
stime=`date +%s`

awsuri=`$topmedpath wherepath $nwdid awsupload`
if [ "$awsuri" = "" ]; then
  Fail "Unable to figure out AWS URI"
fi

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
baserecabcrai=`basename $recabcrai`
sizerecabcrai=`$stat $recabcrai`

#======================================================================
#   All the files of interest exist locally
#   Copy the file to the remote destination only if the size differs
#   If cram is sent, copy the index too
#======================================================================
replacelist=''
rlist=''
cmissinglist=''
cmissing=''

studyname=`GetDB $bamid studyname`
piname=`GetDB $bamid piname`
datayear=`GetDB $bamid datayear`
phs=`GetDB $bamid phs`
metadata="--metadata sample-id=$nwdid,year=$datayear,study=$studyname,pi=$piname,phs=$phs" 

lcstudyname=${studyname,,}              # Force to lower case cause AWS wants it that way
copyuri=$awsuri/$lcstudyname
l=(`$aws ls $copyuri/$nwdid | grep cram`)
awscramsize=''
if [ "${l[3]}" = "$nwdid.cram" ]; then  # Found cram in AWS
  awscramsize=${l[2]}                   # Size of cram
fi
awscraisize=''
if [ "${l[7]}" = "$nwdid.cram.crai" ]; then  # Found crai in AWS
  awscraisize=${l[6]}                   # Size of crai
fi

#sizerecabcrai=`expr $sizerecabcrai - 1`     # Force a small file to resend

if [ "$sizerecabcram" != "$awscramsize" ]; then
  echo "$baserecabcram sizes= $sizerecabcram / $awscramsize"
  if [ "$awscramsize" = "" ]; then
    cmissing="$cmissing $baserecabcram $baserecabcrai"
    cmissinglist="$cmissinglist $recabcram $recabcrai"
  else
    replacelist="$replacelist $recabcram $recabcrai"
    rlist="$rlist $baserecabcram $baserecabcrai"
  fi
else
  if [ "$sizerecabcrai" != "$awscraisize" ]; then
    echo "$baserecabcrai sizes= $sizerecabcrai / $awscraisize"
    if [ "$awscraisize" = "" ]; then
      cmissing="$cmissing $baserecabcrai"
      cmissinglist="$cmissinglist $recabcrai"
    else
      replacelist="$replacelist $recabcrai"
      rlist="$rlist $baserecabcrai"
    fi
  fi
fi

#   First copy just the files that are missing
if [ "$cmissing" != "" ]; then
  echo ""
  echo "====> Copying MISSING files to $copyuri:  $cmissing"
  for f in $cmissinglist; do
    basef=`basename $f | sed -e 's/.recab//'`    # Remove recab from our names
    echo "Creating missing $basef"
    $aws cp --quiet $metadata $f $copyuri/$basef
    if [ "$?" != "0" ]; then
      Fail "Failed to copy $basef to AWS"
    fi
  done
fi

#   Now, copy all the files that need replacing
if [ "$replacelist" != "" ]; then
  echo ""
  echo "====> Replace in $copyuri:  $rlist"
  for f in $replacelist; do
    basef=`basename $f | sed -e 's/.recab//'`    # Remove recab from our names
    echo "Replacing $basef"
    $aws cp --quiet $metadata $f $copyuri/$basef
    if [ "$?" != "0" ]; then
      Fail "Failed to copy $basef to AWS"
    fi
  done
fi

echo "AWS files for $bamid $nwdid are up to date"
$aws ls $copyuri/$nwdid

#echo "Here is metadata for the cram"
#bucket=`$topmedpath wherepath $nwdid awsbucket`
#bucketpath=`$topmedpath wherepath $nwdid awsbucketpath`
#${aws}api head-object --bucket $bucket --key $bucketpath/$lcstudyname/$nwdid.cram.crai

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of files to AWS completed in $etime seconds"
Successful
Log $etime
exit
