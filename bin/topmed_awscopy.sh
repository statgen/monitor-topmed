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

#------------------------------------------------------------------
# Subroutine:
#   Copy2AWS(f, awsfile)
#   Copy file to AWS as awsfile
#------------------------------------------------------------------
function Copy2AWS {
  local f=$1
  local awsfile=$2
  echo "====> Copying $f"
  $aws cp --quiet $metadata $f $copyuri/$awsfile
  if [ "$?" != "0" ]; then
    Fail "Failed to copy $f to AWS as $awsfile"
  fi
}

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
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

#   Get bcf file
#recabvcf=`$topmedpath wherefile $bamid bcf`
#if [ "$recabvcf" = "" ]; then
#  Fail "Unable to determine BCF file for '$bamid'"
#fi
#if [ ! -f $recabvcf ]; then
#  Fail "BCF file for '$bamid' does not exist: $recabvcf"
#fi
#sizerecabvcf=`$stat $recabvcf`

#======================================================================
#   All the files of interest exist locally
#   Copy the file to the remote destination only if the size differs
#======================================================================
studyname=`GetDB $bamid studyname`
piname=`GetDB $bamid piname`
datayear=`GetDB $bamid datayear`
phs=`GetDB $bamid phs`
metadata="--metadata sample-id=$nwdid,year=$datayear,study=$studyname,pi=$piname,phs=$phs" 

lcstudyname=${studyname,,}              # Force to lower case cause AWS wants it that way
copyuri=$awsuri/$lcstudyname

#   Now get the sizes of files of interest in AWS
tmpf=/run/shm/$$
$aws ls $copyuri/$nwdid > $tmpf
l=(`grep -e $nwdid.cram\$ $tmpf`)
awscramsize=${l[2]}
l=(`grep $nwdid.cram.crai $tmpf`)
awscraisize=${l[2]}
#l=(`grep $nwdid.cram.vcf $tmpf`)
#awsvcfsize=${l[2]}
rm -f $tmpf

#   Send files that have changed size
if [ "$sizerecabcram" != "$awscramsize" ]; then
  Copy2AWS $recabcram $nwdid.cram
else
  echo "No need to send unchanged $nwdid.cram"
fi
if [ "$sizerecabcrai" != "$awscraisize" ]; then
  Copy2AWS $recabcram.crai $nwdid.cram.crai
else
  echo "No need to send unchanged $nwdid.cram.crai"
fi
#if [ "$sizerecabvcf" != "$awsvcfsize" ]; then
#  Copy2AWS $recabvcf $nwdid.cram.vcf
#else
#  echo "No need to send unchanged $nwdid.cram.vcf"
#fi

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
