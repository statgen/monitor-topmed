#!/bin/bash
#
#   topmed_awscopy.sh -submit| bamid
#
#	Copy files from local storage to AWS for release to the world
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
aws="/usr/local/bin/aws --profile nhlbi-data-commons s3"
partpath=b38.irc.v1
me=awscopy
markverb=$me
cmdstat="stat --printf=%s --dereference"

#------------------------------------------------------------------
# Subroutine:
#   Copy2AWS(f, awsfile)
#   Copy file to AWS as awsfile
#------------------------------------------------------------------
function Copy2AWS {
  local f=$1
  local awsfile=$2
  echo "====> Copying $f"
  $aws cp --quiet $f $awsuri/$awsfile
  if [ "$?" != "0" ]; then
    Fail "Failed to copy $f to AWS as $awsfile"
  fi
}

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='2:00:00'
  SubmitJob $bamid "$PROJECT-aws" '2G' "$0 $*"
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
nwdid=`GetNWDID $bamid`
bamid=`GetDB $nwdid bamid`

Started
stime=`date +%s`

send2aws=`GetDB $bamid send2aws`
if [ "$send2aws" != "Y" ]; then
  Fail "May not copy $bamid to AWS, send2aws=$send2aws"
fi

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
sizerecabcram=`$cmdstat $recabcram`

#   Create the index file as necessary
CreateIndex $bamid $recabcram
recabcrai=$recabcram.crai
sizerecabcrai=`$cmdstat $recabcrai`

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
#   All the files of interest exist locally, copy to AWS if changed
#======================================================================
#   Get size of files in AWS
tmpf=/run/shm/$$
$aws ls $awsuri/$nwdid > $tmpf
l=(`grep -e $nwdid.$partpath.cram\$ $tmpf`)
awscramsize=${l[2]}
l=(`grep $nwdid.$partpath.cram.crai $tmpf`)
awscraisize=${l[2]}
rm -f $tmpf

#   Send files that have changed size
if [ "$sizerecabcram" != "$awscramsize" ]; then
  Copy2AWS $recabcram $nwdid.$partpath.cram
  Copy2AWS $recabcram.crai $nwdid.$partpath.cram.crai
else
  echo "No need to send unchanged $nwdid.$partpath.cram"
  if [ "$sizerecabcrai" != "$awscraisize" ]; then
    echo Copy2AWS $recabcram.crai $nwdid.$partpath.cram.crai "$sizerecabcrai" != "$awscraisize" 
  else
    echo "No need to send unchanged $nwdid.$partpath.cram.crai"
  fi
fi
$aws ls $awsuri/$nwdid

etime=`date +%s`
etime=`expr $etime - $stime`

echo "Copy of files to AWS completed in $etime seconds"
Successful
Log $etime
exit
