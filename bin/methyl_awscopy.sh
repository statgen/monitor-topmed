#!/bin/bash
#
#   methyl_awscopy.sh -submit| batchid
#
#	Copy files from local storage to AWS for release to the world
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
topmedcmd="$topmedcmd -datatype methyl"
topmedpath="$topmedpath -datatype methyl"
aws="/usr/local/bin/aws --profile nhlbi-data-commons s3"
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
#	Must set bamid, legacy decision from genome
Fail "This script is not ready yet"

bamid=$1
batchid=$1
nwdid=`GetNWDID $sampleid`
sampleid=`GetDB $nwdid sampleid`

Started
stime=`date +%s`

send2aws=`GetDB $sampleid send2aws`
if [ "$send2aws" != "Y" ]; then
  Fail "May not copy $bamid to AWS, send2aws=$send2aws"
fi

awsuri=`$topmedpath wherepath $nwdid awsupload`
if [ "$awsuri" = "" ]; then
  Fail "Unable to figure out AWS URI"
fi

#   Get remapped cram file
recabcram=`$topmedpath wherefile $sampleid b38`
if [ "$recabcram" = "" ]; then
  Fail "Unable to determine CRAM file for '$sampleid'"
fi
if [ ! -f $recabcram ]; then
  Fail "CRAM file for '$sampleid' does not exist: $recabcram"
fi
sizerecabcram=`$cmdstat $recabcram`

#   Create the index file as necessary
CreateIndex $bamid $recabcram
recabcrai=$recabcram.crai
sizerecabcrai=`$cmdstat $recabcrai`

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
