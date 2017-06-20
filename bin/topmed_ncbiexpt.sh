#!/bin/bash
#
#   topmed_ncbiexpt.sh -submit bamid
#
#	Send experiment XML for a bamid to NCBI
#
. /usr/cluster/topmed/bin/topmed_actions.inc
topmedxml="/usr/cluster/topmed/bin/topmed_xml.pl"
ascpcmd="$topmedcmd send2ncbi"
me=ncbiexpt
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  #   Do not allow this to play with anything except year one data
  year=`GetDB $1 datayear`
  if [ "$year" != "1" ]; then
    Fail "$0 $* must be year ONE data, not '$year'"
  fi
  MayIRun $me  $bamid
  MyRealHost $bamid 'bam'
  SubmitJob $bamid "topmed-redo" '2G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Send experiment XML to NCBI"
  exit 1
fi
bamid=$1

Started

#   Do not allow this to play with anything except year one data
year=`GetDB $bamid datayear`
if [ "$year" != "1" ]; then
  Fail "$0 $* must be year ONE data, not '$year'"
fi

nwdid=`GetNWDID $bamid`

#   Go to our working directory
cd $console
if [ "$?" != "0" ]; then
  Fail "Unable to CD to '$console' to create XML file"
fi
mkdir XMLfiles 2>/dev/null
cd XMLfiles
here=`pwd`

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -type expt $bamid
if [ "$?" != "0" ]; then
  Fail "Unable to create experiment XML files"
fi

#   Check files we created
files=''
for f in $nwdid-expt.submit.xml $nwdid.expt.xml; do
  if [ ! -f $f ]; then
    Fail "Missing XML file '$f'"   
  fi
  files="$files $f"
done

#   Make a TAR of the XML files to send
stime=`date +%s`
tar cf $nwdid-expt.tar $files
if [ "$?" != "0" ]; then
  Fail "Unable to create TAR of XML files"
fi
files=$nwdid-expt.tar

echo "Sending XML files to NCBI - $files"
$ascpcmd $files
if [ "$?" != "0" ]; then
  Fail "FAILED to send XML files to NCBI - $files"
fi
echo "XML files '$files' sent to NCBI"
etime=`date +%s`
etime=`expr $etime - $stime`
$topmedcmd -persist -emsg "" mark $bamid $markverb delivered    # Cannot use Successful
Log $etime
