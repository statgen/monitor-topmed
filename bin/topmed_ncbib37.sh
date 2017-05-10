#!/bin/bash
#
#   topmed_ncbib37.sh -xmlonly -submit bamid
#
#	Send the proper set of files to NCBI for the remapped bam build 37
#
. /usr/cluster/topmed/bin/topmed_actions.inc
topmedxml="/usr/cluster/topmed/bin/topmed_xml.pl"
ascpcmd="$topmedcmd send2ncbi"
me=ncbib37
markverb=$me
build=37
version=secondary

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  #   Do not allow this to play with anything except year one data
  year=`$topmedcmd show $1 datayear`
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
  echo "Send Remapped Build $build CRAM file to NCBI"
  exit 1
fi
bamid=$1

Started

#   Do not allow this to play with anything except year one data
year=`$topmedcmd show $bamid datayear`
if [ "$year" != "1" ]; then
  Fail "$0 $* must be year ONE data, not '$year'"
fi

GetNWDID $bamid

center=`$topmedcmd show $bamid center`
if [ "$center" = "" ]; then
  Fail "Invalid bamid '$bamid'. CENTER not known"
fi
piname=`$topmedcmd show $bamid piname`
if [ "$piname" = "" ]; then
  Fail "Invalid bamid '$bamid'. PINAME not known"
fi

#   Go to our working directory
cd $console
if [ "$?" != "0" ]; then
  Fail "Unable to CD to '$console' to create XML file"
fi
mkdir XMLfiles 2>/dev/null
cd XMLfiles
here=`pwd`

#   Figure out what file to send
$topmedpath wherefile $bamid b$build
if [ "$?" = '1' ]; then
  Fail "Multiple b$build crams were found - cannot send to NCBI"
fi

sendfile=`$topmedpath wherefile $bamid b$build`
if [ "$sendfile" = '' ]; then
  Fail "Remapped CRAM for Build $build for bamid '$bamid' not found ($sendfile)"
fi
checksum=`$topmedcmd show $bamid b${build}cramchecksum`
if [ "$checksum" = "" ]; then
  echo "Calculating MD5 for CRAM"
  stime=`date +%s`
  checksum=`md5sum $sendfile | awk '{print $1}'`
  if [ "$checksum" = "" ]; then
    Fail "Unable to calculate the MD5 for CRAM: md5sum $sendfile"
  fi
  $topmedcmd set $bamid b${build}cramchecksum $checksum
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "MD5 calculated for CRAM in $etime seconds"
else
  echo "Obtained b${build}cramchecksum for bamid $bamid from database ($checksum)"
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid file=$sendfile ========="

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -type $version -build $build $bamid $sendfile $checksum
if [ "$?" != "0" ]; then
  Fail "Unable to create $version run XML files"
fi

#   Check files we created
files=''
for f in $nwdid-$version.$build.submit.xml $nwdid-$version.$build.run.xml; do
  if [ ! -f $f ]; then
    Fail "Missing XML file '$f'"   
  fi
  files="$files $f"
done

#   Make a TAR of the XML files to send
tar cf $nwdid-$version.$build.tar $files
if [ "$?" != "0" ]; then
  Fail "Unable to create TAR of XML files"
fi
files=$nwdid-$version.$build.tar

if [ "$xmlonly" != "Y" ]; then
  echo "Sending data file to NCBI - $sendfile"
  ls -l $sendfile

  stime=`date +%s`
  $ascpcmd $sendfile
  rc=$?
  #rm -f $sendfile
  if [ "$rc" != "0" ]; then
    Fail "FAILED to send data file '$sendfile'"
  fi
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "File sent in $etime seconds"
else
  echo "XMLONLY specified, did not send CRAM file"
fi

echo "Sending XML files to NCBI - $files"
$ascpcmd $files
if [ "$?" != "0" ]; then
  Fail "FAILED to send data XML files to NCBI - $files"
fi

echo "XML files '$files' sent to NCBI"
etime=`date +%s`
etime=`expr $etime - $stime`
$topmedcmd -persist -emsg "" mark $bamid $markverb delivered    # Cannot use Successful
Log $etime
