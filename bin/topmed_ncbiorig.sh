#!/bin/bash
#
#   topmed_ncbiorig.sh -submit bamid
#
#	Send the proper set of files to NCBI for the original/secondary files
#
. /usr/cluster/topmed/bin/topmed_actions.inc
topmedxml="/usr/cluster/topmed/bin/topmed_xml.pl"
ascpcmd="$topmedcmd send2ncbi"
me=ncbiorig
markverb=$me
version=secondary

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
  echo "Send original files to NCBI"
  exit 1
fi
bamid=$1
bamid=`$topmedcmd show $1 bamid`

Started

#   Do not allow this to play with anything except year one data
year=`GetDB $bamid datayear`
if [ "$year" != "1" ]; then
  Fail "$0 $* must be year ONE data, not '$year'"
fi

nwdid=`GetNWDID $bamid`

center=`GetDB $bamid center`
if [ "$center" = "" ]; then
  Fail "Invalid bamid '$bamid'. CENTER not known"
fi
build=`GetDB $bamid build`
if [ "$build" = "" ]; then
  Fail "Invalid build '$build'. BUILD not known"
fi
origbam=`GetDB $bamid bamname`
if [ "$origbam" = "" ]; then
  Fail "Invalid bamid '$bamid'. BAMNAME not known"
fi

#   Go to our working directory
cd $console
if [ "$?" != "0" ]; then
  Fail "Unable to CD to '$console' to create XML file"
fi
mkdir XMLfiles 2>/dev/null
cd XMLfiles
here=`pwd`

#   Figure out what file to send. Either a BAM or a CRAM
sendcram='N'
datayear=`GetDB $bamid datayear`  # What year data is this
if [ "$center" = "broad" -a "$datayear" != '1' ]; then
  sendcram='Y'
fi

if [ "$sendcram" = "Y" ]; then
  l=(`$topmedpath wherefile $bamid cram`)        # Get backupdir and backupfile and host
  sendfile="${l[0]}"
  sf=$nwdid.src.cram
  ln -sf $sendfile $sf
  checksum=`$topmedcmd -persist show $bamid cramchecksum`
  sendfile=$sf
  echo "Sending backup cram instead of original file ($sendfile)"
  # It's not supposed to happen, but sometimes the CRAM checksum is missing
  if [ "$checksum" = "" ]; then
    echo "Calculating MD5 for CRAM"
    stime=`date +%s`
    checksum=`calcmd5 $sendfile | awk '{print $1}'`
    if [ "$checksum" != "" ]; then
      $topmedcmd set $bamid cramchecksum $checksum
      etime=`date +%s`
      etime=`expr $etime - $stime`
      echo "MD5 calculated for CRAM in $etime seconds"
    fi
  fi
else
  sendfile=`$topmedpath wherefile $bamid bam`
  if [ ! -f "$sendfile" ]; then
      sendcram=Y                        # No BAM, send CRAM
      echo "Attempting to send the CRAM, rather than the missing BAM"
      sendfile=`$topmedpath wherefile $bamid cram`
      ln -sf $sendfile $nwdid.src.cram
      checksum=`$topmedcmd -persist show $bamid cramchecksum`
      sendfile=$nwdid.src.cram
  else
    ln -sf $sendfile $nwdid.src.bam     # BAM file exists, send it
    sendfile=$nwdid.src.bam
    checksum=`$topmedcmd -persist show $bamid checksum`
  fi
fi

#   NCBI cannot process crams anymore, so we convert this CRAM to a BAM
rmbam=''
ext=${sendfile##*.}
if [ "$ext" = "cram" ]; then
  newbam=/tmp/$nwdid.src.bam
  if [ ! -f $newbam ]; then
    echo "Converting CRAM $sendfile to $newbam"
    $samtools view -b -@ 2 $sendfile > $newbam
    if [ "$?" != "0" ]; then
      Fail "Unable to convert $sendfile to $newbam"
    fi
  else
    echo "Using existing file $newbam"
  fi
  echo "Calculating checksum for converted CRAM"
  md5=(`md5sum $newbam`)
  checksum=${md5[0]}
  if [ "$checksum" = "" ]; then
    Fail "Unable to calculate MD5 for $newbam"
  fi
  echo "Checksum for converted CRAM is $checksum"
  sendfile=`basename $newbam`
  ln -sf $newbam $sendfile 
  ls -l $sendfile $newbam
  rmbam=$newbam              # Remove this file
fi

if [ "$checksum" = "" ]; then
  Fail "Invalid bamid '$bamid' ($sendfile). CHECKSUM not known"
fi

if [ ! -f "$sendfile" ]; then
  Fail "Original file for bamid '$bamid' ($sendfile) not found"
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid file=$sendfile ========="

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -build $build -type $version $bamid $sendfile $checksum
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

echo "Sending data file to NCBI - $sendfile"
ls -l $sendfile
stime=`date +%s`
$ascpcmd $sendfile
rc=$?
rm -f $sendfile
if [ "$rmbam" != '' ]; then
  rm -f $rmbam
fi
if [ "$rc" != "0" ]; then
  Fail "FAILED to send data file '$sendfile' (rc=$rc)"
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "File sent in $etime seconds"

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
