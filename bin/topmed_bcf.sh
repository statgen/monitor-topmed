#!/bin/bash
#
#   topmed_bcf.sh -submit| bamid
#
#	Create BCF file for a remapped CRAM
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
me=bcf
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  cramfile=`$topmedpath wherefile $bamid b$build`
  if [ "$cramfile" = "" ]; then
    Fail "Failed to find '$cramfile'"
  fi
  filesize=`stat --printf=%s $cramfile`
  if [ "$filesize" = "" ]; then
    Fail "Unable to get filesize for '$cramfile'"
  fi
  timeout='12:00:00'
  mem=32G
  #   Large crams can require LOTS of memory and can take a long time
  if [ "$filesize" -gt "40255183256" ]; then
    mem=132G
    timeout='28:00:00'
  elif [ "$filesize" -gt "35044692321" ]; then
    mem=96G
    timeout='26:00:00'
  elif [ "$filesize" -gt "31848365056" ]; then
    mem=84G
    timeout='24:00:00'
  elif [ "$filesize" -gt "23023087355" ]; then
    mem=69G
	timeout='22:00:00'
  elif [ "$filesize" -gt "1566897485" ]; then
    mem=40G
	timeout='20:00:00'
  fi
  if [ "$FORCEMEMORY" != "" ]; then
    mem=$FORCEMEMORY
  fi
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  SubmitJob $bamid "$PROJECT-bcf" $mem "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Create BCF file for a remapped CRAM"
  exit 1
fi
bamid=$1
nwdid=`GetNWDID $bamid`
bamid=`GetDB $nwdid bamid`

Started
stime=`date +%s`

#======================================================================
#   Post process the remapped CRAM from Google Cloud
#======================================================================
#   CD to remapped cram file location
crampath=`$topmedpath wherepath $bamid b$build`
if [ "$crampath" = "" ]; then
  Fail  "Unable to determine where remapped CRAM file should be for '$bamid'"
fi
cd $crampath
if [ "$?" != "0" ]; then
  Fail "Unable to CD to directory for remapped CRAM file '$bamid'"
fi

cramfile=`$topmedpath wherefile $bamid b$build`
if [ ! -f $cramfile ]; then
  Fail "Unable to find remapped cram file '$cramfile'"
fi

#   Create crai for cram if it does not exist
CreateIndex $bamid $cramfile

#   Create the BCF file
bcfdir=`$topmedpath wherepath $bamid bcf`
mkdir -p $bcfdir
if [ "$?" != "0" ]; then
  Fail "Unable to create BCF output directory for '$bamid' [$nwdid] - $bcfdir"
fi
bcffile=$bcfdir/$nwdid.bcf

set -o pipefail
$samtools view -uh $cramfile | $bam clipOverlap --poolSize 100000000 --in -.ubam --out -.ubam | $vt discover2 -z -q 20 -b + -r $vtref -s $nwdid -o $bcffile
if [ "$?" != "0" ]; then
  Fail "Unable to run VT DISCOVER for bamid '$bamid' [$nwdid] on $cramfile creating $bcffile"
fi
$bcftools index $bcffile
if [ "$?" != "0" ]; then
  Fail "Unable to run BCFTOOLS for bamid '$bamid' [$nwdid] on $bcffile"
fi

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Created BCF file for remapped CRAM ($crampath) completed in $etime seconds"
SetDB $bamid 'state_gce38bcf' 20
if [ "$PROJECT" = "topmed" ]; then
  SetDB $bamid 'state_gce38cpbcf' 1
fi

Successful
Log $etime
exit
