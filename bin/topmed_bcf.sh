#!/bin/bash
#
#   topmed_bcf.sh -submit| bamid
#
#	Create BCF file for a remapped CRAM
#
. /usr/cluster/topmed/bin/topmed_actions.inc
vt=/usr/cluster/software/trusty/topmed-year1-freeze3a/master/vt/vt
vtref=/net/topmed/working/mapping/gotcloud/ref/hg38/hs38DH.fa
bcftools=/usr/cluster/bin/bcftools
bam=/usr/cluster/bin/bam

me=bcf
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me $bamid
  MyRealHost $bamid b$build
  SubmitJob $bamid "topmed-$me" '24G' "$0 $*"
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

Started
GetNWDID $bamid
stime=`date +%s`
crampath=`$topmedpath wherepath $bamid b$build`
if [ "$crampath" = "" ]; then
  Fail  "Unable to determine where remapped CRAM file should be for '$bamid'"
fi

#======================================================================
#   Post process the remapped CRAM from Google Cloud
#======================================================================
#   CD to remapped cram file location
cd $crampath
if [ "$?" != "0" ]; then
  $topmedcmd -persist emsg "Unable to CD to directory for remapped CRAM file '$bamid'"
  exit 2
fi

cramfile=`$topmedpath wherefile $bamid b$build`
if [ ! -f $cramfile ]; then
  Fail "Unable to find remapped cram file '$cramfile'"
fi

#   Create crai for cram if it does not exist
crai="$cramfile.crai"
if [ ! -f $crai ]; then
  echo "Creating index file '$cramfile'"
  $samtools index $cramfile 2>&1
  if [ "$?" != "0" ]; then
    Fail "Unable to create index file for bamid '$bamid' [$nwdid]"
  fi
  echo "Created CRAI file for bamid '$bamid' [$nwdid]"
fi

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
SetDB $bamid 'state_gce38bcf_push' 20
SetDB $bamid 'state_gce38bcf_pull' 20
#markverb=$me
SetDB $bamid 'state_gce38bcf' 20

Successful
Log $etime
exit
