#!/bin/bash
#
#   topmed_bcf.sh -submit| bamid
#
#	Create BCF file for a remapped CRAM
#
. /usr/cluster/topmed/bin/topmed_actions.inc
vt=/net/wonderland/home/lefaivej/bin/vt
vtref=/net/topmed/working/mapping/gotcloud/ref/hg38/hs38DH.fa

me=bcf
mem=2G
markverb="${me}ed"
constraint="--constraint eth-10g"
qos="--qos=topmed-$me"
slurmp=topmed-working
realhost=''
build=38

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test $me $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  # Run this on node where remapped cram lives
  #h=`$topmedpath -fallback whathost $1 b$build`
  #if [ "$h" != "" ]; then
  #  realhost="--nodelist=$h"
  #  #qos="--qos=$h-$me"
  #fi

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` $me ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
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

d=`date +%Y/%m/%d`
s=`hostname`
p=`pwd`

crampath=`$topmedpath -fallback wherepath $bamid b$build`
if [ "$crampath" = "" ]; then
  echo "Unable to determine where remapped CRAM file should be for '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid crampath=$crampath pwd=$p ========="

#======================================================================
#   Post process the remapped CRAM from Google Cloud
#======================================================================
#   Mark this as started
$topmedcmd mark $bamid $markverb started

#   CD to remapped cram file location
cd $crampath
if [ "$?" != "0" ]; then
  echo "Unable to CD to directory for remapped CRAM file '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Unable to get expt_sampleid for bamid '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
cramfile=`$topmedpath -fallback wherefile $bamid b$build`
if [ ! -f $cramfile ]; then
  echo "Unable to find remapped cram file '$cramfile'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

#   Create crai for cram if it does not exist
crai="$cramfile.crai"
if [ ! -f $crai ]; then
  echo "Creating index file '$cramfile'"
  $samtools index $cramfile 2>&1
  if [ "$?" != "0" ]; then
    echo "Unable to create index file for bamid '$bamid' [$nwdid]"
    $topmedcmd -persist mark $bamid $markverb failed
    exit 2
  fi
  echo "Created CRAI file for bamid '$bamid' [$nwdid]"
fi

#   Create the BCF file
bcfdir=`$topmedpath wherepath $bamid bcf`
mkdir -p $bcfdir
if [ "$?" != "0" ]; then
  echo "Unable to create BCF output directory for '$bamid' [$nwdid] - $bcfdir"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
bcffile=$bcfdir/$nwdid.bcf

$vt discover -b $cramfile -s $nwdid -r $vtref -o $bcffile
if [ "$?" != "0" ]; then
  echo "Unable to run VT DISCOVER for bamid '$bamid' [$nwdid] on $cramfile creating $bcffile"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
bcftools index $bcffile
if [ "$?" != "0" ]; then
  echo "Unable to run BCFTOOLS for bamid '$bamid' [$nwdid] on $bcffile"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Created BCF file for remapped CRAM ($crampath) completed in $etime seconds"
$topmedcmd -persist mark $bamid $markverb completed
echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit
