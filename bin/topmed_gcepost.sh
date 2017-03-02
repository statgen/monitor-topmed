#!/bin/bash
#
#   topmed_gcepost.sh -submit| bamid
#
#	Post-process a remapped CRAM after being pulled from Google Cloud
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcepost
mem=2G
markverb="${me}ed"
constraint="--constraint eth-10g"
qos="--qos=topmed-$me"
slurmp=topmed-working
realhost=''
gsutil='gsutil -o GSUtil:parallel_composite_upload_threshold=150M'
incominguri='gs://topmed-recabs'
build=38

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedpermit permit test $me $1
  if [ "$?" = "0" ]; then
    echo "$me $1 not permitted" | tee $console/$1-$me.out
    exit 4
  fi 

  # Run this on node where remapped cram lives
  #h=`$topmedpath whathost $1 b$build`
  #if [ "$h" != "" ]; then
  #  realhost="--nodelist=$h"
  #  #qos="--qos=$h-$me"
  #fi

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*`)
  if [ "$?" != "0" ]; then
    $topmedcmd mark $1 $markverb failed
    echo "Failed to submit command to SLURM - $l" > $console/$1-$me.out
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*" >> $console/$1-$me.out
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
  echo "Post-process a remapped CRAM after being pulled from Google Cloud"
  exit 1
fi
bamid=$1

d=`date +%Y/%m/%d`
s=`hostname`
p=`pwd`

crampath=`$topmedpath wherepath $bamid b$build`
if [ "$crampath" = "" ]; then
  echo "Unable to determine where remapped CRAM file should go for '$bamid'"
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

#   Start by looking for file in GCE.  If there, get the MD5 for it and compare with ours
stime=`date +%s`
cramfile=$nwdid.recab.cram
if [ ! -f $cramfile ]; then
  echo "Unable to find the remapped data for '$bamid' [$nwdid] cramfile=$cramfile"
  pwd
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

#   Compare MD5 of file in GCE and local file. If this is rerun, the GCE file might not exist
echo "Calculating MD5 for local file ($cramfile)"
md5=(`md5sum $cramfile`)
md5=${md5[0]}
if [ "$md5" = "" ]; then
  echo "Unable to calculate MD5 for remapped '$bamid' [$nwdid] cramfile=$cramfile"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
$topmedcmd -persist set $bamid b38cramchecksum $md5
echo "Set checksum for b$build file"

export BOTO_CONFIG=/net/topmed/working/shared/tpg_gsutil_config.txt
gce_checksum="$(gsutil stat $incominguri/$nwdid/$cramfile |grep 'md5'|awk {'print $3'}|base64 -d|hexdump -ve '/1 "%02x"')"
if [ "$gce_checksum" != "" ]; then        # File does exist in gcepost.  Check MD5
    if [ "$md5" != "$gce_checksum" ]; then
        echo "Local file ($cramfile) checksum ($md5) does not match $gce_checksum for bamid '$bamid' [$nwdid]"
        $topmedcmd -persist mark $bamid $markverb failed
        exit 3
    fi
    echo "MD5 for local and GCE file match ($md5)"
else
  echo ""
  echo "Cannot actually verify the GCE checksum and the MD5 of the file. Assume it is OK"
fi

#   Save date of file in database
$topmedcmd setdate $bamid datemapping_b38 $cramfile

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Post processing of remapped CRAM ($crampath) completed in $etime seconds"
$topmedcmd -persist mark $bamid $markverb completed
$topmedcmd -persist mark $bamid mapped$build completed
$gsutil rm $incominguri/$nwdid/$cramfile.flagstat 
$gsutil rm $incominguri/$nwdid/$cramfile
echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit
