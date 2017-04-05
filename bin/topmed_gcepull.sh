#!/bin/bash
#
#   topmed_gcepull.sh -submit| bamid
#
#	Copy remapped CRAM for a sample from GCE
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=gcepull
mem=2G
markverb=$me
build=38
qos="--qos=topmed-$me"
realhost=''

cores="--cpus-per-task=2"         # Cores here should be same as gsutil
gsutil='gsutil -o GSUtil:parallel_composite_upload_threshold=150M -o GSUtil:parallel_process_count=2'
incominguri='gs://topmed-recabs'
bcfuri='gs://topmed-bcf'

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedpermit permit test $me $1
  if [ "$?" = "0" ]; then
    echo "$me $1 not permitted" | tee $console/$1-$me.out
    exit 4
  fi 

  # Run this on node where remapped cram will live
  h=`$topmedpath whathost $1 b$build`
  if [ "$h" != "" ]; then
    realhost="--nodelist=$h"
    #qos="--qos=$h-$me"
  fi

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $cores $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*`)
  if [ "$?" != "0" ]; then
    $topmedcmd mark $1 $markverb failed
    echo "Failed to submit command to SLURM - $l" > $console/$1-$me.out
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $cores $qos --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $sq $*" >> $console/$1-$me.out
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
  echo "Copy remapped CRAM for a sample from GCE"
  exit 1
fi
bamid=$1

#   Get remapped cram file location
crampath=`$topmedpath wherepath $bamid b$build`
if [ "$crampath" = "" ]; then
  echo "Unable to determine where remapped CRAM file should go for '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Unable to get expt_sampleid for bamid '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
cramfile=$crampath/$nwdid.recab.cram

d=`date +%Y/%m/%d`
s=`hostname`
p=`pwd`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid cramfile=$cramfile pwd=$p ========="

#   Mark this as started
$topmedcmd mark $bamid $markverb started

mkdir -p $crampath
if [ "$?" != "0" ]; then
  echo "Unable to create '$crampath' for remapped CRAM for '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

#======================================================================
#   Copy remapped CRAM from GCE, check flagstat, fix up database
#======================================================================
cramflagstat=`$topmedcmd show $bamid cramflagstat`
if [ "$cramflagstat" = "" ]; then
  echo "Unable to get cramflagstat for bamid '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

stime=`date +%s`
export BOTO_CONFIG=/net/topmed/working/shared/tpg_gsutil_config.txt
mkdir -p $crampath
if [ "$?" != "0" ]; then
  echo "Unable to create directory for remapped data '$crampath' for bamid '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi

#   Remapped cram could be > one place (arrgh!)  Figure out where it is
inuri=''
for i in $incominguri $bcfuri; do
  p="$i/$nwdid/$nwdid.recab.cram"
  $gsutil stat "$p"
  if [ "$?" = "0" ]; then
    inuri=$i
    break
  fi
done
if [ "$inuri" = "" ]; then
  echo "Unable to find $nwdid/$nwdid.recab.cram in: $incominguri $bcfuri"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi

#   Now know where to look for data. Check flagstat
echo "Checking if flagstat is as we expect from $inuri"
$gsutil cp  $inuri/$nwdid/$nwdid.recab.cram.flagstat $crampath
if [ "$?" != "0" ]; then
  echo "Failed to copy flagstat from GCE: $inuri/$nwdid/$nwdid.recab.cram.flagstat"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
#   Get number of interest from flagstat file and check it
n=`grep 'paired in sequencing' $crampath/$nwdid.recab.cram.flagstat | awk '{print $1}'`
if [ "$n" != "$cramflagstat" ]; then
  echo "Flagstat '$n' did not match cramflagstat '$cramflagstat' for bamid '$bamid' nwdid $nwdid"
  echo "URL=$inuri"
  $topmedcmd -persist mark $bamid $markverb failed
  # Renaming the flagstat file stops pull from happening again
  $gsutil mv $inuri/$nwdid/$nwdid.recab.cram.flagstat $inuri/$nwdid/$nwdid.recab.cram.flagstat.nomatch
  exit 3
fi
echo "Flagstat value is correct: $n"

#   See if we have already done this
f=$crampath/$nwdid/$nwdid.recab.cram
if [ -f $f ]; then
  echo "Replacing existing CRAM $f"
fi

echo "Copying remapped CRAM to local file $crampath"
$gsutil cp $inuri/$nwdid/$nwdid.recab.cram $crampath
if [ "$?" != "0" ]; then
  echo "Failed to copy file from GCE $inuri/$nwdid/$nwdid.recab.cram $crampath"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi

#   Post processing needed here
echo "Calculating MD5 for local file ($cramfile)"
md5=(`md5sum $cramfile`)
md5=${md5[0]}
if [ "$md5" = "" ]; then
  echo "Unable to calculate MD5 for remapped '$bamid' [$nwdid] cramfile=$cramfile"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
$topmedcmd -persist set $bamid cramb38checksum $md5
echo "Set checksum for b$build file"
$topmedcmd -persist set $bamid b${build}flagstat $cramflagstat

#   Save date of file in database
$topmedcmd setdate $bamid datemapping_b38 $cramfile

#   Clean up data in GCE if data found in incoming.  Move remapped data to bcf bucket for now
if [ "$inuri" = "$incominguri" ]; then
  $gsutil mv $incominguri/$nwdid/$cramfile.flagstat  $bcfuri/$nwdid/$cramfile.flagstat
  $gsutil mv $incominguri/$nwdid/$cramfile           $bcfuri/$nwdid/$cramfile
  echo "Moved $incominguri/$nwdid to $bcfuri/$nwdid"
  #   Remove any left over cruft in recabs bucket
  echo "Removing $incominguri/$nwdid"
  $gsutil rm -r $incominguri/$nwdid
else
  echo "Data was not found in $incominguri so we leave it where it was"
fi

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Copy of remapped CRAM from GCE to $crampath completed in $etime seconds"

$topmedcmd -persist mark $bamid $markverb completed
$topmedcmd -persist mark $bamid gcepost completed
echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit
