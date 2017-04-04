#!/bin/bash
#
#   topmed_ncbiexpt.sh -submit bamid
#
#	Send experiment XML for a bamid to NCBI
#
. /usr/cluster/topmed/bin/topmed_actions.inc
topmedxml="/usr/cluster/topmed/bin/topmed_xml.pl"
ascpcmd="$topmedcmd send2ncbi"

mem=2G
markverb=sentexpt
slurmp=topmed-working
qos=topmed-redo
realhost=''
realhost="--nodelist=topmed"       # Force to machine with external interface

#   Do not allow this to play with anything except year one data
year=`$topmedcmd show $bamid datayear`
if [ "$year" != "1" ]; then
  echo "$* must be year ONE data"
  exit 4
fi

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  #$topmedcmd permit test sexpt $1
  #if [ "$?" = "0" ]; then
  #  exit 4
  #fi 

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos $realhost --workdir=$console -J $1-expt --output=$console/$1-sexpt.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos $realhost --workdir=$console -J $1-expt --output=$console/$1-sexpt.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` expt ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
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
shift

$topmedcmd mark $bamid $markverb started

nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Invalid bamid '$bamid'. NWDID not known"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid files=$* ========="
#   Go to our working directory
cd $console
if [ "$?" != "0" ]; then
  echo "Unable to CD to '$console' to create XML file"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
mkdir XMLfiles 2>/dev/null
cd XMLfiles
here=`pwd`

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -type expt $bamid
if [ "$?" != "0" ]; then
  echo "Unable to create experiment XML files"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

#   Check files we created
files=''
for f in $nwdid-expt.submit.xml $nwdid.expt.xml; do
  if [ ! -f $f ]; then
    echo "Missing XML file '$f'"   
    $topmedcmd -persist mark $bamid $markverb failed
    exit 2
  fi
  files="$files $f"
done

#   Make a TAR of the XML files to send
stime=`date +%s`
tar cf $nwdid-expt.tar $files
if [ "$?" != "0" ]; then
  echo "Unable to create TAR of XML files"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
files=$nwdid-expt.tar

echo "Sending XML files to NCBI - $files"
$ascpcmd $files
if [ "$?" = "0" ]; then
  echo "XML files '$files' sent to NCBI"
  $topmedcmd -persist mark $bamid $markverb delivered
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo `date` expt $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
  exit
fi

echo "FAILED to send XML files to NCBI - $files"
$topmedcmd -persist mark $bamid $markverb failed 
exit 1

