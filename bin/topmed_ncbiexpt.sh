#!/bin/bash
#
#   topmed_ncbiexpt.sh -submit bamid
#
#	Send experiment XML for a bamid to NCBI
#
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
ascpcmd="$topmedcmd send2ncbi"
topmedxml="/usr/cluster/monitor/bin/topmed_xml.pl -master_email ''"
mem=2G
if [ "$TOPMED_MEMORY" != "" ]; then mem=$TOPMED_MEMORY; fi
console=/net/topmed/working/topmed-output
markverb=sentexpt
qos=bai                      # This queue is free most of the time
if [ "$TOPMED_QOS" != "" ]; then qos=$TOPMED_QOS; fi

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test expt $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #   This is short and sweet, can run anywhere
  realhost=topmed
  if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
  slurmp="$realhost-incoming"
  slurmqos="$realhost-$qos"

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-expt --output=$console/$1-sexpt.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-expt --output=$console/$1-sexpt.out $0 $*"
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
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid files=$* ========="
#   Go to our working directory
cd $console
if [ "$?" != "0" ]; then
  echo "Unable to CD to '$console' to create XML file"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
mkdir XMLfiles 2>/dev/null
cd XMLfiles
here=`pwd`

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -type expt $bamid
if [ "$?" != "0" ]; then
  echo "Unable to create experiment XML files"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi

#   Check files we created
files=''
for f in $nwdid-expt.submit.xml $nwdid.expt.xml; do
  if [ ! -f $f ]; then
    echo "Missing XML file '$f'"   
    $topmedcmd mark $bamid $markverb failed
    exit 2
  fi
  files="$files $f"
done

#   Make a TAR of the XML files to send
stime=`date +%s`
tar cf $nwdid-expt.tar $files
if [ "$?" != "0" ]; then
  echo "Unable to create TAR of XML files"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
files=$nwdid-expt.tar

echo "Sending XML files to NCBI - $files"
$ascpcmd $files
if [ "$?" = "0" ]; then
  echo "XML files '$files' sent to NCBI"
  $topmedcmd mark $bamid $markverb delivered
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo `date` expt $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
  exit
fi

echo "FAILED to send XML files to NCBI - $files"
$topmedcmd mark $bamid $markverb failed 
exit 1

