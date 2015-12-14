#!/bin/bash
#
#   topmed_ncbiexpt.sh -submit bamid
#
#	Send experiment XML for a bamid to NCBI
#
bindir=/usr/cluster/bin
ascpcmd="$bindir/ascp -i /net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri -l 800M -k 1"
ascpdest='asp-um-sph@gap-upload.ncbi.nlm.nih.gov:protected'
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
topmedxml=/usr/cluster/monitor/bin/topmed_xml.pl
mem=2G
console=/net/topmed/working/topmed-output

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test expt $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #   This is short and sweet, can run anywhere
  realhost=topmed
  slurmp="$realhost-incoming"
  slurmqos="$realhost-ncbi"

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-expt --output=$console/$1-sexpt.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-expt --output=$console/$1-sexpt.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 sentexpt submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobid
    echo `date` expt ${l[3]} >> $console/$1.jobids
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
nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Invalid bamid '$bamid'. NWDID not known"
  $topmedcmd mark $bamid sentexpt failed
  exit 2
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid files=$* ========="
#   Go to our working directory
cd $console
if [ "$?" != "0" ]; then
  echo "Unable to CD to '$console' to create XML file"
  $topmedcmd mark $bamid sentexpt failed
  exit 2
fi
mkdir XMLfiles 2>/dev/null
cd XMLfiles
here=`pwd`

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -type expt $bamid
if [ "$?" != "0" ]; then
  echo "Unable to create experiment XML files"
  $topmedcmd mark $bamid sentexpt failed
  exit 2
fi

#   Check files we created
files=''
for f in $nwdid-expt.submit.xml $nwdid.expt.xml; do
  if [ ! -f $f ]; then
    echo "Missing XML file '$f'"   
    $topmedcmd mark $bamid sentexpt failed
    exit 2
  fi
  files="$files $f"
done

#   Make a TAR of the XML files to send
tar cf $nwdid-expt.tar $files
if [ "$?" != "0" ]; then
  echo "Unable to create TAR of XML files"
  $topmedcmd mark $bamid sentexpt failed
  exit 2
fi
files=$nwdid-expt.tar

echo "Sending XML files to NCBI - $files"
$ascpcmd $files $ascpdest
if [ "$?" = "0" ]; then
  echo "XML files '$files' sent to NCBI"
  $topmedcmd mark $bamid sentexpt delivered
  exit
fi

echo "FAILED to send XML files to NCBI - $files"
$topmedcmd mark $bamid sentexpt failed 
exit 1

