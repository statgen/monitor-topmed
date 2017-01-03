#!/bin/bash
#
#   topmed_ncbiorig.sh -submit bamid
#
#	Send the proper set of files to NCBI for the original/secondary files
#
samtools=/usr/cluster/bin/samtools
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
topmedpath=/usr/cluster/monitor/bin/topmedpath.pl
ascpcmd="$topmedcmd send2ncbi"
topmedxml=/usr/cluster/monitor/bin/topmed_xml.pl
medir=`dirname $0`
calcmd5=/usr/cluster/monitor/bin/topmed_calcmd5.sh
mem=2G
console=/net/topmed/working/topmed-output
tmpconsole=/net/topmed/working/topmed-output
topmeddir=/net/topmed/incoming/topmed
build=37
version=secondary
markverb=sentorig
jobname=orig
slurmp=topmed
qos=topmed-ncbi
realhost=''
realhost="--nodelist=topmed"       # Force to machine with external interface

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test sorig $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos $realhost --workdir=$console -J $1-$jobname --output=$console/$1-$jobname.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos $realhost --workdir=$console -J $1-$jobname --output=$console/$1-$jobname.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobdetails
    echo `date` $jobname ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
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
shift

$topmedcmd mark $bamid $markverb started
#   Get some values from the database
nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Invalid bamid '$bamid'. NWDID not known"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
center=`$topmedcmd show $bamid center`
if [ "$center" = "" ]; then
  echo "Invalid bamid '$bamid'. CENTER not known"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
build=`$topmedcmd show $bamid build`
if [ "$build" = "" ]; then
  echo "Invalid build '$build'. BUILD not known"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
origbam=`$topmedcmd show $bamid bamname`
if [ "$origbam" = "" ]; then
  echo "Invalid bamid '$bamid'. BAMNAME not known"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

#   Go to our working directory
cd $tmpconsole
if [ "$?" != "0" ]; then
  echo "Unable to CD to '$tmpconsole' to create XML file"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
mkdir XMLfiles 2>/dev/null
cd XMLfiles
here=`pwd`

#   Figure out what file to send. Either a BAM or a CRAM
sendcram='N'
offsite=`$topmedcmd show $bamid offsite`    # Original file is offsite or local
datayear=`$topmedcmd show $bamid datayear`  # What year data is this
if [ "$offsite" = "D" ]; then
  sendcram='Y'
else
  if [ "$center" = "broad" -a "$datayear" != '1' ]; then
    sendcram='Y'
  fi
fi

if [ "$sendcram" = "Y" ]; then
  l=(`$topmedpath wherefile $bamid backup`)        # Get backupdir and backupfile and host
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
  l=(`$topmedpath wherefile $bamid bam`)       # Get pathofbam for bam
  sendfile=$nwdid.src.bam
  ln -sf ${l[0]} $sendfile
  checksum=`$topmedcmd -persist show $bamid checksum`
fi

if [ "$checksum" = "" ]; then
  echo "Invalid bamid '$bamid' ($sendfile). CHECKSUM not known"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

if [ ! -f "$sendfile" ]; then
  echo "Original BAM for bamid '$bamid' ($sendfile) not found"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid file=$sendfile ========="

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -build $build -type $version $bamid $sendfile $checksum
if [ "$?" != "0" ]; then
  echo "Unable to create $version run XML files"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi

#   Check files we created
files=''
for f in $nwdid-$version.$build.submit.xml $nwdid-$version.$build.run.xml; do
  if [ ! -f $f ]; then
    echo "Missing XML file '$f'"   
    $topmedcmd -persist mark $bamid $markverb failed
    exit 2
  fi
  files="$files $f"
done

#   Make a TAR of the XML files to send
tar cf $nwdid-$version.$build.tar $files
if [ "$?" != "0" ]; then
  echo "Unable to create TAR of XML files"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
files=$nwdid-$version.$build.tar

echo "Sending data file to NCBI - $sendfile"
ls -l $sendfile
stime=`date +%s`
$ascpcmd $sendfile
rc=$?
rm -f $sendfile
if [ "$rc" != "0" ]; then
  echo "FAILED to send data file '$sendfile' (rc=$rc)"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 2
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "File sent in $etime seconds"

echo "Sending XML files to NCBI - $files"
$ascpcmd $files
if [ "$?" = "0" ]; then
  echo "XML files '$files' sent to NCBI"
  $topmedcmd -persist mark $bamid $markverb delivered
  echo `date` $jobname $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
  exit
fi

echo "FAILED to send data XML files to NCBI - $files"
$topmedcmd mark $bamid $markverb failed 
exit 1




