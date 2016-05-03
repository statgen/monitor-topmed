#!/bin/bash
#
#   topmed_ncbiorig.sh -submit bamid
#
#	Send the proper set of files to NCBI for the original bams
#
samtools=/net/mario/gotcloud/bin/samtools
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
ascpcmd="$topmedcmd send2ncbi"
topmedxml="/usr/cluster/monitor/bin/topmed_xml.pl"
medir=`dirname $0`
calcmd5=$medir/topmed_calcmd5.sh
mem=8G
if [ "$TOPMED_MEMORY" != "" ]; then mem=$TOPMED_MEMORY; fi
realhost=topmed3
#if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
console=/net/topmed/working/topmed-output
tmpconsole=/working/topmed-output
topmeddir=/net/topmed/incoming/topmed
build=37
version=secondary
markverb=sentorig
jobname=orig
qos=ncbi
if [ "$TOPMED_QOS" != "" ]; then qos=$TOPMED_QOS; fi

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test sorig $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #   Figure where to submit this to run - does not need to be local
  l=(`$topmedcmd where $1 bam`)         # Get pathofbam and host for bam
  h="${l[1]}"
  if [ "$h" != "" ]; then realhost=$h; fi
  if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
  slurmp="$realhost-incoming"
  slurmqos="$realhost-$qos"

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-$jobname --output=$console/$1-$jobname.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-$jobname --output=$console/$1-$jobname.out $0 $*"
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
  echo "Send original BAM files to NCBI"
  exit 1
fi
bamid=$1
shift

$topmedcmd mark $bamid $markverb started
#   Get some values from the database
nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Invalid bamid '$bamid'. NWDID not known"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
center=`$topmedcmd show $bamid center`
if [ "$center" = "" ]; then
  echo "Invalid bamid '$bamid'. CENTER not known"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
origbam=`$topmedcmd show $bamid bamname`
if [ "$origbam" = "" ]; then
  echo "Invalid bamid '$bamid'. BAMNAME not known"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi

#   Go to our working directory
cd $tmpconsole
if [ "$?" != "0" ]; then
  echo "Unable to CD to '$tmpconsole' to create XML file"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
mkdir XMLfiles 2>/dev/null
cd XMLfiles
here=`pwd`

#   Figure out what file to send. Either a BAM or a CRAM
if [ "$center" = "broad" ]; then
  l=(`$topmedcmd where $bamid backup`)      # Get backupdir and backupfile and host
  sendfile="${l[1]}"
  sf=$nwdid.src.cram
  ln -sf $sendfile $sf
  checksum=`$topmedcmd show $bamid cramchecksum`
  sendfile=$sf
  # It's not supposed to happen, but sometimes the CRAM checksum is missing
  if [ "$checksum" = "" ]; then
    echo "Calculating MD5 for CRAM"
    stime=`date +%s`
    checksum=`calcmd5.sh $sendfile | awk '{print $1}'`
    if [ "$checksum" != "" ]; then
      $topmedcmd set $bamid cramb37checksum $checksum
      etime=`date +%s`
      etime=`expr $etime - $stime`
      echo "MD5 calculated for CRAM in $etime seconds"
    fi
  fi
else
  l=(`$topmedcmd where $bamid bam`)         # Get pathofbam and host for bam
  sendfile=$nwdid.src.bam
  ln -sf ${l[0]}/$origbam $sendfile
  checksum=`$topmedcmd show $bamid checksum`
fi
if [ "$checksum" = "" ]; then
  echo "Invalid bamid '$bamid' ($sendfile). CHECKSUM not known"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi

if [ ! -f "$sendfile" ]; then
  echo "Original BAM for bamid '$bamid' ($sendfile) not found"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid file=$sendfile ========="

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -type $version $bamid $sendfile $checksum
if [ "$?" != "0" ]; then
  echo "Unable to create $version run XML files"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi

#   Check files we created
files=''
for f in $nwdid-$version.$build.submit.xml $nwdid-$version.$build.run.xml; do
  if [ ! -f $f ]; then
    echo "Missing XML file '$f'"   
    $topmedcmd mark $bamid $markverb failed
    exit 2
  fi
  files="$files $f"
done

#   Make a TAR of the XML files to send
tar cf $nwdid-$version.$build.tar $files
if [ "$?" != "0" ]; then
  echo "Unable to create TAR of XML files"
  $topmedcmd mark $bamid $markverb failed
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
  echo "FAILED to send data file '$sendfile'"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "File sent in $etime seconds"

echo "Sending XML files to NCBI - $files"
$ascpcmd $files
if [ "$?" = "0" ]; then
  echo "XML files '$files' sent to NCBI"
  $topmedcmd mark $bamid $markverb delivered
  echo `date` $jobname $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
  exit
fi

echo "FAILED to send data XML files to NCBI - $files"
$topmedcmd mark $bamid $markverb failed 
exit 1




