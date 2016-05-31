#!/bin/bash
#
#   topmed_ncbib37.sh -xmlonly -submit bamid
#
#	Send the proper set of files to NCBI for the remapped bam build 37
#
samtools=/net/mario/gotcloud/bin/samtools
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
ascpcmd="$topmedcmd send2ncbi"
topmedxml="/usr/cluster/monitor/bin/topmed_xml.pl"
medir=`dirname $0`
calcmd5=/usr/cluster/monitor/bin/topmed_calcmd5.sh
mem=2G
console=/net/topmed/working/topmed-output
tmpconsole=/working/topmed-output
topmeddir=/net/topmed/incoming/topmed
remapdir="/net/topmed/working/schelcj/results /net/topmed2/working/schelcj/results /net/topmed3/working/schelcj/results /net/topmed4/working/schelcj/results /net/topmed5/working/schelcj/results  /net/topmed6/working/schelcj/results"
build=37
version=remap
markverb=sentb$build
jobname=b$build
xmlonly=N
slurmp=topmed
qos=topmed-ncbi
realhost=''

if [ "$1" = "-xmlonly" ]; then shift; xmlonly=Y; fi    # Force just XML to be sent to NCBI
if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test sb37 $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos --workdir=$console -J $1-$jobname --output=$console/$1-$jobname.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$qos --workdir=$console -J $1-$jobname --output=$console/$1-$jobname.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` $jobname ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-xmlonly] [-submit] bamid"
  echo ""
  echo "Send Remapped Build $build CRAM file to NCBI"
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
piname=`$topmedcmd show $bamid piname`
if [ "$piname" = "" ]; then
  echo "Invalid bamid '$bamid'. PINAME not known"
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

#   Figure out what file to send
sendfile=doesnotexist
for d in $remapdir; do
  dd="$d/$center/$piname/$nwdid/bams"
  if [ -f "$dd/$nwdid.recal.cram" ]; then
    ln -sf $dd/$nwdid.recal.cram $nwdid.recal.$build.cram
    sendfile=$nwdid.recal.$build.cram
    break;
  fi
done
if [ ! -f "$sendfile" ]; then
  echo "Remapped CRAM for Build $build for bamid '$bamid' not found ($sendfile)"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
checksum=`$topmedcmd show $bamid cramb37checksum`
if [ "$checksum" = "" ]; then
  echo "Calculating MD5 for CRAM"
  stime=`date +%s`
  checksum=`$calcmd5 $sendfile | awk '{print $1}'`
  if [ "$checksum" = "" ]; then
    echo "Unable to calculate the MD5 for CRAM: md5sum $sendfile"
    $topmedcmd mark $bamid $markverb failed
    exit 3
  fi
  $topmedcmd set $bamid cramb37checksum $checksum
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "MD5 calculated for CRAM in $etime seconds"
else
  echo "Obtained cramb37checksum for bamid $bamid from database ($checksum)"
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

if [ "$xmlonly" != "Y" ]; then
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
else
  echo "XMLONLY specified, did not send CRAM file"
fi

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




