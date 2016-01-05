#!/bin/bash
#
#   topmed_ncbib37.sh -submit [-xmlonly] bamid
#
#   e.g. topmed_ncbib37.sh -xmlonly 11303   (for nygc 11302 NWD792235)
#
#	Send the proper set of files to NCBI for the remap b37 bam
#
bindir=/usr/cluster/bin
samtools=/net/mario/gotcloud/bin/samtools
ascpcmd="$bindir/ascp -i /net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri -l 800M -k 1"
ascpdest='asp-um-sph@gap-submit.ncbi.nlm.nih.gov:protected'
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
topmedxml=/usr/cluster/monitor/bin/topmed_xml.pl
mem=8G
if [ "$TOPMED_MEMORY" != "" ]; then mem=$TOPMED_MEMORY; fi
console=/net/topmed/working/topmed-output
tmpconsole=/working/topmed-output
topmeddir=/net/topmed/incoming/topmed
remapdir=/net/topmed/working/schelcj/results
build=37
version=primary
markverb=sentb$build
jobname=b$build
qos=ncbi
if [ "$TOPMED_QOS" != "" ]; then qos=$TOPMED_QOS; fi

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test sb37 $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #   This will usually have files on both topmeds, so it does not matter where it runs
  l=(`$topmedcmd where $1`)     # Get bampath backuppath bamname realhost realhostindex
  realhost="${l[3]}"
  realhost=topmed3
  if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
  realhostindex="${l[4]}"
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
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` $jobname ${l[3]} $slurmp $slurmqos $mem >> $console/$1.jobids
  fi
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] [-xmlonly] bamid"
  echo ""
  echo "Send Remapped Build $build BAM file to NCBI"
  exit 1
fi
send=all
if [ "$1" = "-xmlonly" ]; then      # We need to resend the XML, not the file
  shift
  send=xmlonly
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
fqremappedcram="$remapdir/$center/$piname/$nwdid/bams/$nwdid.recal.cram"
if [ ! -f "$fqremappedcram" ]; then
  echo "Remapped CRAM for bamid '$bamid' ($fqremappedcram) not found"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid file=$fqremappedcram ========="
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

#   Create the BAM file to be sent
#   samtools view by default writes to std out.  Could use option "-o *.bam" to
#   suppress this and specify the output file name in the options to samtools.
#   Option "-@" tells how many extra threads to use for compression when writing
#   a .bam file.  I'm just guessing 3 -- meaning 4 threads total, including the main one.
sendbam="$nwdid.recal.$build.bam"
if [ "$send" != "xmlonly" ]; then
  echo "Creating BAM $sendbam from CRAM $fqremappedcram"
  stime=`date +%s`
  $samtools view -b -@3 -o $sendbam $fqremappedcram 
  if [ "$?" != "0" ]; then
    echo "Unable to create BAM: $samtools .. $sendbam $fqremappedcram"
    $topmedcmd mark $bamid $markverb failed
    exit 2
  fi
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "BAM created from CRAM in $etime seconds"

  #   Calculate the MD5 for the bam
  echo "Calculating MD5 for new bam"
  stime=`date +%s`
  checksum=`md5sum $sendbam | awk '{print $1}'`
  if [ "$checksum" = "" ]; then
    echo "Unable to calculate the MD5 for the new BAM: md5sum $sendbam"
    $topmedcmd mark $bamid $markverb failed
    exit 3
  fi
  $topmedcmd set $bamid b37bamchecksum $checksum
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "MD5 calculated from CRAM in $etime seconds"
else
  echo "XMLONLY - $sendbam not created"
  checksum=`$topmedcmd show $bamid b37bamchecksum`
  if [ "$checksum" = "" ]; then
    echo "Invalid bamid '$bamid'. b37bamchecksum not known"
    $topmedcmd mark $bamid $markverb failed
    exit 2
  fi
fi

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -type remap $bamid $here/$sendbam $checksum
if [ "$?" != "0" ]; then
  echo "Unable to create remap run XML files"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi

#   Check files we created
files=''
for f in $nwdid-remap.$build.submit.xml $nwdid-remap.$build.run.xml; do
  if [ ! -f $f ]; then
    echo "Missing XML file '$f'"   
    $topmedcmd mark $bamid $markverb failed
    exit 2
  fi
  files="$files $f"
done

#   Make a TAR of the XML files to send
tar cf $nwdid-remap.$build.tar $files
if [ "$?" != "0" ]; then
  echo "Unable to create TAR of XML files"
  $topmedcmd mark $bamid sentexpt failed
  exit 2
fi
files=$nwdid-remap.$build.tar

if [ "$send" != "xmlonly" ]; then
  echo "Sending data file to NCBI - $sendbam"
  ls -l $sendbam
  stime=`date +%s`
  $ascpcmd $sendbam $ascpdest
  rc=$?
  rm -f $sendbam
  if [ "$rc" != "0" ]; then
    echo "FAILED to send data file '$sendbam'"
    $topmedcmd mark $bamid $markverb failed
    exit 2
  fi
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "BAM file sent in $etime seconds"
else
  echo "XMLONLY - $sendbam not sent"
fi

echo "Sending XML files to NCBI - $files"
$ascpcmd $files $ascpdest
if [ "$?" = "0" ]; then
  echo "XML files '$files' sent to NCBI"
  $topmedcmd mark $bamid $markverb delivered
  exit
fi

echo "FAILED to send data XML files to NCBI - $files"
$topmedcmd mark $bamid $markverb failed 
exit 1




