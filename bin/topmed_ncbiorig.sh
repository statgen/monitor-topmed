#!/bin/bash
#
#   topmed_ncbiorig.sh -submit [-xmlonly] bamid
#
#	Send the proper set of files to NCBI for the original bams
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

  #   This will usually have files on both topmeds, so it does not matter where it runs
  l=(`$topmedcmd where $1`)     # Get bampath backuppath bamname realhost realhostindex
  realhost="${l[3]}"
  realhost=topmed
  if [ "$TOPMED_HOST" != "" ]; then realhost=$TOPMED_HOST; fi
  realhostindex="${l[4]}"
  slurmp="$realhost-incoming"
  slurmqos="$realhost-$qos"

  #  Submit this script to be run
  l=(`$bindir/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-$jobname --output=$console/$1-$jobname.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=$bindir/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-$jobname --output=$console/$1-$jobname.out $0 $*"
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
  echo "Usage: $me [-submit] [-xmlonly] bamid"
  echo ""
  echo "Send original BAM files to NCBI"
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
run=`$topmedcmd show $bamid run`
if [ "$run" = "" ]; then
  echo "Invalid bamid '$bamid'. RUN not known"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
origbam=`$topmedcmd show $bamid bamname`
if [ "$origbam" = "" ]; then
  echo "Invalid bamid '$bamid'. BAMNAME not known"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
checksum=`$topmedcmd show $bamid checksum`
if [ "$checksum" = "" ]; then
  echo "Invalid bamid '$bamid'. CHECKSUM not known"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi
fqorigbam=$topmeddir/$center/$run/$origbam
if [ ! -f "$fqorigbam" ]; then
  echo "Original BAM for bamid '$bamid' ($fqorigbam) not found"
  $topmedcmd mark $bamid $markverb failed
  exit 2
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid file=$origbam ========="
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

#   If necessary, create the BAM file to be sent
#   samtools view by default writes to std out.  Could use option "-o *.bam" to
#   suppress this and specify the output file name in the options to samtools.
#   Option "-@" tells how many extra threads to use for compression when writing
#   a .bam file.  I'm just guessing 3 -- meaning 4 threads total, including the main one.
sendbam="$nwdid.src.bam"
if [ "$send" != "xmlonly" ]; then
  if [ "$center" = "broad" ]; then
    echo "Creating BAM $sendbam from CRAM $fqorigbam"
    $samtools view -b -@3 -o $sendbam $fqorigbam 
    if [ "$?" != "0" ]; then
      echo "Unable to create BAM: $samtools .. $sendbam $fqorigbam"
      $topmedcmd mark $bamid $markverb failed
      exit 2
    fi
  else
    ln -fs $fqorigbam $sendbam
    l=`head -1 $sendbam`
    if [ "$l" = "" ]; then
      echo "Unable to find original BAM file: $fqorigbam"
      $topmedcmd mark $bamid $markverb failed
      exit 2
    fi
  fi
else
  echo "XMLONLY - $sendbam not created"
fi

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -type $version $bamid $fqorigbam $checksum
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




