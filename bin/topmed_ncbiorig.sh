#!/bin/bash
#
#   topmed_ncbiorig.sh -submit bamid
#
#	Send the proper set of files to NCBI for the original bams
#
bindir=/usr/cluster/bin
samtools=/net/mario/gotcloud/bin/samtools
ascpcmd="$bindir/ascp -i /net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri -l 800M -k 1"
ascpdest='asp-um-sph@gap-upload.ncbi.nlm.nih.gov:protected'
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
topmedxml=/usr/cluster/monitor/bin/topmed_xml.pl
mem=8G
console=/net/topmed/working/topmed-output
topmeddir=/net/topmed/incoming/topmed

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
  realhost=topmed2              # Force to lesser used machine
  realhostindex="${l[4]}"
  slurmp="$realhost-incoming"
  slurmqos="$realhost-ncbi"

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-ncbiorig --output=$console/$1-ncbiorig.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-ncbiorig --output=$console/$1-ncbiorig.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 sentorig submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobid
    echo `date` ncbiorig ${l[3]} >> $console/$1.jobids
  fi
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid nwdid center path2originalbamfile"
  echo ""
  echo "Send original BAM files to NCBI"
  exit 1
fi
bamid=$1
shift
nwdid=`$topmedcmd show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Invalid bamid '$bamid'. NWDID not known"
  $topmedcmd mark $bamid sentorig failed
  exit 2
fi
center=`$topmedcmd show $bamid center`
if [ "$center" = "" ]; then
  echo "Invalid bamid '$bamid'. CENTER not known"
  $topmedcmd mark $bamid sentorig failed
  exit 2
fi
run=`$topmedcmd show $bamid run`
if [ "$run" = "" ]; then
  echo "Invalid bamid '$bamid'. RUN not known"
  $topmedcmd mark $bamid sentorig failed
  exit 2
fi
origbam=`$topmedcmd show $bamid bamname`
if [ "$origbam" = "" ]; then
  echo "Invalid bamid '$bamid'. BAMNAME not known"
  $topmedcmd mark $bamid sentorig failed
  exit 2
fi
checksum=`$topmedcmd show $bamid checksum`
if [ "$checksum" = "" ]; then
  echo "Invalid bamid '$bamid'. CHECKSUM not known"
  $topmedcmd mark $bamid sentorig failed
  exit 2
fi
fqorigbam=$topmeddir/$center/$run/$origbam
if [ ! -f "$fqorigbam" ]; then
  echo "Original BAM for bamid '$bamid' ($fqorigbam) not found"
  $topmedcmd mark $bamid sentorig failed
  exit 2
fi

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid file=$origbam ========="
#   Go to our working directory
cd $console
if [ "$?" != "0" ]; then
  echo "Unable to CD to '$console' to create XML file"
  $topmedcmd mark $bamid sentorig failed
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
if [ "$center" = "broad" ]; then
  echo "Creating BAM $sendbam from CRAM $fqorigbam"
  $samtools view -b -@3 -o $sendbam $fqorigbam 
  if [ "$?" != "0" ]; then
    echo "Unable to create BAM: $samtools .. $sendbam $fqorigbam"
    $topmedcmd mark $bamid sentorig failed
    exit 2
  fi
else
  ln -fs $fqorigbam $sendbam
  l=`head -1 $sendbam`
  if [ "$l" = "" ]; then
    echo "Unable to find original BAM file: $fqorigbam"
    $topmedcmd mark $bamid sentorig failed
    exit 2
  fi
fi

#   Create the XML to be sent
$topmedxml -xmlprefix $here/ -type secondary $bamid $fqorigbam $checksum
if [ "$?" != "0" ]; then
  echo "Unable to create secondary run XML files"
  $topmedcmd mark $bamid sentorig failed
  exit 2
fi

#   Check files we created.  Note hardcoded 37 here. Fix later
files=''
for f in $nwdid-secondary.37.submit.xml $nwdid-secondary.37.run.xml; do
  if [ ! -f $f ]; then
    echo "Missing XML file '$f'"   
    $topmedcmd mark $bamid sentorig failed
    exit 2
  fi
  files="$files $f"
done

#   Make a TAR of the XML files to send
tar cf $nwdid-secondary.37.tar $files
if [ "$?" != "0" ]; then
  echo "Unable to create TAR of XML files"
  $topmedcmd mark $bamid sentexpt failed
  exit 2
fi
files=$nwdid-secondary.37.tar

echo "Sending data file to NCBI - $sendbam"
ls -l $sendbam
stime=`date +%s`
$ascpcmd $sendbam $ascpdest
rc=$?
rm -f $sendbam
if [ "$rc" != "0" ]; then
  echo "FAILED to send data file '$sendbam'"
  $topmedcmd mark $bamid sentorig failed
  exit 2
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "BAM file sent in $etime seconds"

echo "Sending XML files to NCBI - $files"
$ascpcmd $files $ascpdest
if [ "$?" = "0" ]; then
  echo "XML files '$files' sent to NCBI"
  $topmedcmd mark $bamid sentorig delivered
  exit
fi

echo "FAILED to send data XML files to NCBI - $files"
$topmedcmd mark $bamid sentorig failed 
exit 1




