#!/bin/bash
#
#   topmed_qplot.sh [-submit ] bamid
#
#	Run QPLOT on a BAM file
#
. /usr/cluster/topmed/bin/topmed_actions.inc

me=qplot
markverb=$me

#------------------------------------------------------------------
# Subroutine:
#   QplotCheck(statsfile, bamid)   #  (statsfile is from Qplot output)
#   Sanity check.  Qplot can easily be fooled if samtools truncates when reading 
#   the bamfile (NFS surprise).  Qplot prints the total number of reads in millions 
#   to 2 decimal places.  We add 5001 in case the last digit was rounded down.
#   This should fail if the TotalReads < bamflagstat value - 5000
#------------------------------------------------------------------
function QplotCheck {
  local statsfile=$1
  local bamid=$2

  bamflagstat=`GetDB $bamid bamflagstat`   # Same for cram or bam
  tr=`grep TotalReads $statsfile | awk '{ print $2 }'`
  if [ "$tr" = "" ]; then
    Fail "QPLOT data $statsfile not found or TotalReads not found"
  fi
  tr=`perl -E "print $tr*1000000"`    # Man, is it hard to do mutiply of float in shell!
  tr=`expr $tr + 5001`
  if [ "$tr" -lt "$bamflagstat" ]; then
    Fail "QPLOT data must have been truncated. TotalReads=$tr Flagstat=$bamflagstat"
  fi
  echo "Qplot output seems reasonable: TotalReads=$tr Flagstat=$bamflagstat"
}

#------------------------------------------------------------------
# Subroutine:
#   RC_Check(rc, basenam, msg)    Check return code and fail with msg
#------------------------------------------------------------------
function RC_Check {
  local rc=$1
  local basenam=$2
  local msg=$3

  if [ "$rc" != "0" ]; then
      rm -f $basenam.*	 	#  we're in $outdir, so this will remove 
      Fail $msg	 	 	    #  bad Qplot or verifyBamId output files
  fi
}  

if [ "$1" = "-submit" ]; then	#  subroutine SubmitJob will run sbatch
  shift	 	 	 	#  with the current shell script and args
  bamid=`GetDB $1 bamid`
  RandomRealHost $bamid
  MayIRun $me $bamid $realhost
  timeout='9:00:00'
  SubmitJob $bamid "topmed" '8G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then	 	#  if no arguments, print usage note
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Run qplot on a bam file and update database"
  exit 1
fi

bamid=$1  	 	 	#  Here starts the actual execution
nwdid=`GetNWDID $bamid`
bamid=`GetDB $nwdid bamid`

#   Mark this as started
SetDB $bamid state_fix 3
Started 

#   Figure out index file name
bamfile=`$topmedpath wherefile $bamid bam`
if [ ! -f $bamfile ]; then      # BAM should exist, but if not, try CRAM
  bamfile=`$topmedpath wherefile $bamid cram`
  echo "BAM not found, using CRAM for '$bamid'"
fi
if [ "$bamfile" = "" ]; then
  Fail "Unable to get source file for '$bamid'; $bamfile"
fi

build=`GetDB $bamid build`
if [ "$build" != "37" -a "$build" != "38" ]; then
  Fail "Unknown build '$build', cannot continue with qplot for '$bamfile'"
fi

#   Create the index file as necessary
CreateIndex $bamid $bamfile

#   Create output directory and CD there
outdir=`$topmedpath wherepath $bamid qcresults`
if [ "$outdir" = "" ]; then
  Fail "Unable to get QCRESULTS directory for '$bamid' - $outdir"
fi
mkdir -p $outdir
cd $outdir
if [ "$?" != "0" ]; then
  Fail "Unable to CD to qplot output directory for '$bamid' - $outdir"
fi
echo "Files will be created in $outdir"
stime=`date +%s`

extension="${bamfile##*.}"
basebam=`basename $bamfile .$extension`
echo "Running qplot for build '$build' extension '$extension'"

#  Define shell variables for Qplot and Fan Zhang's new VerifyBamId
gcref=/net/mario/nodeDataMaster/local/ref/gotcloud.ref
fixverifybamid=/usr/cluster/topmed/bin/nhlbi.1648.vbid.rewrite.awk
verifybamid_dir=/usr/cluster/software/trusty/verify-bam-id
verifybamid=$verifybamid_dir/1.0.0-b38/bin/VerifyBamID
qplotnew=/usr/cluster/topmed/bin/qplot	 	# <== New qplot from Tom

#   Assign resource variables specific to each build
#   The logic here is that we set shell variables specific to build 37 or build 38, then 
#   the same command will run the new qplot for either build and either file format.
#   VerifyBamID requires an 'if-elif-fi' construction, since I want to preserve exactly 
#   the old behavior with Goo Jun's code for any build 37 .bam files that we might re-run.
if [ "$build" = "37" ]; then
   reference=${gcref}/hs37d5.fa
   qplot_snp=${gcref}/dbsnp_142.b37.vcf.gz
   resources=$verifybamid_dir/1.0.0-b38/resource/1000g.100k.b37.vcf.gz.dat
   region="-"
   b37_sites=${gcref}/hapmap_3.3.b37.sites.vcf.gz
elif [ "$build" = "38" ]; then
   reference=${gcref}/hg38/hs38DH.fa
   qplot_snp=${gcref}/hg38/dbsnp_142.b38.vcf.gz
   resources=$verifybamid_dir/1.0.0-b38/resource/1000g.100k.b38.vcf.gz.dat
   region=/net/topmed/incoming/study.reference/study.reference/nhlbi.3066.hs38DH.autosome.list
else
  Fail "Unknown build '$build', cannot continue with qplot for '$bamfile'"
fi

#   Run qplot, either build, output written to current working directory.
#   Same code works for either .bam or .cram
   rm $nwdid.*                          # Remove existing files (not all nicely named)
   $samtools view  -u -T $reference  $bamfile |	 	\
   $qplotnew --reference $reference --dbsnp $qplot_snp	\
   --regions $region  --bamLabels $nwdid	 	\
   --stats $nwdid.qp.stats	 	 	 	\
   --Rcode $nwdid.qp.R  -.ubam  > /run/shm/$nwdid.qp.log 2>&1
   RC_Check $? $nwdid "QPLOT failed for '$bamfile'"
   QplotCheck "$nwdid.qp.stats" $bamid               #  Fails if stats are bad

#   Run verifybamid, output written to current working directory
#   For .cram or build 38 this uses Fan Zhang's new VerifyBamID.
#   For build 37 .bams this uses Goo Jun's original verifyBamId.
#   That writes output file  NWD*.vb.selfSM  directly.  No need 
#   for an awk script.  I don't expect ever to see a .bam format 
#   file on build 38, but Fan's new code will handle it.
if [ "$extension" = "cram" -o "$build" = "38" ]; then	#  Run Fan Zhang's VerifyBamId
  ( $verifybamid	 	 \
     --UDPath ${resources}.UD	 \
     --BedPath ${resources}.bed  \
     --MeanPath ${resources}.mu  \
     --Reference $reference	 \
     --BamFile $bamfile > $nwdid.vb ) 2>/run/shm/$nwdid.vb.log
  RC_Check $? $nwdid "VerifyBAMID failed for '$bamfile'"

  #  Convert this into standard verifybamid output ( $basebam.vb.selfSM ).
  awk -f $fixverifybamid  $nwdid.vb	    # Do not pass NWD identifier
  RC_Check $? $nwdid "$fixverifybamid failed for '$bamfile'"
elif [ "$extension" = "bam" -a "$build" = "37" ]; then	#  Run Goo Jun's old code
  /net/mario/gotcloud/bin/verifyBamID --bam  $bamfile --vcf $b37_sites	 	 	\
    --site --free-full --chip-none --ignoreRG --precise --maxDepth 80	\
    --grid 0.02	--out $nwdid.vb  > /run/shm/$nwdid.vb.log 2>&1
  RC_Check $? $nwdid "VerifyBAMID failed for '$bamfile'"
fi
rm /run/shm/$nwdid.??.log       # Get rid of unused logs

etime=`date +%s`
etime=`expr $etime - $stime`
echo "QPLOT on '$bamfile' successful (at second $etime)"

echo "Created files:"
ls -la $nwdid.*

#   Now attempt to put the QPLOT data into the database
$topmedqplot $outdir $nwdid
if [ "$?" != "0" ]; then
  echo "Maybe try again later with: $topmedqplot $outdir $nwdid"
  Fail "Unable to update the database with the QCPLOT results for '$bamid' [$outdir $nwdid]"
fi

SetDB $bamid state_fix 20           # Keep track of what is redone. Remove someday
Successful
Log $etime
exit 0
