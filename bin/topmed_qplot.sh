#!/bin/bash
#
#   topmed_qplot.sh [-submit [-mem NG]] bamid
#
#	Run QPLOT on a BAM file
#
. /usr/cluster/topmed/bin/topmed_actions.inc

gcbin=/net/mario/gotcloud/bin
gcref=/net/mario/nodeDataMaster/local/ref/gotcloud.ref
topoutdir=/net/topmed/incoming/qc.results
fixverifybamid=/usr/cluster/topmed/bin/nhlbi.1648.vbid.rewrite.awk

me=qplot
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me $bamid
  MyRealHost $bamid 'bam'
  SubmitJob $bamid "$realhost-$me" '16G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid"
  echo ""
  echo "Run qplot on a bam file and update database"
  exit 1
fi
bamid=$1
bamfile=`$topmedpath wherefile $bamid bam`
extension="${bamfile##*.}"

#   Mark this as started
Started 
GetNWDID $bamid

#   Figure out index file name
bai=$bamfile.bai
if [ "$extension" = "cram" ]; then
  bai=$bamfile.crai
fi

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

#   If necessary, create the index file
if [ -f $bai ]; then
  echo "Using existing index file '$bai'"
else
  echo "Creating index file '$bai'"
  $samtools index $bamfile 2>&1
  if [ "$?" != "0" ]; then
    Fail "Unable to create index file for '$bamfile' in '$outdir'"
  fi
fi

#   Run qplot, output written to current working directory
basebam=`basename $bamfile .$extension`
build=`$topmedcmd -persist show $bamid build`
rc=none
echo "Running qplot for build '$build'"
if [ "$build" = "37" ]; then
  $gcbin/qplot --reference  $gcref/hs37d5.fa --dbsnp $gcref/dbsnp_142.b37.vcf.gz \
  --label $basebam --stats $basebam.qp.stats --Rcode $basebam.qp.R $bamfile 2>&1
  rc=$?
fi
if [ "$build" = "38" ]; then        # Someday this will fail because it got moved
  hmkangfiles=/net/fantasia/home/hmkang/code/working/gotcloud_topmed_tmp
  $samtools view -uh -T /data/local/ref/gotcloud.ref/hg38/hs38DH.fa $bamfile | $hmkangfiles/gotcloud/bin/qplot --reference /data/local/ref/gotcloud.ref/hg38/hs38DH.fa --dbsnp /data/local/ref/gotcloud.ref/hg38/dbsnp_142.b38.vcf.gz --stats $basebam.qp.stats --Rcode $basebam.qp.R -.ubam 2>&1
  rc=$?
fi
if [ "$rc" = "none" ]; then
  Fail "Unknown build '$build', cannot continue with qplot for '$bamfile'"
fi
if [ "$rc" != "0" ]; then
  rm -f $basebam.*
  Fail "QPLOT failed for '$bamfile'"
fi

#   One last sanity check. Qplot can easily be fooled if samtools truncates
#   when reading the bamfile (NFS surprise).
#   This should fail if the TotalReads < bmflagstat value + 5000
statsfile=
bamflagstat=`$topmedcmd show $bamid bamflagstat`   # Same for cram or bam
tr=`grep TotalReads $basebam.qp.stats | awk '{ print $2 }'`
if [ "$tr" = "" ]; then
  Fail "QPLOT data $nwdid.src.qp.stats not found or TotalReads not found"
fi
tr=`perl -E "print $tr*1000000"`    # Man is hard to do mutiply of float in shell !
tr=`expr $tr + 5001`
if [ "$tr" -lt "$bamflagstat" ]; then
  Fail "QPLOT data must have been truncated. TotalReads=$tr  Flagstat=$bamflagstat"
fi
echo "Qplot output seems reasonable: TotalReads=$tr  Flagstat=$bamflagstat"

etime=`date +%s`
etime=`expr $etime - $stime`
echo "QPLOT on '$bamfile' successful (at second $etime)"

#   Run verifybamid, output written to current working directory
#   Notice the special case for a cram.  Very non-production looking
if [ "$extension" = "cram" ]; then
  hmkangfiles=/net/fantasia/home/hmkang/code/working/gotcloud_topmed_tmp
  $hmkangfiles/contamination-finder/build/ContaminationFinder --UDPath $hmkangfiles/resources/1000g.10k.vcf.gz.dat.UD --BedPath $hmkangfiles/resources/10k.build38.bed --MeanPath $hmkangfiles/resources/1000g.10k.vcf.gz.dat.mu --Reference /data/local/ref/gotcloud.ref/hg38/hs38DH.fa --BamFile $bamfile > $basebam.vb
  rc=$?
  #   Convert this verifybamid output ($basebam.vb.selfSM ) like all others
  if [ "$rc" = "0" ]; then
    awk -f $fixverifybamid -v NWD=$nwdid $basebam.vb
    rc=$?
  fi
else
  $gcbin/verifyBamID --bam  $bamfile --vcf $gcref/hapmap_3.3.b37.sites.vcf.gz	\
    --site --free-full --chip-none --ignoreRG --precise --maxDepth 80 \
    --grid 0.02	--out  $basebam.vb 2>&1
    rc=$?
fi
if [ "$rc" != "0" ]; then
  rm -f $basebam.*
  Fail "VerifyBAMID failed for '$bamid'"
fi
etime=`date +%s`
echo "VerifyBAMID on '$bamfile' successful (at second $etime)"
etime=`expr $etime - $stime`
echo "Command completed in $etime seconds. Created files:"
ls -la $basebam.*

#   Now attempt to put the QPLOT data into the database
$topmedqplot $outdir $nwdid
if [ "$?" != "0" ]; then
  echo "Maybe try again later with: $topmedqplot $outdir $nwdid"
  Fail "Unable to update the database with the QCPLOT results for '$bamid' [$outdir $nwdid]"
fi

Successful
Log $etime
exit 0
