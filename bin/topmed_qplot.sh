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
verifybamid_dir=/usr/cluster/software/trusty/verify-bam-id
verifybamid=$verifybamid_dir/1.0.0-b38/bin/VerifyBamID
qplot=/usr/cluster/bin/qplot
me=qplot
markverb=$me

#------------------------------------------------------------------
# Subroutine:
#   QplotCheck(statsfile, bamid)
#   Sanity check. Qplot can easily be fooled if samtools truncates
#   when reading the bamfile (NFS surprise).
#   This should fail if the TotalReads < bmflagstat value + 5000
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
      rm -f $basebam.*
      Fail $msg
  fi
}  

if [ "$1" = "-submit" ]; then
  shift
  bamid=`GetDB $1 bamid`
  MayIRun $me $bamid
  RandomRealHost $bamid 'bam'
  SubmitJob $bamid "topmed" '8G' "$0 $*"
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

#   Mark this as started
Started 
nwdid=`GetNWDID $bamid`

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

#   Run qplot.  Right now we have different processes for cram vs bam and which build
extension="${bamfile##*.}"
basebam=`basename $bamfile .$extension`
echo "Running qplot for build '$build' extension '$extension'"

ref37=$gcref/hs37d5.fa
dbsnp37=$gcref/dbsnp_142.b37.vcf.gz
ref38=/data/local/ref/gotcloud.ref/hg38/hs38DH.fa
dbsnp38=/data/local/ref/gotcloud.ref/hg38/dbsnp_142.b38.vcf.gz
hmkangfiles=/net/fantasia/home/hmkang/code/working/gotcloud_topmed_tmp

#===================================================
#  Build 37, bam or cram
#===================================================
if [ "$build" = "37" ]; then
  if [ "$extension" = "bam" ]; then
    #   Run qplot on 37 bam
    $qplot --reference $ref37 --dbsnp $dbsnp37 --label $basebam --stats $basebam.qp.stats \
      --Rcode $basebam.qp.R $bamfile 2>&1
    RC_Check $? $basebam "QPLOT failed for '$bamfile'"
    QplotCheck "$basebam.qp.stats" $bamid               # Fails if stats are bad
    #   Run verifybamid - use old version
    $gcbin/verifyBamID --bam  $bamfile --vcf $gcref/hapmap_3.3.b37.sites.vcf.gz	\
      --site --free-full --chip-none --ignoreRG --precise --maxDepth 80 \
      --grid 0.02 --out $basebam.vb 2>&1
    #   New verifybamid does not work with b37 bam
    #verifybamid_res=$verifybamid_dir/src/VerifyBamID/resource
    #verifybamid_dat=$verifybamid_res/1000g.100k.b37.vcf.gz.dat
    #$verifybamid --UDPath ${verifybamid_dat}.UD --BedPath ${verifybamid_dat}.bed \
    #  --MeanPath ${verifybamid_dat}.mu --Reference $ref37 --BamFile $bamfile > $basebam.vb.raw 2>&1
echo "Not doing: $verifybamid --UDPath ${verifybamid_dat}.UD --BedPath ${verifybamid_dat}.bed --MeanPath ${verifybamid_dat}.mu --Reference $ref37 --BamFile $bamfile > $basebam.vb.raw"
    RC_Check $? $basebam "VerifyBAMID failed for '$bamfile'"
  else
    #   Run qplot on 37 cram
    $samtools view -uh -T $ref37 $bamfile | $hmkangfiles/gotcloud/bin/qplot --reference $ref37 \
      --stats $basebam.qp.stats --Rcode $basebam.qp.R -.ubam 2>&1
    RC_Check $? $basebam "QPLOT failed for '$bamfile'"
    QplotCheck "$basebam.qp.stats" $bamid               # Fails if stats are bad
    #   Run verifybamid
    verifybamid_res=$verifybamid_dir/src/VerifyBamID/resource
    verifybamid_dat=$verifybamid_res/1000g.100k.b37.vcf.gz.dat
    $verifybamid --UDPath ${verifybamid_dat}.UD --BedPath ${verifybamid_dat}.bed \
      --MeanPath ${verifybamid_dat}.mu --Reference $ref37 --BamFile $bamfile > $basebam.vb
    RC_Check $? $basebam "VerifyBAMID failed for '$bamfile'"
      #   Convert this verifybamid output ($basebam.vb.selfSM) like all others
    awk -f $fixverifybamid -v NWD=$nwdid $basebam.vb
    RC_Check $? $basebam "$fixverifybamid failed for '$bamfile'"
  fi
fi

#===================================================
#  Build 38, bam or cram
#===================================================
if [ "$build" = "38" ]; then
  if [ "$extension" = "bam" ]; then
    #   Run qplot on 38 bam
    $qplot --reference $ref38 --dbsnp $dbsnp38 --label $basebam --stats $basebam.qp.stats \
      --Rcode $basebam.qp.R $bamfile 2>&1
    RC_Check $? $basebam "QPLOT failed for '$bamfile'"
    QplotCheck "$basebam.qp.stats" $bamid               # Fails if stats are bad
    #   Run verifybamid
    verifybamid_res=$verifybamid_dir/src/VerifyBamID/resource
    verifybamid_dat=$verifybamid_res/1000g.100k.b38.vcf.gz.dat
    $verifybamid --UDPath ${verifybamid_dat}.UD --BedPath ${verifybamid_dat}.bed \
      --MeanPath ${verifybamid_dat}.mu --Reference $ref38 --BamFile $bamfile > $basebam.vb
    RC_Check $? $basebam "VerifyBAMID failed for '$bamfile'"
  else
    #   Run qplot on 38 cram
    $samtools view -uh -T $ref38 $bamfile | $qplot --reference $ref38 --dbsnp $dbsnp38 --stats $basebam.qp.stats --Rcode $basebam.qp.R -.ubam 2>&1
    RC_Check $? $basebam "QPLOT failed for '$bamfile'"
    QplotCheck "$basebam.qp.stats" $bamid               # Fails if stats are bad
    #   Run verifybamid
    verifybamid_res=$verifybamid_dir/src/VerifyBamID/resource
    verifybamid_dat=$verifybamid_res/1000g.100k.b38.vcf.gz.dat
    $verifybamid --UDPath ${verifybamid_dat}.UD --BedPath ${verifybamid_dat}.bed \
      --MeanPath ${verifybamid_dat}.mu --Reference $ref38 --BamFile $bamfile > $basebam.vb
    RC_Check $? $basebam "VerifyBAMID failed for '$bamfile'"
      #   Convert this verifybamid output ($basebam.vb.selfSM) like all others
    awk -f $fixverifybamid -v NWD=$nwdid $basebam.vb
    RC_Check $? $basebam "$fixverifybamid failed for '$bamfile'"
  fi
fi

etime=`date +%s`
etime=`expr $etime - $stime`
echo "QPLOT on '$bamfile' successful (at second $etime)"

echo "Created files:"
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
