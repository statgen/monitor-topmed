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
mem=16G                         # Really should be 8G, increase to avoid too many at once
markverb="${me}ed"
constraint=''                   # "--constraint eth-10g"
constraint="--constraint eth-10g"       # Force to major nodes, not r63xx nodes
qos=''
slurmp=topmed-working
realhost=''
nodelist=''

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test $me $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  # Run this on node where bam lives
  h=`$topmedpath whathost $1 bam`
  if [ "$h" != "" ]; then
    qos="--qos=$h-$me"
    nodelist="--nodelist=$h"
  fi

  #   Can run anywhere. Low rate of access to cram, small output
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos $nodelist --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem $realhost $constraint $qos $nodelist --workdir=$console -J $1-$me --output=$console/$1-$me.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 $markverb submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save job details
    echo `date` qplot ${l[3]} $slurmp $qos $realhost $mem $realhost >> $console/$1.jobids
  fi
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

#   Is this a cram or bam
extension="${bamfile##*.}"

#   Mark this as started
$topmedcmd mark $bamid $markverb started
d=`date +%Y/%m/%d`
stime=`date +%s`
s=`hostname`
echo "#========= $d host=$s $SLURM_JOB_ID $0 bamid=$bamid bamfile=$bamfile ========="
bai=$bamfile.bai
if [ "$extension" = "cram" ]; then
  bai=$bamfile.crai
fi
if [ ! -f $bai ]; then
  echo "Index '$bai' does not exist"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi

nwdid=`$topmedcmd -persist show $bamid expt_sampleid`
if [ "$nwdid" = "" ]; then
  echo "Unable to find the NWDID for '$bamid'"
  exit 6
fi

#   Create output directory and CD there
outdir=`$topmedpath wherepath $bamid qcresults`
if [ "$outdir" = "" ]; then
  echo "Unable to get QCRESULTS directory for '$bamid' - $outdir"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
mkdir -p $outdir
cd $outdir
if [ "$?" != "0" ]; then
  echo "Unable to CD to qplot output directory for '$bamid' - $outdir"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
echo "Files will be created in $outdir"

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
if [ "$build" = "38" ]; then
  hmkangfiles=/net/fantasia/home/hmkang/code/working/gotcloud_topmed_tmp
  $samtools view -uh -T /data/local/ref/gotcloud.ref/hg38/hs38DH.fa $bamfile | $hmkangfiles/gotcloud/bin/qplot --reference /data/local/ref/gotcloud.ref/hg38/hs38DH.fa --dbsnp /data/local/ref/gotcloud.ref/hg38/dbsnp_142.b38.vcf.gz --stats $basebam.qp.stats --Rcode $basebam.qp.R -.ubam 2>&1
  rc=$?
fi
if [ "$rc" = "none" ]; then
  echo "Unknown build '$build', cannot continue with qplot for '$bamfile'"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
fi
if [ "$rc" != "0" ]; then
  echo "QPLOT failed for '$bamfile'"
  $topmedcmd -persist mark $bamid $markverb failed
  rm -f $basebam.*
  exit 4
fi

#   One last sanity check. Qplot can easily be fooled if samtools truncates
#   when reading the bamfile (NFS surprise).
#   This should fail if the TotalReads < bmflagstat value + 5000 
bamflagstat=`$topmedcmd show $bamid bamflagstat`   # Same for cram or bam
tr=`grep TotalReads $nwdid.src.qp.stats | awk '{ print $2 }'`
tr=`perl -E "print $tr*1000000"`    # Man is hard to do mutiply of float in shell !
tr=`expr $tr + 5001`
if [ "$tr" -lt "$bamflagstat" ]; then
  echo "QPLOT data must have been truncated. TotalReads=$tr  Flagstat=$bamflagstat"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 3
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
  echo "VerifyBAMID failed for '$bamid'"
  $topmedcmd -persist mark $bamid $markverb failed
  rm -f $basebam.*
  exit 5
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "VerifyBAMID on '$bamfile' successful (at second $etime)"
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Command completed in $etime seconds. Created files:"
ls -la $basebam.*

#   Now attempt to put the QPLOT data into the database
$topmedqplot $outdir $nwdid
if [ "$?" != "0" ]; then
  echo "Unable to update the database with the QCPLOT results for '$bamid' [$outdir $nwdid]"
  echo "Maybe try again later with: $topmedqplot $outdir $nwdid"
  $topmedcmd -persist mark $bamid $markverb failed
  exit 7
fi

$topmedcmd -persist mark $bamid $markverb completed
  echo `date` $me $SLURM_JOB_ID ok $etime secs >> $console/$bamid.jobids
exit 0
