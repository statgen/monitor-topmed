#!/bin/bash
#
#   topmed_fix.sh -submit bamid|nwdid
#
#	Fix some sort of problem. This script changes all the time
#
. /usr/cluster/topmed/bin/topmed_actions.inc
me=fix
markverb=fix

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  #MayIRun $me  $bamid
  #h=(`$topmedpath wherepath $bamid b38 | sed -e 's:/: :g'`)
  #realhost=${h[1]}
  #SubmitJob $bamid "$realhost-fix" '8G' "$0 $*"
  RandomRealHost $bamid
  SubmitJob $bamid "topmed-fix" '8G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid|nwdid"
  echo ""
  echo "Fix some sort of problem."
  exit 1
fi

bamid=$1
build=b38

Started

echo "# $0 $bamid $nwdid  -- redo header for mis-mapped b38 crams"
recab=`$topmedpath wherefile $bamid b38`
tmpdir=/tmp
fixlog=$console/fix.log

#  Create a corrected cram
b=`basename $recab`
newrecab=$tmpdir/$b
nwdid=`GetDB $bamid expt_sampleid`

rgmap=`ls /net/topmed8/working/call_sets/rehdr/rgmaps/*/$nwdid.rgmap`
if [ "$rgmap" = "" ]; then
  Fail "Unable to find $rgmap"
fi

stime=`date +%s`
/net/fantasia/home/hmkang/tools/apigenome.test/bin/cramore cram-update-rg --ref /net/mario/nodeDataMaster/local/ref/gotcloud.ref/hg38/hs38DH.fa --rg-map $rgmap --sam $recab --out $newrecab
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Corrected recab completed in $etime seconds"
if [ "$rc" != "0" ]; then
  echo "Corrected recab for $bamid build $build failed" >> $fixlog
  Fail "Corrected recab failed" 
fi

#   Get flagstat values for new file, must match original flagstat
echo "Calculating flagstat for $newrecab"
b=`basename $newrecab`
outfile=/tmp/$b.flagstat
$samtools flagstat  $newrecab > $outfile
if [ "$?" != "0" ]; then
  Fail "Unable to calculate flagstat for $newrecab"
fi
a=(`grep 'paired in sequencing' $outfile`)   #   Get paired reads from flagstat
newflagstat=${a[0]}
rm -f $outfile
if [ "$newflagstat" = "" ]; then
  Fail "Unable to get reads of paired in sequencing from '$output' for '$newrecab'"
fi

origflagstat=`GetDB $bamid b38flagstat`
if [ "$origflagstat" != "$newflagstat" ]; then
  Fail "New remapped file flagstat incorrect $origflagstat != $newflagstat"
fi
echo "New remapped file is valid, flagstats match ($newflagstat)"

#  Create a new index for the file
$topmedmakeindex $newrecab $build $console/$bamid-$me.out
if [ "$?" != "0" ]; then
  Fail "Unable to create index file for '$recab'"
fi

#  Calculate the new MD% for this new file
echo "Calculating MD5 for $newrecab"
stime=`date +%s`
md5=(`md5sum $newrecab`)
b38cramchecksum=${md5[0]}
if [ "$b38cramchecksum" = "" ]; then
  Fail "MD5sum for recab failed"
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM for new recab completed in $etime seconds"

#  All is good, put this fie into production
mv $newrecab $recab
if [ "$?" != "0" ]; then
  echo "Unable to move remapped recab for $bamid build $b failed" >> $fixlog
  Fail "Unable to move remapped recab failed" 
fi
mv $newrecab.crai $recab.crai
if [ "$?" != "0" ]; then
  echo "Unable to move remapped crai for $bamid build $b failed" >> $fixlog
  Fail "Unable to move remapped craiß failed" 
fi

#  Update database
SetDB $bamid b38cramchecksum $b38cramchecksum
SetDB $bamid state_gce38copy 1
SetDB $bamid b38flagstat $newflagstat

echo "Corrected recab for $bamid successful" >> $fixlog
Successful
Log $etime
exit

#================= Old runs =================#
exit

. /usr/cluster/topmed/bin/topmed_actions.inc
me=fix
markverb=$me

if [ "$1" = "-submit" ]; then
  shift
  bamid=`$topmedcmd show $1 bamid`
  MayIRun $me  $bamid
  RandomRealHost $bamid
  SubmitJob $bamid "topmed-fix" '2G' "$0 $*"
  exit
fi

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid|nwdid"
  echo ""
  echo "Fix some sort of problem."
  exit 1
fi

Started
bamid=$1
b=b37
echo "# $0 $bamid $nwdid  -- set ${b}cramchecksum for $b remapped cram"

#   Get file to process
file=`$topmedpath wherefile $bamid $b`
if [ "$file" = "" ]; then
  echo "Unable to find file for '$b'"
  exit 4
fi 

#   Verify the MD5
stime=`date +%s`
md5=`md5sum < $file`
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM  completed in $etime seconds"
if [ "$rc" != "0" ]; then
  echo "MD5sum for $bamid build $b failed" >> $console/fix.log
  Fail "MD5sum failed" 
fi

SetDB $bamid ${b}cramchecksum $md5

echo "MD5sum for $bamid build $b successful" >> $console/fix.log
Successful
Log $etime
exit


nwdid=`GetNWDID $1`
bamid=`$topmedcmd show $1 bamid`

#   Get file to process
bcffile=`$topmedpath wherefile $bamid bcf`
if [ "$bcffile" = "" ]; then
  echo "Unable to find file for '$b'"
  exit 4
fi
if [ ! -f $bcffile ]; then
  echo "No BCF file found: $bcffile"
  exit 2
fi

vt=/usr/cluster/software/trusty/topmed-year1-freeze3a/master/vt/vt
vtref=/net/topmed/working/mapping/gotcloud/ref/hg38/hs38DH.fa
bcftools=/usr/cluster/bin/bcftools

$bcftools index $bcffile
if [ "$?" != "0" ]; then
  echo "Unable to run BCFTOOLS for bamid '$bamid' [$nwdid] on $bcffile"
  exit 4
fi
echo "Created index BCF file for $bcffile"
exit
#================= Old runs =================#
Started
bamid=$1
b=b37
echo "# $0 $bamid $nwdid  -- set ${b}cramchecksum for $b remapped cram"

#   Get file to process
file=`$topmedpath wherefile $bamid $b`
if [ "$file" = "" ]; then
  echo "Unable to find file for '$b'"
  exit 4
fi 

#   Verify the MD5
stime=`date +%s`
md5=`md5sum < $file`
rc=$?
etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM  completed in $etime seconds"
if [ "$rc" != "0" ]; then
  echo "MD5sum for $bamid build $b failed" >> $console/fix.log
  Fail "MD5sum failed" 
fi

SetDB $bamid ${b}cramchecksum $md5

echo "MD5sum for $bamid build $b successful" >> $console/fix.log
Successful
Log $etime
exit
