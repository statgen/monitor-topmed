#!/bin/bash
#
#   topmed_rgmmap.sh bamid|nwdid [markverb]
#
#	RG headers in remapped cram are wrong, rebuild fixed cram
#
. /usr/cluster/$PROJECT/bin/topmed_actions.inc
me=rgmap
markverb=''                 # Avoids setting database unless provided by

if [ "$1" = "" ]; then
  me=`basename $0`
  echo "Usage: $me bamid|nwdid  [markverb]"
  echo ""
  echo "RG headers in remapped cram are wrong, rebuild fixed cram."
  exit 1
fi

bamid=$1
if [ "$2" != "" ]; then
  markverb=$2               # Allow database state for an action to be set
fi
build=b38

#Started

echo "# $0 $bamid $nwdid  -- redo header for mis-mapped b38 crams"
recab=`$topmedpath wherefile $bamid b38`
tmpdir=/tmp

#  Create a corrected cram
b=`basename $recab`
newrecab=$tmpdir/$b
nwdid=`GetDB $bamid expt_sampleid`

#   Build correct rg mapping file. No check of return code as we already know it
rgmap=$tmpdir/$nwdid.rgmap
efile=/tmp/$bamid.err
/usr/cluster/topmed/bin/topmedrgmap.pl $bamid $rgmap 2> $efile
if [ "$?" != "0" ]; then
  e=`cat $efile`
  rm -f $efile
  Warn "Create rgmap failed: $e" 
  exit 2
fi
rm -f $efile

if [ ! -f $newrecab ]; then
  stime=`date +%s`
  /net/fantasia/home/hmkang/tools/apigenome.test/bin/cramore cram-update-rg --ref /net/mario/nodeDataMaster/local/ref/gotcloud.ref/hg38/hs38DH.fa --rg-map $rgmap --sam $recab --out $newrecab
  rc=$?
  rm -f $rgmap
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "Corrected recab completed in $etime seconds"
  if [ "$rc" != "0" ]; then
    Fail "Corrected recab failed" 
  fi
else
  echo "Using existing recab $newcab"
fi
#   Get flagstat values for new file, must match original flagstat
newflagstat=`CalcFlagstat $bamid $newrecab`
if [ "$newflagstat" = "" ]; then
  Fail "Unable to get reads of paired in sequencing for '$newrecab'"
fi

origflagstat=`GetDB $bamid cramflagstat`
if [ "$origflagstat" != "$newflagstat" ]; then
  Fail "New remapped file flagstat incorrect $origflagstat != $newflagstat"
fi
echo "New remapped file is valid, flagstats match ($newflagstat)"

#  Create a new index for the file
CreateIndex $bamid $newrecab

#  Calculate the new MD5 for this new file
echo "Calculating MD5 for $newrecab"
stime=`date +%s`
b38cramchecksum=`CalcMD5 $bamid $newrecab`
etime=`date +%s`
etime=`expr $etime - $stime`
echo "MD5SUM for new recab completed in $etime seconds"

#  All is good, put this file into production
echo "Replacing $recab with $newrecab"
mv $newrecab $recab
if [ "$?" != "0" ]; then
  Fail "Unable to move remapped recab failed" 
fi
mv $newrecab.crai $recab.crai
if [ "$?" != "0" ]; then
  Fail "Unable to move remapped craiß failed" 
fi

#  Update database
echo "Updating database values"
SetDB $bamid b38cramchecksum $b38cramchecksum
SetDB $bamid state_gce38copy 1
SetDB $bamid b38flagstat $newflagstat

echo "Corrected recab for $bamid $nwdid was successful"
Successful
exit
