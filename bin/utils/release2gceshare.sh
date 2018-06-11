#!/bin/bash
#
#  release2gceshare.sh go
#     Given a list of sampleids (NWDID), update the database so these are copied to GCE
#
topmedcmd=/usr/cluster/topmed/bin/topmedcmd.pl
me=`basename $0`

if [ "$1" != "go" ]; then
  echo "$me go"
  echo ""
  echo "Set the monitor database flag so a sample will be copied to the GCE requestor-pays bucket"
  echo "You must pipe a list of sampleids to this script"
  echo ""
  echo "E.g.  cat my-nwdid-list.txt | $0 go"
  exit 1
fi

tmpf=/tmp/$me.tmp
cat > $tmpf
echo "Catpured your list of sampleids"
n=`head -1 $tmpf|grep NWD`
if [ "$n" = "" ]; then
  echo "You list of NWDIDs did have NWD"
  head -3 $tmpf
fi

k=0
for n in `cat $tmpf`; do
  $topmedcmd set $n state_gce38copy 1
  if [ "$?" = "0" ]; then
    k=`expr $k + 1`
  fi
done
echo "Set gce38copy flag for $k samples" 
rm -f $tmpf
