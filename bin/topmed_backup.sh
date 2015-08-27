#!/bin/bash
#
#   topmed_backup.sh -submit bamid bamfile
#
#	Backup the original BAM as a simple copy
#
bindir=/usr/cluster/bin
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
backupdir=/net/topmed/working/backups
mem=${TOPMED_MEM:-1G}
console=/net/topmed/working/topmed-output
slurmp=${TOPMED_PARTITION:-topmed-incoming}
slurmqos=topmed-backup

if [ "$1" = "-submit" ]; then
  shift 
  homehost=`echo $2 | cut -d / -f 3`    # Should be topmed/topmed2
  if [ "$homehost" != "" ]; then
    console=/net/$homehost/working/topmed-output
    #slurmp="$homehost-incoming"
  fi
 
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-backup --output=$console/$1-backup.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-backup --output=$console/$1-backup.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 backedup submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobid
    $topmedcmd set $1 jobidbackup ${l[3]}
  fi
  exit
fi

if [ "$2" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid bamfile"
  echo ""
  echo "Backup a bam file and update database"
  exit 1
fi
bamid=$1
bamfile=$2

#   This is a bit of a hack to insure backups of data on topmed
#   go to topmed2 and vice versa
homehost=`echo $2 | cut -d / -f 3`      # Should be topmed/topmed2
if [ "$homehost" != "" ]; then
  if [ "$homehost" = "topmed" ]; then
    homehost=topmed2
  else
    if [ "$homehost" = "topmed2" ]; then
      homehost=topmed
    fi
  fi
  backupdir=/net/$homehost/working/backups
fi

#   Mark this as started  Calc dest directory in a poor manner
d=`echo $bamfile | sed -e s:/net/topmed/::`
d=`echo $d | sed -e s:/net/topmed2/::`
d=`echo $d | sed -e s:/net/topmed3/::`
dest=`dirname $d`
d=`date +%Y/%m/%d`
cd $backupdir
if [ "$?" != "0" ]; then
  echo "Unable to CD $backupdir"
  $topmedcmd mark $bamid backedup failed
  exit 2
fi
mkdir -p $dest
cd $dest
if [ "$?" != "0" ]; then
  echo "Unable to CD $backupdir/$dest"
  $topmedcmd mark $bamid backedup failed
  exit 2
fi
s=`hostname`
echo "#========= '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid bamfile=$bamfile dest=$backupdir/$dest ========="

#   Mark this as started
$topmedcmd mark $bamid backedup started
stime=`date +%s`
#   Try to force everything as readonly
chmod 555 $bamfile 2> /dev/null

#   Get NWDID
nwdid=`$bindir/samtools view -H $bamfile | grep '^@RG' | grep -o 'SM:\S*' | sort -u | cut -d \: -f 2`
if [ "$nwdid" = "" ]; then
  echo "Unable to extract NWDID from header of '$bamfile'"
  $topmedcmd mark $bamid backedup failed
  exit 2
fi

#   Backup using a simple copy
cp -p $bamfile .
rc=$?
if [ "$rc" != "0" ]; then
  echo "Copy failed: cp -p $bamfile $backupdir/$dest"
  $topmedcmd mark $bamid backedup failed
  exit 3
fi

etime=`date +%s`
etime=`expr $etime - $stime`
echo "Backup completed in $etime seconds to $backupdir/$dest"
$topmedcmd set $bamid nwdid $nwdid
if [ "$?" != "0" ]; then
  echo "Command failed: $topmedcmd set $bamid nwdid $nwdid"
  $topmedcmd mark $bamid backedup failed
  exit 3
fi

#   All was good
$topmedcmd mark $bamid backedup completed
