#!/bin/bash
#
#   topmed_ncbi.sh -submit bamid + paris of cramfile colname
#   topmed_ncbi.sh -submit bamid XML  xml1  xml2 xml3
#
#	Send the proper set of files to NCBI
#
bindir=/usr/cluster/bin
ascpcmd="$bindir/ascp -i /net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri -l 800M -k 1"
#ascpcmd="$bindir/ascp -i /net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.ppk -Q -l 800M -k 1"
ascpdest='asp-um-sph@gap-upload.ncbi.nlm.nih.gov:protected'
topmedcmd=/usr/cluster/monitor/bin/topmedcmd.pl
mem=2G
console=/net/topmed/working/topmed-output

if [ "$1" = "-submit" ]; then
  shift
  #   May I submit this job?
  $topmedcmd permit test ncbi $1
  if [ "$?" = "0" ]; then
    exit 4
  fi 

  #   This will usually have files on both topmeds, so it does not matter where it runs
  l=(`$topmedcmd where $1`)     # Get bampath backuppath bamname realhost realhostindex
  realhost="${l[3]}"
  realhostindex="${l[4]}"
  slurmp="$realhost-incoming"
  slurmqos="$realhost-ncbi"

  #  Submit this script to be run
  l=(`/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-ncbi --output=$console/$1-ncbi.out $0 $*`)
  if [ "$?" != "0" ]; then
    echo "Failed to submit command to SLURM"
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp --mem=$mem --qos=$slurmqos --workdir=$console -J $1-ncbi --output=$console/$1-ncbi.out $0 $*"
    exit 1
  fi
  $topmedcmd mark $1 cped2ncbi submitted
  if [ "${l[0]}" = "Submitted" ]; then      # Job was submitted, save jobid
    $topmedcmd set $1 jobidcp2ncbi ${l[3]}
  fi
  exit
fi

if [ "$3" = "" ]; then
  me=`basename $0`
  echo "Usage: $me [-submit] bamid cramfile colname [cramfile colname]"
  echo "Usage: $me [-submit] bamid XML xml1 xml2 xml3"
  echo ""
  echo "Copy CRAM files to NCBI"
  exit 1
fi
bamid=$1
shift

d=`date +%Y/%m/%d`
echo "#========= $d $SLURM_JOB_ID $0 bamid=$bamid files=$* ========="

#   Well hidden documentation for Aspera command line syntax comes from:
#       http://download.asperasoft.com/download/docs/ascp/3.5.2/html/
#       http://download.asperasoft.com/download/docs/ascp/2.7/html/
#       http://download.asperasoft.com/download/docs/ascp/2.7/html/fasp/ascp-usage.html.

#   Special case - if bamid = XML, then send all XML files at once
if [ "$1" = "XML" ]; then
  shift
  echo "Sending XML files to NCBI - $*"
  $ascpcmd $* $ascpdest                 # We count on this working, no retries
  if [ "$?" = "0" ]; then
    echo "XML files sent to NCBI"       # Nothing to mark in the database as success
    exit
  fi
  echo "FAILED sending XML files to NCBI"       # We have no real bamid for the XML files
  $topmedcmd mark $bamid cped2ncbi failed       # Mark first bamid we sent as failed
  exit 1
fi

#   Here we send one file at a time (usually just two files)
#   Columns must come in pairs - path to file and database column to set to Y
while [ "$1" != "" ]; do
  f=$1
  dbcol=$2
  shift
  shift
  if [ "$dbcol" = "" ]; then
    echo "Not enough parameters: $0 $bamid $*"
    exit 2
  fi

  echo "Sending '$f' to NCBI  bamid=$bamid"
  stime=`date +%s`

  # We really really want this to complete, cause recovery is a giant pain
  retries=10
  waitsec=60
  until [ "$retries" = "0" ]; do
    $ascpcmd $f $ascpdest
    rc=$?
    if [ "$rc" = "0" ]; then
      $topmedcmd set $bamid $dbcol Y        # Mark this CRAM as delivered
      retries=0
    else 
      echo "FAILED sending '$f' to NCBI, retrying in $waitsecs seconds"
      retries=`expr $retries - 1`
      sleep $waitsec
      waitsec=`expr $waitsec + 60`
    fi
  done
  if [ "$rc" != "0" ]; then
    echo "I've tried about all I can try"
    $topmedcmd mark $bamid cped2ncbi failed
    exit $rc
  fi
 
  etime=`date +%s`
  etime=`expr $etime - $stime`
  echo "Sent '$f' in $etime seconds for bamid=$bamid"
done

#   All files for this bamid have been delivered
$topmedcmd mark $bamid cped2ncbi delivered
exit

#   Somehow someday this we also need to do
$topmedcmd mark $bamid cped2ncbi completed
exit 2

