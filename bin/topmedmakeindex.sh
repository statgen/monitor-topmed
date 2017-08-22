#!/bin/bash
#
#   topmedmakeindex.sh bamorcramfile [build] [consolelog]
#
#	samtools can fail in multiple ways and never tell us
#   This script attempts to centralize all known attempts
#   to check for samtools errors for an index file
#   so we can possibly insure the index file will work  (POS)
#
#   This script should be called by the various shell scripts
#   which need a valid index
#
#   If the path to a consolelog is provided and the index file
#   fails, we can check for a broken samtools cache file
#
samtools=/usr/cluster/bin/samtools
topmedcmd=/usr/cluster/topmed/bin/topmedcmd.pl
b37="3:148,100,000-148,110,000"         # Search for this
b38="chr$b37"                           # Except this was changed for 38   arrgh

if [ "$1" = "" ]; then
  echo "Invalid call sequence:  $0 $*"
  exit 11
fi
file=$1
build="$2"
consolelog="$3"
buildstr=$b37
if [ "$build" = "38" -o "$build" = "b38" ]; then
  buildstr=$b38
fi

extension="${file##*.}"
if [ "$extension" = "bam" ]; then
  indexfile=$file.bai
elif [ "$extension" = "cram" ]; then
  indexfile=$file.crai
else
  echo "Unknown extension '$extension' for $file"
  exit 3
fi
basebam=`basename $file .$extension`

#   Check for a valid index file
if [ -z $indexfile ]; then
  rm -f $indexfile
fi
if [ -f $indexfile ]; then          # Index exists, let's see if it is any good
  n=`$samtools view $file $buildstr | wc -l`
  if [ "$n" = "" -o "$n" -lt "1000" ]; then
    echo "Index file '$indexfile' is invalid"
    #   This might be a trashed reference index for samtools, if so remove it
    if [ "$3" != "" ]; then
      a=`grep 'cram_ref_load: Assertion' $consolelog`
      if [ "$a" != "" ]; then
        rm -rf $HOME/.cache/hts-ref/*/*
        echo "Samtools cram_ref_load error, removed dirty reference cache"
      fi
    fi
    rm -f $indexfile
  else
    echo "Index file '$indexfile' is valid apparently"
    exit 0
  fi
fi

#   Create the index file
echo "Creating index file '$indexfile'"
$samtools index $file 2>&1
if [ "$?" = "0" ]; then
  echo "Created index file '$indexfile'"
  exit 0
fi
echo "Unable to create index file for '$file'"
exit 3

