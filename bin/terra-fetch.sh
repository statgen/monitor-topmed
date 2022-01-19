#!/bin/bash
#
#   Fetch data from the BROAD using Terra (new Firecloud)
#
#   Usage:  terra-fetch.sh [-v]  URI  localdir
#
#   The URL sent by Broad probably looks like:
#         https://app.terra.bio/#workspaces/broad-genomics-delivery/...
#         Using a browser to visit this URL, find Files where the Google bucket 
#         URI can be seen.  Copy the bucket URI into the command line above.
#         Should look like  gs://fc-c41c35f9-a9a0-4f34-9783-31be6b07487e
#         and it can include directory levels inside the bucket.
#   localdir  is either an intended directory name on 
#         /net/ddn/incoming/topmed/broad or else an absolute 
#         pathname to somewhere else
#   This is not interactive, so this can be run in the background
#

me=`basename $0`
convert_checksums=/net/topmed/incoming/study.reference/study.reference/nhlbi.4275.gsutil.checksums.awk
export PROJECT=topmed
ddndir=/net/ddn/incoming/topmed/broad       #  Copy bucket data into here

#------------------------------------------------------------------
# Subroutine:  LogMsg msg
#------------------------------------------------------------------
function LogMsg {
  local msg="$1"
  d=`date`
  if [ "$verbose" = "1" ]; then echo "$msg"; fi
  echo "$d  $msg" >> $logfile
}

if [ "$#" -lt 2 ]; then
  echo "$me - Fetch data from the BROAD using Terra"
  echo ""
  echo "Usage:  $me  [-v] uri  localdir"
  echo "  -v  - verbose, tell me more as steps are completed"
  echo "  uri - Google storage URI (gs://fc...) needed to access the data tree"
  echo "  localdir - name of directory where DDN incoming data should go,"
  echo "    either a relative path in  /net/ddn/incoming/topmed/broad"
  echo "    or an absolute pathname for somewhere else"
  echo ""
  exit 1
fi
verbose=0
if [ "$1" = "-v" ]; then
  verbose=1
  shift
fi
uri=$1
localdir=$2

#   Remove trailing / in uri so we have the correct GS URI
i=$((${#uri}-1))
c=`echo "${uri:$i:1}"`
if [ "$c" = "/" ]; then
  i=${#uri}
  uri=`echo ${uri::$i-1}`
fi

#   Test whether localdir is an absolute path name
abs=`echo $localdir | cut -c 1`
if [ "$abs" != "/" ]; then
  dir=${ddndir}/$localdir
else
  dir=$localdir
fi

mkdir -p $dir
chmod 0770 $dir || exit 5       # Be sure permissions are set properly
cd $dir || exit 4

logfile=$dir/fetch.log	 	# Use an absolute path name for logfile
echo "Log file at $logfile"
LogMsg "$me"
LogMsg $localdir
LogMsg "Using Google Store $uri"

#   As of February 2021, Google evidently no longer guarantees to provide md5 
#   checksums as part of its "-L" long file listing.  So instead we must use 
#   the files "NWD*.cram.md5" provided by Broad to construct "Manifest.txt" 
#   at the end, after all of the files have been downloaded.

#   Old version -- Get Google data to build the file Manifest.txt
LogMsg "Retrieve Google file lists"
gsutil ls -r -l ${uri}/'**' > gsutil.file.list
gsutil ls -r -L ${uri}/'**' > gsutil.complete
awk -f $convert_checksums  gsutil.complete > gsutil.checksums
#   grep -v cram.crai gsutil.checksums | grep -v cram.md5 | cut -f 1 > Manifest.txt
#   Checksums are space-separated, but 'cut -f' looks for tab.
echo `date` >> $logfile
LogMsg "File lists retrieved"

#   Finally know what we are supposed to do
mkdir incoming
cd    incoming || exit 5

#	Set the Boto values to meet your needs
stime=`date +%s`
opts="-o Boto:parallel_process_count=10 -o Boto:parallel_thread_count=60"
LogMsg "Fetching bucket $uri to `pwd`"
echo "This will take quite some time. Using $opts"
gsutil -m $opts -q rsync -r $uri . >> $logfile 2>&1
if [ "$?" != "0" ]; then
  LogMsg "$me -FAILED to fetch all the data from $uri - See $logfile"
  exit 6
fi
etime=`date +%s`
etime=`expr $etime - $stime`
LogMsg "Fetch of CRAMS completed in $etime seconds"
echo "See files in `pwd` and $logfile"

#   Files are in holding directory 'incoming' - copy the stuff we want to $dir
#   Using 'find' below makes this independent of the Broad directory structure.
#   Usually, this script has already copied files directly from the bucket into DDN.

cd ..
stime=`date +%s`
LogMsg "Now copying files of interest into $dir"
for f in `find . -name \*.cram -print`; do
    nwd=${f%.*}
    mv ${nwd}.* .
done
for file in `ls NWD*.md5`; do
    echo `cat $file` ' '`echo $file | sed -e 's/.md5//'` >> Manifest.txt
done
etime=`date +%s`
etime=`expr $etime - $stime`
LogMsg "Copy of CRAMS and Manifest.txt completed in $etime seconds"
#   Partial check of results:
numa=`cat Manifest.txt | wc -l`
numc=`ls *.cram | wc -l`
if [ "$numc" != "$numa" -o "$numc" = "0" ]; then
  LogMsg "$me - something went wrong, no CRAMS found or Manifest incomplete"
  exit 9
else
#   Clean out cruft from the Broad (Terry's former action)
#   rm -rf incoming	#  leave the empty Broad directory tree.
ls . | wc	 	#  have to DO something inside the "else"
fi
#   Insure the CRAI files are younger than the source file. Avoid re-making indexes. Correct permissions
touch *.crai
chmod 0640 *.cram *.crai

LogMsg "SUCCESSFUL  Fetched $numc $uri files to $dir"
LogMsg "Manifest.txt shows  $numa  .cram files"

#   Register this new run in the monitor database
#   This part of the code is hard-wired with center = broad
linkdir=`echo $dir | sed 's:/net::'`
center=broad
cd /net/topmed/incoming/topmed/$center
ln -s ../../../../$linkdir .
cd $HOME
/usr/cluster/topmed/bin/topmed_init.pl -center $center updatedb >> $logfile
rc=$?
echo '' >> $logfile
exit $rc

