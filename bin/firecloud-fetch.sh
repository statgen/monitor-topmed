#!/bin/bash
#
#   Fetch data from the BROAD from firecloud
#
#   Usage:  firecloud-fetch.sh  URL  localdir
#
#   where URL was sent by Broad and probably looks like:
#         https://portal.firecloud.org/#workspaces/...
#   and   localdir  is either an intended directory name on 
#         /net/ddn/incoming/topmed/broad or else an absolute 
#         pathname to somewhere else
#   Not interactive, so this can be run in the background
#
me=`basename $0`
convert_checksums=/net/topmed/incoming/study.reference/study.reference/nhlbi.4275.gsutil.checksums.awk
export PROJECT=topmed

if [ "$#" -lt 2 ]; then
  echo "$me - Fetch data from the BROAD from firecloud"
  echo ""
  echo "Usage:  $me  url  localdir"
  echo "  url - firecloud URL provided by Broad"
  echo "  localdir - name of directory where DDN incoming data should go,"
  echo "  either a relative path in  /net/ddn/incoming/topmed/broad"
  echo "  or an absolute pathname for somewhere else"
  echo ""
  exit 1
fi

url=$1
localdir=$2
#   Copy bucket data into here
ddndir=/net/ddn/incoming/topmed/broad

abs=`echo $localdir | cut -c 1`
if [ "$abs" != "/" ]; then
  dir=${ddndir}/$localdir
else
  dir=$localdir
fi
mkdir -p $dir
chmod 0770 $dir || exit 5       # Be sure permissions are set properly
cd $dir || exit 4

logfile=fetch.log
echo "Log file at $dir/$logfile"
echo ''   >  $logfile
echo "`date`	`uname -n`" >> $logfile
echo $me  >> $logfile
echo $url >> $logfile
echo $localdir >> $logfile
echo '' >> $logfile

#   Get bucket name
token=`gcloud auth print-access-token`
if [ "$token" = "" ]; then
  echo "Unable to get gcloud token with 'gcloud auth print-access-token'" >> $logfile
  exit 2
fi

#   We can't trust broad to send us exactly the same URL format
#   so here we try to construct the URL we really need
#   We want: https://api.firecloud.org/api/workspaces/broad-genomics-data/Gabriel_TOPMED_etc
#   but maybe we get: https://portal.firecloud.org/#workspaces/broad-genomics-data/Gabriel_TOPMED_etc

x=`echo $url | grep 'https://portal.firecloud.org/'`
if [ "$x" != "" ]; then
  echo "Converting URL into something that will work for us" >> $logfile
  x=`echo $url | tr -d \# | sed -e 's+https://portal.firecloud.org/++'`
  if [ "$x" = "" ]; then
    echo "Failed to construct URL like: https://api.firecloud.org/api/workspaces/broad-genomics-data/.." >> $logfile
    exit 3
  fi
  url="https://api.firecloud.org/api/$x"
  echo "We will use URL=$url" >> $logfile
fi
echo '' >> $logfile
echo $url >> $logfile
echo '' >> $logfile
curl -s -X GET -H "Authorization:Bearer $token " --header 'Accept: text/html' "$url" >> $logfile
echo '' >> $logfile
echo '' >> $logfile

#   Build a perl script which picks the actual Google bucket name 
#   out of the reply from curl that's already in the logfile.

cat > extract.pl <<EOF
\$a=\`cat $logfile\`;
if (\$a =~ /.bucketName.:.([^\"]+)/) { print \$1; }
EOF

bucket=`perl extract.pl`
if [ "$bucket" = "" ]; then
  echo "$me - Unable to find bucketName from URL=$url" >> $logfile
  echo "See curl output in `pwd`/$logfile" >> $logfile
  exit 2
fi
rm -f extract.pl
uri="gs://$bucket"
echo $uri   >> $logfile
echo  ''    >> $logfile
echo "`date`\t \tBuild Manifest.txt" >> $logfile

#   Get Google data to build the file Manifest.txt
#   Checksums are space-separated, but 'cut -f' looks for tab.

gsutil ls -r -l ${uri}/'**' > gsutil.file.list
gsutil ls -r -L ${uri}/'**' > gsutil.complete
awk -f $convert_checksums  gsutil.complete > gsutil.checksums
grep -v cram.crai gsutil.checksums | cut -f 1 > Manifest.txt
echo `date` >> $logfile
rm -f gsutil.checksums

#   Finally know what we are supposed to do
mkdir incoming
cd    incoming || exit 5

stime=`date +%s`
echo '' >> ../$logfile
echo "Fetching bucket $uri to `pwd`"  >> ../$logfile
echo "This will take quite some time" >> ../$logfile
echo '' >> ../$logfile
gsutil -m -q rsync -r $uri . >> ../$logfile 2>&1
if [ "$?" != "0" ]; then
  echo "$me -FAILED to fetch all the data from $uri" >> ../$logfile
  exit 6
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Fetch of CRAMS completed in $etime seconds" >> ../$logfile

cd ..

#   Files are in holding directory 'incoming' - copy the stuff we want to $dir
#   This script already copies files directly from the bucket to DDN.
stime=`date +%s`
echo `date` >> $logfile
echo "Now copying files of interest into $dir" >> $logfile
for f in `find . -name \*.cram -print`; do
    nwd=${f%.*}
    mv ${nwd}.* .
done
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Copy of CRAMS to DDN completed in $etime seconds" >> $logfile
echo `date` >> $logfile
n=`ls *.cram|wc -l`
if [ "$n" = "0" ]; then
  echo "$me - something went wrong, no CRAMS found" >> $logfile
  exit 9
else
#   Clean out cruft from the Broad
  rm -rf incoming
fi
#   Insure the CRAI files are younger than the source file. Avoid re-making indexes. Correct permissions
touch *.crai
chmod 0660 *.cram *.crai

echo "SUCCESSFUL  Fetched $n $uri files to $ddndir/$dir"   >> $logfile
echo "Manifest.txt shows `ls -l Manifest.txt` .cram files" >> $logfile
echo `date` >> $logfile
echo  ''    >> $logfile
echo -n "$localdir	:  " >> $logfile
#   Register this new run in the monitor database
#   This part of the code is hard-wired with center = broad

linkdir=`echo $dir | sed 's:/net::'`
center=broad
cd /net/topmed/incoming/topmed/$center
ln -s ../../../../$linkdir .
cd $HOME
/usr/cluster/topmed/bin/topmed_init.pl -center $center updatedb >> $dir/$logfile
rc=$?
echo '' >> $dir/$logfile
exit $rc


