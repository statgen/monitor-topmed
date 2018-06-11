#!/bin/bash
#
#   Fetch data from the BROAD from firecloud
#
#   Usage: cd where you want data copied
#       firecloud-fetch.sh URL
#   where URL was sent by Broad and probably looks like:
#        https://portal.firecloud.org/#workspaces/...
#
me=`basename $0`
#   Copy bucket data into here
ddndir=/net/ddn/incoming/topmed/broad
#   After we think DDN works, we download directly into DDN
tmpdir=$ddndir
tmpdir=/net/topmed9/incoming/topmed/year4-tempdata

if [ "$2" = "" ]; then
  echo "$me - Fetch data from the BROAD from firecloud"
  echo ""
  echo "Usage: $me newname url"
  echo "  newname - name of directory where DDN incoming data should go"
  echo "  url - firecloud URL provided by Broad"
  echo ""
  exit 1
fi
dir=$1
url=$2

echo -n "Do you really want to fetch bucket data to '$dir'? (y/n) "
read a
if [ "$a" != "y" ]; then
  echo "Nothing done"
  exit 3
fi

if [ -d $dir ]; then
  echo -n "Directory '$dir' already exists - continue? (y/n) "
  read a
  if [ "$a" != "y" ]; then
    echo "Nothing done"
    exit 3
  fi
else
  mkdir $dir || exit 4
  echo "Created directory '$dir'"
fi

#   Data goes into $ddndir/$dir   For now download to $tmpdir/$dir
cd $ddndir || exit 4
mkdir -p $dir
cd $tmpdir || exit 4            # Someday remove these two lines
mkdir -p $dir
cd $dir || exit 4

#   Get bucket name
token=`gcloud auth print-access-token`
if [ "$token" = "" ]; then
  echo "Unable to get gcloud token with 'gcloud auth print-access-token'"
  exit 2
fi
echo ""

#   We can't trust broad to send us exactly the same URL format
#   so here we try to construct the URL we really need
#   We want: https://api.firecloud.org/api/workspaces/broad-genomics-data/TOPMED_etc
#   but maybe we get: https://portal.firecloud.org/#workspaces/broad-genomics-data/TOPMED_etc
x=`echo $url | grep 'https://portal.firecloud.org/'`
if [ "$x" != "" ]; then
  echo "Converting URL into something that will work for us"
  x=`echo $url | tr -d \# | sed -e 's+https://portal.firecloud.org/++'`
  if [ "$x" = "" ]; then
    echo "Failed to construct URL like: https://api.firecloud.org/api/workspaces/broad-genomics-data/.."
    exit 3
  fi
  url="https://api.firecloud.org/api/$x"
  echo "We will use URL=$url"
fi

tmpf=firecloud.log.txt
date > $tmpf
echo $url >> $tmpf
echo "" >> $tmpf
curl -s -X GET -H "Authorization:Bearer $token " --header 'Accept: text/html' "$url" >> $tmpf
echo "" >> $tmpf

cat > $tmpf.pl <<EOF
\$a=\`cat $tmpf\`;
if (\$a =~ /.bucketName.:.([^\"]+)/) { print \$1; }
EOF

bucket=`perl $tmpf.pl`
if [ "$bucket" = "" ]; then
  echo "$me - Unable to find buckName from URL=$url" | tee -a $tmpf
  echo "See curl output in `pwd`/$tmpf" | tee -a $tmpf
  exit 2
fi
uri="gs://$bucket"
echo $uri >> $tmpf

#   Finally know what we are supposed to do
stime=`date +%s`
echo ""
echo "Fetching bucket $uri to `pwd`"
echo "This will take quite some time"
echo ""
gsutil -m rsync -r $uri .
if [ "$?" != "0" ]; then
  echo "$me -FAILED to fetch all the data from $uri" | tee -a $tmpf
  exit 6
fi
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Fetch of CRAMS completed in $etime seconds" | tee -a $tmpf

#   Must now get MD5 for each cram file
stime=`date +%s`
echo ""
echo "Now get the MD5s for CRAM files"
for cram in `find . -name \*.cram -print`; do
  c=`echo $cram | sed -e 's:./::'`
  gsutil hash -h $uri/$c | grep md5 > $c.md5
done
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Fetch of MD5s completed in $etime seconds" | tee -a $tmpf

#   Files are in holding directory - copy the stuff we want to DDN
#   Must also get MD5 for the cram file
#   Someday we copy files directly from the bucket to DDN.
stime=`date +%s`
echo ""
echo "Now copying files of interest into $ddndir/$dir"
for f in `find . -name \*.cram -print`; do
    nwd=${f%.*}
    rsync -a ${nwd}.* $ddndir/$dir 
done
etime=`date +%s`
etime=`expr $etime - $stime`
echo "Copy of CRAMS to DDN completed in $etime seconds" | tee -a $tmpf
cp -p $tmpf $ddndir/$dir

#   Clean out cruft from the Broad. For now the cruft is left in $tmpdir/$dir
cd $ddndir/$dir || exit 8
n=`ls *.cram|wc -l`
if [ "$n" = "0" ]; then
  echo "$me - something went wrong, no CRAMS found"
  exit 9
fi
echo ""
echo "SUCCESSFUL  Fetched $n $uri files to $ddndir/$dir" | tee -a $tmpf

echo "Creating Manifest.txt from the MD5 files" | tee -a $tmpf
mf=Manifest.txt
rm -f $mf
for f in *.cram; do
  md5=(`cat $f.md5`)
  echo "${md5[2]}   $f" >> $mf
done
echo "Created in $ddndir/$dir/$mf" | tee -a $tmpf
ls -l $mf
