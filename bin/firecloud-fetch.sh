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
tmpdir=bucket-incoming

if [ "$2" = "" ]; then
  echo "$me - Fetch data from the BROAD from firecloud"
  echo ""
  echo "Usage: $me newname url"
  echo "  newname - name of directory where incoming data should go"
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

#   Data goes into $dir eventually, but we really download to incoming initially
mkdir $dir/$tmpdir 2> /dev/null
cd $dir || exit 4
cd $tmpdir || exit 4

#   Get bucket name
token=`gcloud auth print-access-token`
if [ "$token" = "" ]; then
  echo "Unable to get gcloud token with 'gcloud auth print-access-token'"
  exit 2
fi
echo ""

tmpf=curldata
curl -s -X GET -H "Authorization:Bearer $token " --header 'Accept: text/html' "$url" > $tmpf

cat > $tmpf.pl <<EOF
\$a=\`cat $tmpf\`;
if (\$a =~ /.bucketName.:.([^\"]+)/) { print \$1; }
EOF

bucket=`perl $tmpf.pl`
if [ "$uri" = "" ]; then
  echo "$me - Unable to find buckName from the URL you provided"
  echo "See curl output in $tmpf"
  exit 2
fi
uri="gs://$uri"

#   Finally know what we are supposed to do
echo "Fetching data from bucket $uri"
gsutil -m rsync $uri .
if [ "$?" != "0" ]; then
  echo "$me -FAILED to fetch all the data from $uri"
  exit 6
fi

echo "Now merging files of interest into $dir"
cd $dir || exit 4
for x in cram crai md5; do
  echo -n "$x "
  find $tmpdir -name \*.$x -exec mv {} . \;
done
echo ""
n=`ls *.cram|wc -l`
if [ "$n" = "0" ]; then
  echo "$me - something went wrong, no CRAMS found"
  exit 9
fi
echo "Cleaning out detritus from the bucket"
rm -rf $tmpdir
echo "SUCCESSFUL  Fetched $n $uri files to $dir"

echo "Creating Manifest.txt from the MD5 files"
mf=Manifest.txt
rm -f $mf
for f in *.cram; do
  md5=`cat $f.md5`
  echo "$md5   $f" >> $mf
done
ls -l $mf

  

