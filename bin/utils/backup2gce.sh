#!/bin/bash 
####################################################################
#
# backup2gce.sh
#
#   Use this script to backup software and meta-data to GCE
#   PROJECT must be set
#
#   Input:
#    dir  directory list to copy to GCE
#
####################################################################
me=`basename $0`
p=${PROJECT-noproject}
gceuri=gs://topmed-irc-working/software-backups/MYSQLBAK
if [ "$1" = "" ]; then
  echo "$me dir-to-backup-to-GCE"
  echo ""
  echo "You must set the environment variable PROJECT"
  exit 1
fi

#	Do actual copy
gsutil cp -r $1/* $gceuri

