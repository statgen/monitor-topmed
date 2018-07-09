#!/bin/bash 
####################################################################
#
# backupmysql.sh
#
#   Use this script to backup a MySQL database and put the result
#   in a specified directory
#   PROJECT must be set
#
#   Input:
#    realm    realm to use
#    destdir  directory
#    name     optional name of backup file
#
#   E.g.  dsc-backupmysql.sh archives $HOME/MYSQLBAK
#
####################################################################
me=`basename $0`
realm=$LOGNAME                          # Also project
if [ "$HOME" != "" ]; then              # Sometimes cron does not set HOME
  HOME=/net/topmed/working/home/topmed
fi
realmdir=$HOME/.db_connections
pgm=/usr/bin/mysqldump
destdir=$HOME/BACKUPS/MYSQLBAK
name=""

#	Allow realm and database to be overridden
if [ $# -ge 1 ]; then
  realm=$1
fi
if [ $# -ge 2 ]; then
  destdir=$2
fi
if [ $# -ge 3 ]; then
  name=$3
fi

#   Figure out connection details
dbconnect="$realmdir/$realm"
u=`grep USER= $dbconnect | sed -e 's/USER=//'`
p=`grep PASS= $dbconnect | sed -e 's/PASS=//'`
h=`grep SERVER=host= $dbconnect | sed -e 's/SERVER=host=//'`
db=`grep DATABASE= $dbconnect | sed -e 's/DATABASE=//'`
pgmopts="--opt --host=$h --user=$u --password=$p --single-transaction"

#	Do actual backup
if [ "$name" = "" ]; then
  name=`date +%d`
fi
f=$destdir/$realm-$name.$db.sql.gz

echo "Backup of database $db ..."
$pgm $pgmopts $db | gzip -c > $f
if [ $? != 0 ]; then
  echo "$me - Error backing up database '$db' (realm=$realm)"
  ls -la $f
  exit 1
fi
chmod 600 $f
ls -la $f
