#!/bin/bash

#
# ascp download from shares
#
# Broad Institute, BITS
#
# Fragment borrowed from ascp-shares.sh Authored by Aspera
#
#   This is a modified version by Tom and is a version of the script
#   used by Tom to download broad files manually.  June 2016
#
#   This was modified to provide hooks so we can fetch the Manifest.txt
#   file which tells us what files can be downloaded.
#   This version will write key information to /tmp so we can extract
#   the raw command that needs to be used for all the other files to fetch.
#   It is used in conjunction with topmed_aspera.pl
#

DBG=false
#DBG=true
HOOK=false
LOGDIR=''

#
# Internals
#

_usage() {
  echo " To decrypt files, please set the  ASPERA_SCP_FILEPASS environemnt variable."
  echo
  echo "  Usage: $0 [-hook] [-L dir] SHARES_URL  USER:PASSWORD  /SHARES/PATH/TO/SOURCE/FILE  LOCAL_TARGET"
  echo
  echo "  Example: ascp-shares-down.sh https://shares.broadinstitute.org SN0000000:ABCD1234efgh /SN00000000/foo.bar /home/user/me/my_aspera_downloads"
  echo
  echo "  This version strips source pathnames for downloading multiple files within a directory."
  echo
  exit 1
}

_prep() {
  MKTEMP=mktemp
}
  
_get_response_from_json() {
  JSON=\{\"transfer_requests\":\[\{\"transfer_request\":\{\"paths\":\[\{\"source\":\"$SHARES_SOURCE\"\}\]\}\}\]\}
  JSONFILE=`$MKTEMP`
  echo $JSON > $JSONFILE
  
  if ($DBG) 
  then
    echo "--------------------------"
    echo
    echo JSON:
    echo
    echo `cat $JSONFILE`
    echo
  fi
  
  RESPONSE=`curl -ki -sS -H "Accept: application/json" -H "Content-type: application/json" -u $USERPASS -X POST -d @$JSONFILE $URL_DOWNLOAD_SETUP`
  rm $JSONFILE

  if ($DBG)
  then
    echo "--------------------------"
    echo
    echo RESPONSE: 
    echo
    echo "$RESPONSE"
    echo
  fi

}

_download_setup() {
_get_response_from_json


  TOKEN=`echo $RESPONSE | sed 's/\\\\\//\//g' | sed 's/[{}]//g' | awk -v k="text" '{n=split($0,a,","); for (i=1; i<=n; i++) print a[i]}'   | sed 's/\"//g' | grep -w token | cut -d: -f2`
  SOURCE=`echo $RESPONSE | sed 's/\\\\\//\//g' | sed 's/[{}]//g' | awk -v k="text" '{n=split($0,a,","); for (i=1; i<=n; i++) print a[i]}'   | sed 's/\"//g' | grep -w source | cut -d: -f8 | sed 's/\]//g'`
  USER=`echo $RESPONSE | sed 's/\\\\\//\//g' | sed 's/[{}]//g' | awk -v k="text" '{n=split($0,a,","); for (i=1; i<=n; i++) print a[i]}'   | sed 's/\"//g' | grep -w remote_user | cut -d: -f2`
  HOST=`echo $RESPONSE | sed 's/\\\\\//\//g' | sed 's/[{}]//g' | awk -v k="text" '{n=split($0,a,","); for (i=1; i<=n; i++) print a[i]}'   | sed 's/\"//g' | grep -w remote_host | cut -d: -f2`

  if ($DBG)
  then
    echo "--------------------------"
    echo
    echo VALUES FROM RESPONSE: 
    echo TOKEN: $TOKEN
    echo SOURCE: $SOURCE
    echo USER: $USER
    echo HOST: $HOST
    echo
  fi


}

_temp_key_make() {
  KEY=`$MKTEMP`
    echo "-----BEGIN DSA PRIVATE KEY-----
MIIBuwIBAAKBgQDkKQHD6m4yIxgjsey6Pny46acZXERsJHy54p/BqXIyYkVOAkEp
KgvT3qTTNmykWWw4ovOP1+Di1c/2FpYcllcTphkWcS8lA7j012mUEecXavXjPPG0
i3t5vtB8xLy33kQ3e9v9/Lwh0xcRfua0d5UfFwopBIAXvJAr3B6raps8+QIVALws
yeqsx3EolCaCVXJf+61ceJppAoGAPoPtEP4yzHG2XtcxCfXab4u9zE6wPz4ePJt0
UTn3fUvnQmJT7i0KVCRr3g2H2OZMWF12y0jUq8QBuZ2so3CHee7W1VmAdbN7Fxc+
cyV9nE6zURqAaPyt2bE+rgM1pP6LQUYxgD3xKdv1ZG+kDIDEf6U3onjcKbmA6ckx
T6GavoACgYEAobapDv5p2foH+cG5K07sIFD9r0RD7uKJnlqjYAXzFc8U76wXKgu6
WXup2ac0Co+RnZp7Hsa9G+E+iJ6poI9pOR08XTdPly4yDULNST4PwlfrbSFT9FVh
zkWfpOvAUc8fkQAhZqv/PE6VhFQ8w03Z8GpqXx7b3NvBR+EfIx368KoCFEyfl0vH
Ta7g6mGwIMXrdTQQ8fZs
-----END DSA PRIVATE KEY-----
" > $KEY
    }

_download() {
  
_temp_key_make

  #    Add exclusion of single character directories if we are hooked
  EXCLUDE=''
  if ($HOOK); then
    EXCLUDE="-E '?'"
  fi

  if [ -z ${ASPERA_SCP_FILEPASS}"" ]; then
    echo "No ASPERA_SCP_FILEPASS environment variable found. Downloading without decryption"
    DECRYPT=''
  else
    DECRYPT='--file-crypt=decrypt'
  fi

  #   Original commands were
  #   CMD="ascp --ignore-host-key -k 1 -Q -l 1500M -d -i $KEY -W $TOKEN --src-base=$SOURCE -L $LOGDIR  $USER@$HOST:$SOURCE $TARGET >> $RECORD  2>&1"
  #      or
  #   CMD="ascp --ignore-host-key -k 1 -Q -l 1500M -d --file-crypt=decrypt -i $KEY -W $TOKEN --src-base=$SOURCE -L $LOGDIR  $USER@$HOST:$SOURCE $TARGET >> $RECORD  2>&1"
  BASECMD="ascp --ignore-host-key -k 1 -Q -l 1500M -d -i $KEY -W $TOKEN $DECRYPT"
  #   Capture important information so we can fetch individual files
  if ($HOOK); then
    b=`basename $0`
    hook=/tmp/$b.hooks
    cp -p $KEY $hook.key
    echo "keyfile=$hook.key" > $hook
    echo "token=$TOKEN" >> $hook
    echo "decryptpasswd=$ASPERA_SCP_FILEPASS" >> $hook
    echo "basecmd=$BASECMD" >> $hook
    echo "Wrote hooks to '$hook'"  
  fi
  CMD="$BASECMD $EXCLUDE $LOGDIR --src-base=$SOURCE $USER@$HOST:$SOURCE $TARGET 2>&1"

  if ($DBG); then
    echo $CMD
  fi
  eval $CMD
  rm $KEY
}

#
# Main program
#

if [ "$1" = "-hook" ]; then
  HOOK=true
  shift
fi
if [ "$1" = "-L" ]; then
  LOGDIR="-L $2"
  shift
  shift
fi

if test $# -lt 3; then
  _usage
fi

URL_DOWNLOAD_SETUP=$1/node_api/files/download_setup
USERPASS=$2
SHARES_SOURCE=$3

if test $# -lt 4; then
   TARGET=`pwd`
else
   TARGET=$4
fi

_prep
_download_setup
_download