#
#	Define variables that all TOPMED action scripts need defined
#
topmedcmd=/usr/cluster/topmed/bin/topmedcmd.pl
topmedpermit=/usr/cluster/topmed/bin/topmedpermit.pl
topmedpath=/usr/cluster/topmed/bin/topmedpath.pl
topmedrename=/usr/cluster/topmed/bin/topmedrename.pl
topmedqplot=/usr/cluster/topmed/bin/topmedqplot.pl
bam=/usr/cluster/bin/bam
samtools=/usr/cluster/bin/samtools
vt=/usr/cluster/software/trusty/topmed-year1-freeze3a/master/vt/vt
vtref=/net/topmed/working/mapping/gotcloud/ref/hg38/hs38DH.fa
bcftools=/usr/cluster/bin/bcftools
console=/net/topmed/working/topmed-output
slurmp=topmed-working
build=38
realhost=''

cores=''                            # Set to same as gsutil if it will be used
gsutil="gsutil -m"
gsutilbig='gsutil -o GSUtil:parallel_composite_upload_threshold=150M'
#   No need to parallel processing 
#cores="--cpus-per-task=2"         # Cores here should be same as gsutil
#gsutil='gsutil -o GSUtil:parallel_composite_upload_threshold=150M -o GSUtil:parallel_process_count=2'
incominguri='gs://topmed-recabs'
bcfuri='gs://topmed-bcf'
export BOTO_CONFIG=/net/topmed/working/shared/tpg_gsutil_config.txt

#==================================================================
#   Define useful functions to be used by the action scripts
#==================================================================

#------------------------------------------------------------------
# Subroutine:
#   MayIRun(name, bamid)    Is bamid 'bamid' permitted to run for script 'name'
#
#   Does not return if permission is denied
#------------------------------------------------------------------
function MayIRun {
  local name=$1
  local bamid=$2
  $topmedpermit permit test $name $bamid
  if [ "$?" = "0" ]; then
    echo "$name $bamid not permitted" | tee $console/$bamid-$name.out
    exit 4
  fi 
}

#------------------------------------------------------------------
# Subroutine:
#   Log(timel)
#   Uses global variables: $SLURM_JOB_ID $console $bamid $me
#------------------------------------------------------------------
function Log {
  local timel=$1
  echo `date` $me $SLURM_JOB_ID ok $timel secs >> $console/$bamid.jobids
}

#------------------------------------------------------------------
# Subroutine:
#   Fail(msg)             Mark database as failed, save error msg
#   Uses global variables: $bamid $markverb
#   Does not return
#------------------------------------------------------------------
function Fail {
  local msg="$1"
  if [ "$markverb" != "" ]; then
    $topmedcmd -persist -emsg "$msg" mark $bamid $markverb failed
  else
    echo $msg
    echo "Task for $bamid FAILED, no markverb set"
  fi
  exit 1
}

#------------------------------------------------------------------
# Subroutine:
#   Started()               Mark database as started
#   Uses global variables: $bamid $markverb
#------------------------------------------------------------------
function Started {
  if [ "$markverb" != "" ]; then
    $topmedcmd -persist mark $bamid $markverb started
  fi
  d=`date +%Y/%m/%d`
  s=`hostname`
  echo "#==Started== '$d' host=$s $SLURM_JOB_ID $0 bamid=$bamid ========="
}

#------------------------------------------------------------------
# Subroutine:
#   Successful()            Mark database as completed
#   Uses global variables: $bamid $markverb
#------------------------------------------------------------------
function Successful {
  if [ "$markverb" != "" ]; then
    $topmedcmd -persist -emsg "" mark $bamid $markverb completed
  else
    echo "Task for $bamid successful, no markverb set"
  fi
}

#------------------------------------------------------------------
# Subroutine:
#   SetDB(bamid, col, val)
#------------------------------------------------------------------
function SetDB {
  local bamid=$1
  local col=$2
  local val=$3
  $topmedcmd -persist set $bamid $col "$val"
  if [ "$?" != "0" ]; then
    Fail "Unable to set column '$col' to '$val' for bamid '$bamid'" 
  fi
}

#------------------------------------------------------------------
# Subroutine:
#   val=`GetDB(bamid, col)`
#   Uses global variables: $bamid
#------------------------------------------------------------------
function GetDB {
  local bamid=$1
  local col=$2
  $topmedcmd show $bamid $col
  if [ "$?" != "0" ]; then
    Fail "Unable to get column '$col' for bamid '$bamid'" 
  fi
}

#------------------------------------------------------------------
# Subroutine:
#   GetNWDID(bamid)         prints $nwdid
#------------------------------------------------------------------
function GetNWDID {
  local bamid=$1
  nwdid=`$topmedcmd -persist show $bamid expt_sampleid`
  if [ "$nwdid" = "" ]; then
    Fail "Unable to find the NWDID for '$bamid'" 
  fi
  echo $nwdid
}

#------------------------------------------------------------------
# Subroutine:
#   MyRealHost(bamid, type)   Set host where type of input file lives
#   Sets $realhost  (--nodelist=hostname)
#------------------------------------------------------------------
function MyRealHost {
  local bamid=$1
  local type=$2
  realhost=`$topmedpath whathost $bamid $type`
}

#------------------------------------------------------------------
# Subroutine:
#   RandomRealHost(bamid)   Set host based on bamid
#   Sets $realhost
#------------------------------------------------------------------
function RandomRealHost {
  local bamid=$1

  realhost=''
  # All hosts are equal, so try to spread it around
  n=`perl -e "{ print $bamid%10 }"`
  h="topmed$n"
  if [ "$h" = "topmed0" ]; then
    h=topmed10
  fi
  if [ "$h" = "topmed8" ]; then     # May not use this
    h=topmed5
  fi
  if [ "$h" = "topmed1" ]; then
    h=topmed
  fi
  realhost=$h
}

#------------------------------------------------------------------
# Subroutine:
#   CheckRGMap(bamid)     rc=0 if rgmap not needed, rc=1 if needed
#------------------------------------------------------------------
function CheckRGMap {
  local bamid=$1
  
  /usr/cluster/topmed/bin/topmedrgmap.pl -check $bamid /dev/null
  return $?
}


#------------------------------------------------------------------
# Subroutine:
#   CreateIndex(bamid,file)    Create CRAM or BAM index
#
#	samtools can fail in multiple ways and never tell us
#   This script attempts to centralize all known attempts
#   to check for samtools errors for an index file
#   so we can possibly insure the index file will work  (POS)
#
#   Attempts to use global variable $me
#------------------------------------------------------------------
function CreateIndex {
  local bamid=$1
  local file=$2

  #  To verify an index, we search for a string - depending on the extension
  ext="${file##*.}"
  buildstr='notset'
  if [ "$ext" = "bam" ]; then
    buildstr="3:148,100,000-148,110,000"
  fi
  if [ "$ext" = "cram" ]; then
    buildstr="chr3:148,100,000-148,110,000"
  fi
  if [ "buildstr" = 'notset' ]; then
    Fail "Unable to determine build for $file"
  fi

  extension="${file##*.}"
  indexfile='notset'
  if [ "$extension" = "bam" ]; then
    indexfile=$file.bai
  fi
  if [ "$extension" = "cram" ]; then
    indexfile=$file.crai
  fi
  if [ "$indexfile" = "notset" ]; then
    Fail "Unknown extension '$extension' for $file"
  fi

  #   If index file is zero length, recreate it
  if [ -z $indexfile ]; then
    rm -f $indexfile
  fi

  #   Index exists, let's see if it is any good
  if [ -f $indexfile ]; then
    n=`$samtools view $file $buildstr | wc -l`
    if [ "$n" = "" -o "$n" -lt "1000" ]; then
      echo "Index file '$indexfile' is invalid"
      rm -f $indexfile                  # Remove bad index
      #   This might be a trashed reference index for samtools, if so remove it
      consfile="$console/$bamid-$me.out"
      a=`grep 'cram_ref_load: Assertion' $consfile`     # Could fail if $me not defined
      if [ "$a" != "" ]; then
        rm -rf $HOME/.cache/hts-ref/*/*
        echo "Samtools cram_ref_load error, removed dirty reference cache and retrying"
      fi
    else
      echo "Index file '$indexfile' is valid apparently"
      return
    fi
  fi

  #   Create the index file
  echo "Creating index file '$indexfile'"
  $samtools index $file 2>&1
  if [ "$?" = "0" ]; then
    echo "Created index file '$indexfile'"
  else
    Fail "Unable to create index file for '$file'"
  fi
}


#------------------------------------------------------------------
# Subroutine:
#   CalcFlagstat(bamid,file)    Prints flagstat for a bam or cram
#
#	Calculate the flagstat paired in sequencing number for a bam or cram
#------------------------------------------------------------------
function CalcFlagstat {
  local bamid=$1
  local file=$2

  # Special hack for this center
  ref=/net/mario/gotcloud/gotcloud.ref/hs37d5.fa
  c=`$topmedcmd show $bamid center`
  if [ "$c" = "illumina" ]; then
    ref=/net/topmed/incoming/study.reference/study.reference/illumina.hg19.fa
  fi

  #   Get flagstat values for input file
  of=`basename $file`               # Calculate temp file for flagstat results
  of="/tmp/$of.tmp"
  >&2 echo "Calculating flagstat for $file"
  $samtools flagstat  $file > $of
  if [ "$?" != "0" ]; then
    Fail "$samtools flagstat $file failed. Results in $of"
  fi

  #   Outfile has results of flagstat call, get paired reads
  a=(`grep 'paired in sequencing' $of`)
  if [ "${a[0]}" = "" ]; then
    Fail "Unable to get reads of paired in sequencing from '$of' for '$file'"
  fi
  rm -f $outfile
  echo ${a[0]}
}


#------------------------------------------------------------------
# Subroutine:
#   CalcMD5(bamid, file)    Prints md5 for a file
#------------------------------------------------------------------
function CalcMD5 {
  local bamid=$1
  local file=$2

  md5=`md5sum $file`
  checksum=${md5[0]}
  if [ "$checksum" = "" ]; then
    Fail "Unable to calculate MD5 for '$file' [$bamid]"
  fi
  echo $checksum
}

#------------------------------------------------------------------
# Subroutine:
#   SubmitJob(bamid, qos, mem, cmdtorun)    Submit this script to SLURM
#   Uses global variables: $me $markverb $console $slurp $realhost $cores
#------------------------------------------------------------------
function SubmitJob {
  local bamid=$1
  local qos=$2
  local mem=$3
  local cmdtorun="$4"
  if [ "$qos" != "" ]; then        # Use value passed or caller must set it
    q="--qos=$qos"
  fi
  if [ "$3" != "" ]; then          # Use value passed or caller must set it
    mem="--mem=$mem"
  fi
  if [ "$realhost" != "" ]; then
    r="--nodelist=$realhost"
  else
    r=''
  fi
  if [ "$cmdtorun" != "" ]; then
    cmdtorun="$cmdtorun"
  else
    echo "No command to execute provided to SubmitJob"
    exit 3
  fi
  l=(`/usr/cluster/bin/sbatch -p $slurmp $mem $r $cores $q --workdir=$console -J $bamid-$me --output=$console/$bamid-$me.out $cmdtorun`)
  if [ "$?" != "0" ]; then
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp $mem $r $cores $q --workdir=$console -J $bamid-$me --output=$console/$bamid-$me.out $cmdtorun" >> $console/$bamid-$me.out
    Fail "Failed to submit command to SLURM - $l" | tee -a $console/$bamid-$me.out
    exit 1
  fi
  # Job was submitted, save job details
  if [ "${l[0]}" = "Submitted" ]; then      
    echo `date` $me ${l[3]} $slurmp $qos $mem >> $console/$bamid.jobids
  fi
  $topmedcmd mark $bamid $markverb submitted
  exit
}

#------------------------------------------------------------------
# Subroutine:
#   SubmitSRUN(bamid, qos, mem, cmdtorun)    Submit this script to SLURM
#   Uses global variables: $me $markverb $console $slurp $realhost $cores
#------------------------------------------------------------------
function SubmitSRUN {
  local bamid=$1
  local qos=$2
  local mem=$3
  local cmdtorun="$4"
  if [ "$qos" != "" ]; then        # Use value passed or caller must set it
    q="--qos=$qos"
  fi
  if [ "$3" != "" ]; then          # Use value passed or caller must set it
    mem="--mem=$mem"
  fi
  if [ "$realhost" != "" ]; then
    r="--nodelist=$realhost"
  else
    r=''
  fi
  if [ "$cmdtorun" != "" ]; then
    cmdtorun="$cmdtorun"
  else
    echo "No command to execute provided to SubmitJob"
    exit 3
  fi
  cd $console || exit 5
  l=(`/usr/cluster/bin/srun -p $slurmp $mem $r $cores $q -J $bamid-$me --output=$console/$bamid-$me.out $cmdtorun`)
  if [ "$?" != "0" ]; then
    echo "CMD=/usr/cluster/bin/sbatch -p $slurmp $mem $r $cores $q --workdir=$console -J $bamid-$me --output=$console/$bamid-$me.out $cmdtorun" >> $console/$bamid-$me.out
    Fail "Failed to submit command to SLURM - $l" | tee -a $console/$bamid-$me.out
    exit 1
  fi
  # Job was submitted, save job details
  if [ "${l[0]}" = "Submitted" ]; then      
    echo `date` $me ${l[3]} $slurmp $qos $mem >> $console/$bamid.jobids
  fi
  $topmedcmd mark $bamid $markverb submitted
  exit
}


