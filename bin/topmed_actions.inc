#
#	Define variables that all TOPMED action scripts need defined
#	Requires PROJECT to be set
#
source /usr/cluster/$PROJECT/bin/${PROJECT}_defs.inc

#==================================================================
#   Define useful functions to be used by the action scripts
#==================================================================

#------------------------------------------------------------------
# Subroutine:
#   Trap SIGTERM signals forced by a timeout
#   Generate a Fail message so we know this timed out and can restart it
#------------------------------------------------------------------
trap "finish" USR1
function finish {
  #   Derive bamid from jobname
  l=(`echo $SLURM_JOB_NAME|sed -e 's/-/ /'`)
  bamid=${l[0]}
  m="Job timed out  bamid=$bamid jobid=$SLURM_JOB_ID `date`"
  Fail "$m"
}

#------------------------------------------------------------------
# Subroutine:
#   MayIRun(name, bamid)    Is bamid permitted to run for script 'name'
#
#   Does not return if permission is denied
#------------------------------------------------------------------
function MayIRun {
  $topmedpermit test $1 $2
  rc=$?
  if [ "$rc" = "4" ]; then
    exit 4
  fi
  if [ "$rc" = "5" ]; then
    exit 5
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
#   Warn(msg)
#   Writes msg to stderr cause I can never remember how to in bash
#------------------------------------------------------------------
function Warn {
  >&2 echo $*
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
    $topmedcmd -emsg "$msg" mark $bamid $markverb failed
  else
    >&2 echo $msg
    #>&2 echo "Task for $bamid FAILED, no markverb set"
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
    $topmedcmd mark $bamid $markverb started
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
    $topmedcmd -emsg "" mark $bamid $markverb completed
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
  $topmedcmd set $bamid $col "$val"
  if [ "$?" != "0" ]; then
    Fail "Unable to set column '$col' to '$val' for bamid '$bamid'" 
  fi
}

#------------------------------------------------------------------
# Subroutine:
#   val=`GetDB(bamid, col)`
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
#   val=`GetNWDID(bamid)`       # Prints NWDID
#------------------------------------------------------------------
function GetNWDID {
  local bamid=$1
  nwdid=`$topmedcmd show $bamid expt_sampleid`
  if [ "$nwdid" = "" ]; then
    Fail "Unable to find the NWDID for '$bamid'" 
  fi
  echo $nwdid
}

#------------------------------------------------------------------
# Subroutine:
#   MyRealHost(bamid, type)   Set realhost where type of input file lives
#------------------------------------------------------------------
function MyRealHost {
  local bamid=$1
  local type=$2
  realhost=`$topmedpath whathost $bamid $type`
}

#------------------------------------------------------------------
# Subroutine:
#   CheckRGMap(bamid)     rc=0 if rgmap not needed, rc=1 if needed
#------------------------------------------------------------------
function CheckRGMap {
  local bamid=$1
  
  /usr/cluster/$PROJECT/bin/topmedrgmap.pl -check $bamid /dev/null
  return $?
}

#------------------------------------------------------------------
# Subroutine:
#   b=`GetBuild(bamid,file)`    Determine build for a file
#------------------------------------------------------------------
function GetBuild {
  local bamid=$1
  local file=$2

  build=''
  l=(`$samtools view -H $file | grep '@SQ' | cut -f 1-3 2>/dev/null | head -2`)
  if [ "${l[2]}" = "LN:249250621" -a "${l[5]}" = "LN:243199373" ]; then
    build=37
  else
    if [ "${l[2]}" = "LN:248956422" -a "${l[5]}" = "LN:242193529" ]; then
      build=38
    fi
  fi
  if [ "$build" = "" ]; then
    Fail "Unable to determine build for '$file' (L2=${l[2]} L5=${l[5]}"
  fi
  >&2 echo "File '$file' [$bamid] is from build $build"
  echo $build
}

#------------------------------------------------------------------
# Subroutine:
#   ipath=`ValidateIndex(bamid,file)`
#
#   file must exist
#
#	print indexfile_path valid_or_not
#------------------------------------------------------------------
function ValidateIndex {
  local bamid=$1
  local file=$2         # BAM or CRAM

  #  To verify an index, we search for a string - depending on the extension
  ext="${file##*.}"
  buildstr='notset'
  indexfile='notset'
  if [ "$ext" = "bam" ]; then
    indexfile=$file.bai
    buildstr="3:148,100,000-148,110,000"
  fi
  if [ "$ext" = "cram" ]; then
    indexfile=$file.crai
    build=`GetBuild $bamid $file`
    if [ "$build" = "37" ]; then
      buildstr="hs37d5:35,450,000-35,477,943"
      if [ "$project" = "inpsyght" ]; then  # Broad use differ 37 ref for this project
        buildstr="GL000192.1:544,000-547,490"
      fi
    fi
    if [ "$build" = "38" ]; then
      buildstr="chr19_KI270938v1_alt"
    fi
  fi
  if [ "buildstr" = 'notset' ]; then
    Fail "ValidateIndex: Unable to determine string to search for in header of $file"
  fi

  #   Index must exist and not be null
  if [ ! -f $indexfile ]; then
    echo "$indexfile invalid"
    return
  fi
  if [ ! -s $indexfile ]; then
    echo "$indexfile invalid"
    return
  fi

  #   Index exists, let's see if it is any good
  n=`$samtools view $file $buildstr | wc -l`
  if [ "$n" = "" -o "$n" -lt "100" ]; then
    #   This might be a trashed reference index for samtools, if so remove it
    if [ "$me" != "" ]; then
      consfile="$console/$bamid-$me.out"
      a=`grep 'cram_ref_load: Assertion' $consfile`     # Could fail if $me not defined
      if [ "$a" != "" ]; then
        rm -rf $HOME/.cache/hts-ref/*/*
         >&2 echo "Samtools cram_ref_load error, removed dirty reference cache and retrying"
      fi
    fi 
    echo "$indexfile invalid"
    return
  fi
  echo "$indexfile valid"
  return
}

#------------------------------------------------------------------
# Subroutine:
#   CreateIndex(bamid,file)    Create CRAM/BAM index
#
#	samtools can fail in multiple ways and never tell us
#   This script attempts to centralize all known attempts
#   to check for samtools errors for an index file
#   so we can possibly insure the index file will work  (POS)
#------------------------------------------------------------------
function CreateIndex {
  stat="stat --printf=%Y --dereference"

  local bamid=$1
  local file=$2

  craifile=$file.crai           # This is file to be created
  nameoffile=CRAI
  x=`echo $file | grep bam`
  if [ "$x" != "" ]; then
    craifile=$file.bai
    nameoffile=BAI
  fi

  timecram=`$stat $file 2>/dev/null`
  timecrai=`$stat $craifile 2>/dev/null`

  if [ "$timecram" = "" ]; then
    Fail "INPUT file for '$bamid' does not exist: $file"
  fi
  #   Check if existing CRAI is valid
  if [ "$timecrai" != "" ]; then
    s=(`ValidateIndex $bamid $file`)
    craifile="${s[0]}"            # Path to index file
    status="${s[1]}"              # valid or invalid
    if [ "$status" != "valid" ]; then
      echo "Index file '$indexfile' was invalid, attempt to re-create it"
      rm -f $craifile
    fi
  fi

  #   CRAI must be older than CRAM
  calcchecksum=n
  if [ "$timecrai" = "" ]; then
    timecrai=`expr $timecram - 10`
  fi
  if [ "$timecram" -ge "$timecrai" ]; then
    echo "No $nameoffile or old one, creating index file '$craifile'"
    $samtools index $file 2>&1
    if [ "$?" != "0" ]; then
      Fail "Unable to create CRAI index for '$file'"
    fi
    echo "Created $nameoffile index for '$file'"
    calcchecksum=y
  fi
 
  s=(`ValidateIndex $bamid $file`)
  status="${s[1]}"              # valid or invalid
  if [ "$status" != "valid" ]; then
    rm -f $craifile
    Fail "Index file '$craifile' remains invalid. I give up"
  fi

  if [ "$nameoffile" = "BAI" ]; then
    return
  fi

  # Need checksum for CRAM index
  x=`grep recab $file`
  if [ "$x" = "" ]; then
    return
  fi
  md5=`GetDB $bamid b38craichecksum`
  #  Calculate checksum?
  if [ "$calcchecksum" = "n" -a "$md5" != "" ]; then
    return
  fi
  v=`CalcMD5 $bamid $craifile`
  SetDB $bamid b38craichecksum $v
  echo "Calculated checksum for CRAI"
}

#------------------------------------------------------------------
# Subroutine:
#   CreateBCFIndex(bamid,bcffile)    Create BCF index for remapped bcf, update db
#------------------------------------------------------------------
function CreateBCFIndex {
  stat="stat --printf=%Y --dereference"

  local bamid=$1
  local bcffile=$2
  csifile=$bcffile.csi

  timebcf=`$stat $bcffile 2>/dev/null`
  timecsi=`$stat $csifile 2>/dev/null`

  if [ "$timebcf" = "" ]; then
    Fail "BCF file for '$bamid' does not exist: $bcffile"
  fi
  #   CSI must be older than BCF
  if [ "$timecsi" = "" ]; then
    timecsi=`expr $timebcf - 10`
  fi
  if [ "$timebcf" -ge "$timecsi" ]; then
    echo "No CSI or old CSI, creating index file '$csifile'"
    $bcftools index $bcffile
    if [ "$?" != "0" ]; then
      Fail "Unable to create CSI for $bamid' for: $bcffile"
    fi
    echo "Created CSI index for BCF '$bcffile'"
  fi
  if [ ! -s $csifile ]; then
    rm -f $csifile
    Fail "Zero length CSI file detected for $bamid for: $bcffile"
  fi
}

#------------------------------------------------------------------
# Subroutine:
#   v=`CalcFlagstatFromFile(file)`    Prints flagstat for a bam or cram
#
#	Calculate the flagstat from a flagstat file
#------------------------------------------------------------------
function CalcFlagstatFromFile {
  local file=$1

  #   File has results of flagstat call, get fields and calculate 
  a=(`grep 'in total' $file`)
  if [ "${a[0]}" = "" ]; then
    Fail "Unable to get in total from '$file'"
  fi
  intotal=${a[0]}
  a=(`grep 'secondary' $file`)
  if [ "${a[0]}" = "" ]; then
    ${a[0]} = 0
  fi
  secondary=${a[0]}
  a=(`grep 'supplementary' $file`)
  if [ "${a[0]}" = "" ]; then
    ${a[0]} = 0
  fi
  supplementary=${a[0]}
  flagstat=`expr $intotal - $secondary - $supplementary`
  echo $flagstat
}

#------------------------------------------------------------------
# Subroutine:
#   v=`CalcFlagstat(bamid,file)`    Prints flagstat for a bam or cram
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

  #   Outfile has results of flagstat, print flagstat
  CalcFlagstatFromFile $of
  rm -f $of
}

#------------------------------------------------------------------
# Subroutine:
#   v=`CalcMD5(bamid, file)`    Prints md5 for a file
#------------------------------------------------------------------
function CalcMD5 {
  local bamid=$1
  local file=$2

  md5=`md5sum $file`
  if [ "$md5" = "" ]; then
    Fail "Unable to calculate MD5 for '$file' [$bamid]"
  fi
  checksum=($md5)
  echo ${checksum[0]}
}

#------------------------------------------------------------------
# Subroutine:
#   SubmitJob(bamid, qos, mem, cmdtorun)    Submit this script to SLURM
#   Uses global variables: $me $markverb $console $slurp $realhost $cores $timeout
#	Uses environment variable: SBATCH_OPT
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
  if [ "$timeout" != "" ]; then
    t="--time=$timeout --signal=B:USR1@60"
  else
    t=''
  fi
  if [ "$cmdtorun" != "" ]; then
    cmdtorun="$cmdtorun"
  else
    echo "No command to execute provided to SubmitJob"
    exit 3
  fi
  l=(`$slurmpath/sbatch -p $slurmp $mem $r $t $cores $q $SBATCH_OPT --workdir=$console -J $bamid-$me --output=$console/$bamid-$me.out $cmdtorun`)
  if [ "$?" != "0" ]; then
    echo "CMD=$slurmpath/sbatch -p $slurmp $mem $r $cores $q $SBATCH_OPT --workdir=$console -J $bamid-$me --output=$console/$bamid-$me.out $cmdtorun" >> $console/$bamid-$me.out
    Fail "Failed to submit command to SLURM - $l" | tee -a $console/$bamid-$me.out
    exit 1
  fi
  # Job was submitted, save job details
  if [ "${l[0]}" = "Submitted" ]; then      
    echo `date` $me ${l[3]} $slurmp $qos $mem >> $console/$bamid.jobids
  fi
  if [ "$markverb" != "" ]; then
    $topmedcmd mark $bamid $markverb submitted
  fi
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
  l=(`slurmpath/srun -p $slurmp $mem $r $cores $q -J $bamid-$me --output=$console/$bamid-$me.out $cmdtorun`)
  if [ "$?" != "0" ]; then
    echo "CMD=$slurmpath/sbatch -p $slurmp $mem $r $cores $q --workdir=$console -J $bamid-$me --output=$console/$bamid-$me.out $cmdtorun" >> $console/$bamid-$me.out
    Fail "Failed to submit command to SLURM - $l" | tee -a $console/$bamid-$me.out
    exit 1
  fi
  # Job was submitted, save job details
  if [ "${l[0]}" = "Submitted" ]; then      
    echo `date` $me ${l[3]} $slurmp $qos $mem >> $console/$bamid.jobids
  fi
  if [ "$markverb" != "" ]; then
    $topmedcmd mark $bamid $markverb submitted
  fi
  exit
}
