#
#       Define variables for PROJECT TOPMED action scripts
#
topmedcmd=/usr/cluster/$PROJECT/bin/topmedcmd.pl
topmedpermit=/usr/cluster/$PROJECT/bin/topmedpermit.pl
topmedpath=/usr/cluster/$PROJECT/bin/topmedpath.pl
topmedrename=/usr/cluster/$PROJECT/bin/topmedrename.pl
topmedqplot=/usr/cluster/$PROJECT/bin/topmedqplot.pl
bam=/usr/cluster/bin/bam
samtools=/usr/cluster/bin/samtools
vt=/usr/cluster/topmed/bin/vt
vtref=/net/topmed/working/mapping/gotcloud/ref/hg38/hs38DH.fa
bcftools=/usr/cluster/bin/bcftools
console=/net/$PROJECT/working/${PROJECT}-output
slurmp=topmed-working
build=38
realhost=''

cores=''                            # Set to same as gsutil if it will be used
gsutil="gsutil -m"
gsshareoption="-u topmed-1366"      # Used to copy to SHARE
#gsutilbig='gsutil -o GSUtil:parallel_composite_upload_threshold=150M'
export BOTO_CONFIG=/net/topmed/working/shared/tpg_gsutil_config.txt

incominguri='gs://topmed-recabs'
bcfuri='gs://topmed-bcf'

#------------------------------------------------------------------
# Subroutine:
#   RandomRealHost(bamid, realhost)   Set realhost based on bamid
#------------------------------------------------------------------
function RandomRealHost {
  local bamid=$1                # Normal processing, pick realhost from bamid
  declare -A bamid2host
  bamid2host[0]=topmed
  bamid2host[1]=topmed2
  bamid2host[2]=topmed3
  bamid2host[3]=topmed4
  bamid2host[4]=topmed5
  bamid2host[5]=topmed6
  bamid2host[6]=topmed7
  bamid2host[7]=topmed9
  bamid2host[8]=topmed10
  #   See if anything is down
  f=/run/shm/RandomRealHost.down
  if [ -f $f ]; then
    while read l; do
      ll=($l)
      bamid2host[${ll[0]}]=DOWN
      #echo "Marked index ${ll[0]} [${ll[1]}] as DOWN"
    done < $f
  fi
  n=`perl -e "{ print $bamid%9 }"`
  #   See this entry is down, use first entry. If topmed down, cannot submit jobs
  if [ "${bamid2host[$n]}" = "" -o "${bamid2host[$n]}" = "DOWN" ]; then
    n=0
  fi
  if [ "$SLURMNODE" != "" ]; then	# Allow user to force node
    realhost=$SLURMNODE
  else
    realhost=${bamid2host[$n]}
  fi
}
