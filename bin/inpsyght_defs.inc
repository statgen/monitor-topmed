#
#       Define variables for PROJECT INPSYGHT action scripts
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
slurmp=$PROJECT
build=38
realhost=''

cores=''                            # Set to same as gsutil if it will be used
gsutil="gsutil -m"
gsshareoption="-u topmed-1366"      # Used to copy to SHARE
#gsutilbig='gsutil -o GSUtil:parallel_composite_upload_threshold=150M'
export BOTO_CONFIG=/net/topmed/working/shared/tpg_gsutil_config.txt

#incominguri='gs://topmed-recabs'
#bcfuri='gs://topmed-bcf'

#------------------------------------------------------------------
# Subroutine:
#   RandomRealHost(bamid, realhost)   Set realhost based on bamid
#------------------------------------------------------------------
function RandomRealHost {
  #local bamid=$1                # Normal processing, pick realhost from bamid
  realhost=$PROJECT
}
