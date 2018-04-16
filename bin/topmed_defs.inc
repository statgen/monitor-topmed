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
vt=/usr/cluster/topmed/bin/vt
vtref=/net/topmed/working/mapping/gotcloud/ref/hg38/hs38DH.fa
bcftools=/usr/cluster/bin/bcftools
console=/net/topmed/working/topmed-output
slurmp=topmed-working
build=38
realhost=''

cores=''                            # Set to same as gsutil if it will be used
gsutil="gsutil -m"
gsshareoption="-u topmed-1366"      # Used to copy to SHARE
gsutilbig='gsutil -o GSUtil:parallel_composite_upload_threshold=150M'
#   No need to parallel processing 
#cores="--cpus-per-task=2"         # Cores here should be same as gsutil
#gsutil='gsutil -o GSUtil:parallel_composite_upload_threshold=150M -o GSUtil:parallel_process_count=2'
incominguri='gs://topmed-recabs'
bcfuri='gs://topmed-bcf'
export BOTO_CONFIG=/net/topmed/working/shared/tpg_gsutil_config.txt

