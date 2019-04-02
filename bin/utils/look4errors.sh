#!/bin/bash
#  look4errors.sh
#
#    Check database for states in error
#
q="'"
tmp=/tmp/look4errors.tmp
s=''
for state in state_arrive state_verify state_backup state_cram state_qplot state_gce38push state_gce38pull state_b38 state_gce38bcf state_gce38cpbcf state_gce38copy state_gcecleanup state_aws38copy; do
  s="$s $state=99 or ";
done  
s="select bamid,expt_sampleid from bamfiles where $s state_fix=99"
eval dosql.pl -realm topmed -batch -cmd ${q}$s${q} 2>/dev/null | awk '{if(NR>5)print}' | grep -v ' row' | grep -v Bye | grep -v Empty > $tmp
sz=`stat --printf %s $tmp`
if [ "$sz" -gt "3" ]; then
  echo `date '+%Y/%m/%d %H:%M'` `cat $tmp`
fi
rm -f $tmp
exit
