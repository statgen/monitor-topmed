
	#  Parse the output from 'samtools view -H' and return a region specifier 
	#  for the last 30 kb from the last reference sequence name in the .bam 
	#  or .cram file.  Useful for checking whether a .cram file is truncated.

	#  Run this, for example, as:

	#  samtools view <pathto.cram> `samtools view -H <pathto.cram> | awk -f  \
	#	 	 	 	  nhlbi.1530.cram.end.region.awk` | wc

	#  and check for a message on stderr that says:
	#    "[W::hts_close] EOF marker is absent. The input is probably truncated."


BEGIN	     {	FS = OFS = "\t";  }

$1 ~ /@SQ/   {	for (i = 2; i <= NF; i++)
	 	 {  if (sub("SN:", "", $i))
	 	 	chr = $i;
	 	    if (sub("LN:", "", $i))
	 	 	end = 0 + $i;	 	      } }

END	     {	beg = end > 30000 ? end - 30000 : 1;
	 	print chr ":" beg "-" end;	 	}

