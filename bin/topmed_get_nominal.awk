
	#  Return mean and std.dev of mapped insert size, 
	#  number of reads on which this is based, and then for 
	#  end 1 and end 2:  end #, # reads, # base pairs.

	#  Run this as:    foreach file ( *.recal.bam )
	#     echo $file '  ' `samtools view $file 
	#	  | awk -f sardinia.89.awk`  >>  target.lengths

BEGIN	     {	FS = OFS = "\t";
	 	base = base ? base : 300;
	 	rows = rows ? rows : 2e+5;	}

$1 !~ /^@/  &&  $2 < 900  &&  $5 > 22  &&  $6 !~ /[DNHP]/	\
	     {
	 	cdv = split($6, cigar, "[MSIDX]");
	 	cdv = cdv == 1 ? 2 : cdv;
	 	len =  0;
	 	for (i = 1; i < cdv; i++)
	 	     len += cigar[i];
	 	end = ($2 - ($2 % 64)) / 64;
	 	rln[end] += (len / base);
	 	rct[end]++;
	 	isz = ($9 >= 0 ? $9 : -$9) / base;
	 	if (isz > 0  &&  isz < 5)
	 	 {  prr += (isz - 1);
	 	    pss += (isz - 1) * (isz - 1);
	 	    pnn++;	 	     } }

pnn >= rows  {	exit;  }

END	     {	isz = base + base * (prr / pnn);
	 	isd = base * sqrt((pss - prr * prr / pnn) / pnn);

	 	printf "\t%7.2f\t%7.2f\t %6.2f K", isz, isd, pnn / 1e+3;
	 	for (end in rct)
	 	    printf "\t     %d\t%d\t%7.3f", end, rct[end],
	 	 	base * (rln[end] / rct[end]);
	 	printf "\n";	 	 	 	 	}


 0	     {	if (end == 1)	 	raa += (len / base);
	 	else if (end == 2)	rbb += (len / base);
	 	else	 	 	rcc += (len / base);

	 	printf "\t%d\t%d\t%d\n", cdv, end, len;
	 	if (cdv > 3)
	 	print $6;	 	 	}





