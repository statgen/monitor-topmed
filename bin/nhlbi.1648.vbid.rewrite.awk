
	#  Convert Fan Zhang's verifyBamId output format to Hyun's output format.

	#  This finds the NWD identifier either from the input file name or from 
	#  an explicit argument passed as "-v NWD=NWDxxxyyy", with preference for 
	#  the explicit argument if both are available.

	#  Run as:  "awk -f nhlbi.1648.vbid.rewrite.awk  [-v NWD=NWD999296]  <inputfile>"

BEGIN	     {	FS = OFS = "\t";  }

/Alpha:/     {	NWD = length(NWD) < 9 ? substr(FILENAME, 1, 9) : NWD;
	 	outfile = NWD ".vb.selfSM";
	 	k = split($0, val, ":");

	 	if (length(NWD) != 9 || substr(NWD, 1, 3) != "NWD")
	 	 {  outfile = "/dev/null";
	 	    print "vbid.rewrite.awk - cannot find valid NWD identifier."  >> "/dev/stderr";
	 	    print "Either call with argument \"-v NWD=NWDxxxyyy\" or use" >> "/dev/stderr";
	 	    print "an explicit filename starting with NWDxxxyyy."         >> "/dev/stderr";
	 	    exit 1;  } }

END	     {	print "#SEQ_ID", "RG", "CHIP_ID", "#SNPS", "#READS", "AVG_DP",
	 	"FREEMIX", "FREELK1", "FREELK0", "FREE_RH", "FREE_RA", "CHIPMIX",
	 	"CHIPLK1", "CHIPLK0", "CHIP_RH", "CHIP_RA", "DPREF", "RDPHET",
	 	"RDPALT"  > outfile;
	 	print NWD, "-", "-", "-", "-", "-", val[2], "-", "-", "-", "-",
	 	 "-", "-", "-", "-", "-", "-", "-", "-"  >>  outfile;  }



   #  Expected input format:

	#  Hello, World!
	#  PCs in OptimizaLLK():
	#  PC1:-0.0048166	PC2:-0.0134297
	#  Alpha:0.000517336

   #  Format to match on output (long lines have been wrapped):

	#  #SEQ_ID	RG	CHIP_ID	#SNPS	#READS	AVG_DP	FREEMIX	FREELK1	FREELK0	 	 	\
	#  FREE_RH	FREE_RA	CHIPMIX	CHIPLK1	CHIPLK0	CHIP_RH	CHIP_RA	DPREF	RDPHET	RDPALT
	#  NWD778639	ALL	NA	1391791	5688821	40.87	0.07037	17194740.10	19206365.62	\
	#  0.50187	0.00290	NA	NA	NA	NA	NA	NA	NA	NA
