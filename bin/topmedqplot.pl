#!/usr/bin/perl
###################################################################
#
# Name: topmedqplot.pl
#
# Description:
#   Use this program to update the NHLBI database with data created by qplot
#
#   Parts of this were written by HM Kang (hmkang@umich.eduE) in 2015
#
# ChangeLog:
#   $Log: topmedqplot.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use FindBin qw($Bin $Script);
use lib (
  qq($FindBin::Bin),
  qq($FindBin::Bin/../lib),
  qq($FindBin::Bin/../lib/perl5),
  qq($FindBin::Bin/../local/lib/perl5),
  qq(/usr/cluster/topmed/lib/perl5),
  qq(/usr/cluster/topmed/local/lib/perl5),
);
use My_DB;
use Getopt::Long;
use POSIX qw(strftime);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    qcdata_table => 'qc_results',
    bamfiles_table => 'bamfiles',
    password => 'please',               # Password to be entered for --remove  (yes, not secure)
    verbose => 0,
);
my $dbh;

Getopt::Long::GetOptions( \%opts,qw(
    help finddups realm=s remove=s verbose
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 1 || $opts{help}) {
    warn "$Script [options] dirname nwdid\n" .
        "\nUpdate the topmed database with QPLOT data\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $dirname = shift @ARGV;
my $nwdid   = shift @ARGV;

#--------------------------------------------------------------
#   Special admin hook - show all duplicates in the database
#--------------------------------------------------------------
#   Show all the duplicates. There should be none, but ...
if ($opts{finddups}) {
    my $aref = Dups();
    if ($#$aref < 0) { print "$Script - No duplicates found!\n"; exit; }
    print "Found " . ($#$aref+1) . " duplicates\n";
    if ($opts{verbose}) { print join(' ', @$aref) . "\n"; }
    exit;
}

#   Generate the SQL to remove duplicate records
if ($opts{remove}) {
    if ($opts{remove} ne $opts{password}) {
        die "$Script - password '$opts{remove}' invalid\n";
    }

    my $aref = Dups();
    if ($#$aref < 0) { print "$Script - No duplicates found!\n"; exit; }
    foreach my $bamid (@$aref) {
        my $sql = "SELECT id from $opts{qcdata_table} WHERE bam_id=$bamid ORDER by id";
        my $sth = DoSQL($sql);
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) { print "$Script - Oops, now there is no duplicate '$bamid'?\n"; next; }
        my @dupbamids = ();
        for (my $i=0; $i<$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref();
            push @dupbamids,$href->{id};            
        }
        my $keepbamid = pop(@dupbamids);
        if ($opts{verbose}) { print "Keep $keepbamid, remove " . join(',',@dupbamids) . "\n"; }
        else {
            foreach my $id (@dupbamids) {
                print "DELETE FROM $opts{qcdata_table} WHERE id=$id;\n";
            }
        }
    }
    exit;
}

#   See if the NWDID we were provided is valid
DBConnect($opts{realm});
my $sth = DoSQL("SELECT bamid,build FROM $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'");
my $rowsofdata = $sth->rows();
if (! $rowsofdata) { die "$Script - NWDID '$nwdid' is unknown\n"; }
my $href = $sth->fetchrow_hashref;
my $bamid = $href->{bamid};
my $build = $href->{build};

#--------------------------------------------------------------
#   Extract qplot values or die trying
#--------------------------------------------------------------
my $f = $dirname . '/' . $nwdid;
my $self  = `ls $f*.vb.selfSM`;
chomp($self);
my $stats = `ls $f*.qp.stats`;
chomp($stats);
my $r     = `ls $f*.qp.R`;
chomp($r);
my @metrics = parseQCFiles($self, $stats, $r);

if ($opts{verbose}) { print "$nwdid parsed:\n"; foreach my $v (@metrics) { print "$v\n"; } }

my @colnames = parseQCColumns();
if ($#colnames != $#metrics) {
    die "$Script - Number of columns [$#colnames] returned by parseQCFiles [$#metrics] did not match\n";
}
#   Figure out indexes to columns of interest
my ($pct_freemix_index, $pct_genome_dp10_index, $pct_genome_dp20_index, $mean_depth_index) = (-1, -1, -1, -1);
for (my $i=0; $i<=$#colnames; $i++) {
    if ($colnames[$i] eq 'PCT_FREEMIX')     { $pct_freemix_index = $i; next; }
    if ($colnames[$i] eq 'PCT_GENOME_DP10') { $pct_genome_dp10_index = $i; next; }
    if ($colnames[$i] eq 'PCT_GENOME_DP20') { $pct_genome_dp20_index = $i; next; }
    if ($colnames[$i] eq 'MEAN_DEPTH')      { $mean_depth_index = $i; }
}
foreach my $i ($pct_freemix_index, $pct_genome_dp10_index, $pct_genome_dp20_index, $mean_depth_index) {
    if ($i < 0) { die "$Script - Unable to find index to colname I need\n"; }
}

#   The new qplot introduces pct_genome_dp20_index which does not apply when build!=38
#   Now set the quality columns to be included in the SQL
my ($qc_flagged, $qc_pass, $qc_fail) = (0, 0, 0);   # Default values for SQL columns
if ($build >= 38) {
    if ($metrics[$pct_freemix_index]<=3 &&
        $metrics[$pct_genome_dp10_index]>=95 &&
        $metrics[$pct_genome_dp20_index]>=90 &&         # For new qplot
        $metrics[$mean_depth_index]>=25 &&
        $metrics[$mean_depth_index]<30) {
        $qc_flagged = 1;
    }
    if ($metrics[$pct_freemix_index]<=3 &&
        $metrics[$pct_genome_dp10_index]>=95 &&
        $metrics[$pct_genome_dp20_index]>=90 &&         # For new qplot
        $metrics[$mean_depth_index]>=25) {
        $qc_pass = 1;
    }
    if ($metrics[$pct_freemix_index]>3 ||
        $metrics[$pct_genome_dp10_index]<95 ||
        $metrics[$pct_genome_dp20_index]<90 ||         # For new qplot
        $metrics[$mean_depth_index]<25) {
        $qc_fail = 1;
    }
}
else {          # This applies to build 37 etc.
    if ($metrics[$pct_freemix_index]<=3 &&
        $metrics[$pct_genome_dp10_index]>=95 &&
        $metrics[$mean_depth_index]>=25 &&
        $metrics[$mean_depth_index]<30) {
        $qc_flagged = 1;
    }
    if ($metrics[$pct_freemix_index]<=3 &&
        $metrics[$pct_genome_dp10_index]>=95 &&
        $metrics[$mean_depth_index]>=25) {
        $qc_pass = 1;
    }
    if ($metrics[$pct_freemix_index]>3 ||
        $metrics[$pct_genome_dp10_index]<95 ||
        $metrics[$mean_depth_index]<25) {
        $qc_fail = 1;
    }
}
push @colnames,'qc_fail','qc_pass','qc_flagged';
push @metrics, $qc_fail,$qc_pass,$qc_flagged;

#--------------------------------------------------------------
#   Update database with values we obtained
#--------------------------------------------------------------
my $date = shift(@metrics);
$_ = shift(@colnames);              # Ignore name of this first column
my $sql;

#   See if there already is a qc result for this NWDID
$sth = DoSQL("SELECT bam_id FROM $opts{qcdata_table} WHERE bam_id=$bamid");
$rowsofdata = $sth->rows();
if ($rowsofdata) {                  # This already exists, do UPDATE
    $sql = "UPDATE $opts{qcdata_table} SET created_at='$date',";
    for (my $i=0; $i<=$#metrics; $i++) {        # No need to specify bam_id
        $sql .= "$colnames[$i]=$metrics[$i],";
    }
    chop($sql);
    $sql .= " WHERE bam_id=$bamid";
}
else {                              # First time, do INSERT
    #   First value is a string, rest are floats
    $sql = "INSERT INTO $opts{qcdata_table} (bam_id,created_at," . lc(join(',',@colnames)) .
        ") VALUES($bamid,'$date',";
    for (my $i=0; $i<=$#metrics; $i++) {
        $sql .= $metrics[$i] . ',';
    }
    chop($sql);
    $sql .= ')';
}
if ($opts{verbose}) { print "SQL=$sql\n"; }
$sth = DoSQL($sql);
if (! $sth) { die "$Script Failed to update database. SQL=$sql\n"; }

print "Updated QC database data for '$nwdid'\n";
exit;

#==================================================================
# Subroutine:
#   $aref = Dups()
#
#   Return a reference to an array of duplicate qc_result entries
#==================================================================
sub Dups {
    #my ($href, $ncbihref) = @_;
    my @dups = ();
    
    $dbh = DBConnect($opts{realm}) ||
        die "$Script - Unable to connect to database\n";
    my $sth = DoSQL("SELECT count(bam_id) as count,bam_id FROM $opts{qcdata_table} group by bam_id having count > 1");
    my $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        my $href;
        for (my $i=0; $i<$rowsofdata; $i++) {
            $href = $sth->fetchrow_hashref();
            push @dups, $href->{bam_id};
        }
    }
    return \@dups;
}

#==================================================================
# Subroutine:
#   parseQCColumns - Return column names of data returned
#
#   Usage: @colnames = parseColumns();
#
# Arguments: None
#
# Returns:
#   An array of qplot column names (not necessarily database column names)
#==================================================================
sub parseQCColumns {
    my @colnames = ( 
        'DATE',               # date of QC assessment
        'PCT_FREEMIX',        # % contamination estimated from verifyBamID
        'N_READS_M',          # Number of total reads in millions
        'PCT_MAPPED',         # % of mapped reads
        'PCT_MQ0',            # % of mapped reads with zero mapping quality
        'PCT_PAIRED',         # % of paired reads
        'PCT_PROP_PAIRED',    # % of properly paired reads
        'MAPPED_Gb',          # Total gigabases of mapped reads, excluding (1) QC-failed (2) secondary, (3) duplicated reads, or (4) 'N' bases in the reference genome.
        'Q20_Gb',             # Total gigabases of mapped, excluding (1)-(4) above
        'PCT_Q20_BASE',       # % of Q20 bases among all mapped bases
        'MEAN_DEPTH',         # Average mapped depth dividing by non-'N' bases, excluding (1)-(4) above.
        'PCT_GENOME_COV',     # % of genome covered, excluding the regions with 'N' bases
        'ISIZE_MODE',         # Mode of insert size
        'ISIZE_MEDIAN',       # Median insert size
        'PCT_DUPS',           # Percentage of duplicated reads among mapped reads
        'PCT_GENOME_DP5',     # % of genome coverged at depth 5 or greater, excluding the regions with 'N' bases
        'PCT_GENOME_DP10',    # % of genome coverged at depth 10 or greater, excluding the regions with 'N' bases
        'PCT_GENOME_DP20',    # % of genome coverged at depth 20 or greater, excluding the regions with 'N' bases
        'PCT_GENOME_DP30',    # % of genome coverged at depth 30 or greater, excluding the regions with 'N' bases
        'VMR_DEPTH',          # Variance-to-mean ratio of depth, as an indicator of over-disperson of depth (e.g. due to GC bias)
        'DEPTH_Q10',          # Average mapped depth for Q10 bases, excluding (1)-(4) above, dividing by non-'N' bases.
        'DEPTH_Q20',          # Average mapped depth for Q20 bases, excluding (1)-(4) above, dividing by non-'N' bases.
        'DEPTH_Q30',          # Average mapped depth for Q30 bases, excluding (1)-(4) above, dividing by non-'N' bases.
        'RAW_BASE_Gb',        # (Approximate) Total number of sequenced gigabases (N_READS_M x MAX_CYCLE)
        'PCT_OVERLAP_READS',  # % of overlapping reads where insert size is equal or smaller than the cycle.
        'PCT_OVERLAP_BASES',  # % of overlapping bases, counting only overlapping fraction of reads based on insert size
        'ISIZE_IQR',          # Inter-quartile range of insert size
        'ISIZE_STDEV',        # Standard deviation of insert size (capped at 10kb).
        'GC_DEPTH_0_1',       # GC-normalized depth for top 1 percentile
        'GC_DEPTH_1_5',       # GC-normalized depth for top 1-5 percentile
        'GC_DEPTH_5_25',      # GC-normalized depth for top 5-25 percentile
        'GC_DEPTH_25_75',     # GC-normalized depth for top 25-75 percentile
        'GC_DEPTH_75_95',     # GC-normalized depth for top 75-95 percentile
        'GC_DEPTH_95_99',     # GC-normalized depth for top 95-99 percentile
        'GC_DEPTH_99_100',    # GC-normalized depth for top 99-100 percentile
        'LIBRARY_SIZE_M',     # Library size (# of reads in library in million) = (TOTAL_READS) / (0-log(1-PCT_DUPS/100))
    );
    return @colnames;
}

#==================================================================
# Subroutine:
#   parseQCFiles - Extract values of interest from qplot output
#
#   Usage: @metrics = parseQCFiles($genocheckfile, $statfile, $rfile);
#
# Arguments: output files from qplot
#   selfsm - selfSM, e.g. NWD321439.vb.selfSM
#   stats  - states  e.g. NWD321439.qp.stats
#   rfile  - R file, e.g. NWD321439.qp.R
#
#
#   Example:
#       $id = 'NWD321439';
#       @metrics = parseQCFiles("$id.vb.selfSM","$id.qp.stats","$id.qp.R");
#
# Returns:
#   An array of interesting stats
#==================================================================
sub parseQCFiles {
    my ($genocheckfile,$statfile,$rfile) = @_;
    my @metrics = ();    

    #   Get second line of first file
    chomp($genocheckfile);
    open(IN, $genocheckfile) ||
        die "Unable to read file '$genocheckfile': $!\n";
    my $line = <IN>;
    $line = <IN>;
    close(IN);
    
    my @F = split(/\s+/,$line);
    push(@metrics,sprintf("%.3lf",$F[6]*100));

    #   Read second file
    #   @flags indicates the fields in statfile we want
    my @flags = (0,1,1,0,0,1,0,1,1,1,1,1,1,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    chomp($statfile);
    open(IN, $statfile) ||
        die "Unable to read file '$statfile': $!\n";   
    for (my $i=0; <IN>; ++$i) {
	    if ( $flags[$i] == 1 ) {
	        my ($name,$val) = split;
	        push(@metrics,$val);
	    }
    }
    close(IN);

    my $nReadsM = $metrics[1];
    my $meanDepth = $metrics[9];
    my $pctGenomeCov = $metrics[10];
    my $pctDups = $metrics[13];

    #   Read third filem get data from just some lines
    chomp($rfile);
    open(IN,$rfile) ||
        die "Unable to read file '$rfile': $!\n";   
    my @phredX = ();
    my @phredZ = ();
    my @cycleX = ();
    my @cycleZ = ();
    my @gcbiasX = ();
    my @gcbiasY = ();
    my @gcbiasZ = ();
    my @isizeX = ();
    my @isizeY = ();
    my @depthX = ();
    my @depthY = ();
    for (my $i=0; <IN>; ++$i) {
	    if ( $i == 11 ) {
	        @phredX = parseRVector($_);
	    }
	    elsif ( $i == 15 ) {
	        @phredZ = parseRVector($_);
	    }
	    elsif ( $i == 42 ) {
	        @cycleX = parseRVector($_);
	    }
	    elsif ( $i == 46 ) {
	        @cycleZ = parseRVector($_);
	    }
	    elsif ( $i == 67 ) {
	        @gcbiasZ = parseRVector($_);
	    }
	    elsif ( $i == 68 ) {
	        @gcbiasX = parseRVector($_);
	    }
	    elsif ( $i == 69 ) {
	        @gcbiasY = parseRVector($_);
	    }
	    elsif ( $i == 94 ) {
	        @isizeX = parseRVector($_);
	    }
	    elsif ( $i == 95 ) {
	        @isizeY = parseRVector($_);	    
	    }
	    elsif ( $i == 133 ) {
	        @depthX = parseRVector($_);
	    }
	    elsif ( $i == 134 ) {
	        @depthY = parseRVector($_);
	    }
    }

    #   Calculate % DP5, DP10, DP20, DP30 (among DP>0) and VMR
    my @cumDPs = (0) x ($#depthX+1);
    $cumDPs[$#depthX] = $depthY[$#depthX];
    my $cntDP = 0;
    my $sumDP = 0;
    my $sumsqDP = 0;
    for (my $i=$#depthX-1; $i>=0; --$i) {
	    $cumDPs[$i] += $cumDPs[$i+1] + $depthY[$i];
	    $cntDP += $depthY[$i];
	    $sumDP += ( $depthX[$i] * $depthY[$i]);
	    $sumsqDP += ($depthX[$i] * $depthX[$i] * $depthY[$i]);
    }
    $cntDP /= ($pctGenomeCov/100);
    my $meanDP = $sumDP / $cntDP;
    my $varDP = $sumsqDP / $cntDP - $meanDP * $meanDP;
    my $VMR = $varDP / $meanDP;
    
    my $pctGenomeDP5 = $pctGenomeCov * $cumDPs[4] / $cumDPs[0];
    my $pctGenomeDP10 = $pctGenomeCov * $cumDPs[9] / $cumDPs[0];
    my $pctGenomeDP20 = $pctGenomeCov * $cumDPs[19] / $cumDPs[0];    
    my $pctGenomeDP30 = $pctGenomeCov * $cumDPs[29] / $cumDPs[0];

    #   Calculate %Q10, Q20, Q30 bases
    my ($q0, $q10, $q20, $q30) = (0,0,0,0);
    for (my $i=0; $i<@phredX; ++$i) {
	    $q0 += $phredZ[$i];
	    if ( $phredX[$i] >= 10 ) {
	        $q10 += $phredZ[$i];
	        if ( $phredX[$i] >= 20 ) {
		        $q20 += $phredZ[$i];
		        if ( $phredX[$i] >= 30 ) {
		            $q30 += $phredZ[$i];		    
		        }
	        }
	    }
    }
    my $q10DP = $meanDepth * $q10 / $q0;
    my $q20DP = $meanDepth * $q20 / $q0;
    my $q30DP = $meanDepth * $q30 / $q0;
    my $maxCycle = $cycleX[$#cycleX];

    my $rawBases = $metrics[1] * $maxCycle / 1000;

    #   Fraction of pairs with too small insert sizes
    my ($totalPairs,$overlappingPairs,$overlappingFracPairs) = (0,0,0);
    for (my $i=0; $i<@isizeX; ++$i) {
	    $totalPairs += $isizeY[$i];
	    if ( $isizeX[$i] <= $maxCycle ) {
	        $overlappingPairs += $isizeY[$i] ;
	        $overlappingFracPairs += ($isizeY[$i] * ($maxCycle - $isizeX[$i]) / $maxCycle);
	    }
    }
    my $pctOverlapReads = $overlappingPairs/$totalPairs*100;
    my $pctOverlapBases = $overlappingFracPairs/$totalPairs*100;

    #   When calculating SE, cap insert size to 1kb
    my $maxIS = 1000;
    my ($cntIS,$sumIS,$sqSumIS) = (0,0);
    my ($q1IS,$q3IS) = (0,0);
    for (my $i=0; $i<@isizeX; ++$i) {
	    $cntIS += $isizeY[$i];
	    my $x = ($isizeX[$i] > $maxIS) ? $maxIS : $isizeX[$i];
	    $sumIS += ($x * $isizeY[$i]);
	    $sqSumIS += ( $x * $x * $isizeY[$i] );
	    $q1IS = $isizeX[$i] if ( $cntIS < 0.25 * $totalPairs );
	    $q3IS = $isizeX[$i] if ( $cntIS < 0.75 * $totalPairs );	
    }

    my $gcTotal = 0;
    my @gcCutOffs = (0,0.01,0.05,0.25,0.75,0.95,0.99,1);
    my @gcDepthNum =  (0,0,0,0,0,0,0);
    my @gcDepthDen =  (0,0,0,0,0,0,0);    
    for (my $i=0; $i<@gcbiasX; ++$i) {
	    $gcTotal += $gcbiasZ[$i];
	    $gcbiasY[$i] = 1 if ( $gcbiasY[$i] eq "NA" );
    }
    my $gcSum = 0;
    for (my $i=0; $i<@gcbiasX; ++$i) {  
	    for (my $j=0; $j<$#gcCutOffs; ++$j) {  # [gcSum, gcSum + gcbiasZ[i]]
	        #   Not included at all
	        if ( $gcSum + $gcbiasZ[$i] < $gcCutOffs[$j] * $gcTotal ) { }    # no overlap, do nothing
	        elsif ( $gcSum >= $gcCutOffs[$j+1] * $gcTotal) {}               # no overlap, do nothing
	        elsif ( $gcSum >= $gcCutOffs[$j] * $gcTotal) {
		        if ( $gcSum + $gcbiasZ[$i] < $gcCutOffs[$j+1] * $gcTotal ) {    # completely contained
		            $gcDepthDen[$j] += $gcbiasZ[$i];
		            $gcDepthNum[$j] += ( $gcbiasY[$i] * $gcbiasZ[$i] );
		        }
		        else {          # Partially contained - cut at the right side
		            my $frac = ($gcCutOffs[$j+1]*$gcTotal - $gcSum)/($gcbiasZ[$i]);
		            $gcDepthDen[$j] += ( $frac * $gcbiasZ[$i] );
		            $gcDepthNum[$j] += ( $frac * $gcbiasY[$i] * $gcbiasZ[$i] );
		        }
	        }
	        else {          # gcSum < $gcCutOffs[$j] < $gcSum + gcbiasZ[$j]
		        my $frac = ($gcbiasZ[$i]+$gcSum - $gcCutOffs[$j]*$gcTotal)/($gcbiasZ[$i]);
		        $gcDepthDen[$j] += ( $frac * $gcbiasZ[$i] );
		        $gcDepthNum[$j] += ( $frac * $gcbiasY[$i] * $gcbiasZ[$i] );
	        }
	    }
	    $gcSum += $gcbiasZ[$i];
    }

	#   Get creation date of first file
    my $mindate = (stat($genocheckfile))[9];
    $mindate = (stat($rfile))[9] if ( (stat($rfile))[9] < $mindate );
    $mindate = (stat($statfile))[9] if ( (stat($statfile))[9] < $mindate );
    $mindate = sprintf(strftime("%Y/%m/%d", localtime( $mindate )));    

    my @data = ($mindate,@metrics,
	    sprintf("%.2lf",$pctGenomeDP5),
	    sprintf("%.2lf",$pctGenomeDP10),
	    sprintf("%.2lf",$pctGenomeDP20),
	    sprintf("%.2lf",$pctGenomeDP30),
	    sprintf("%.2lf",$VMR),
	    sprintf("%.2lf",$q10DP),
	    sprintf("%.2lf",$q20DP),
	    sprintf("%.2lf",$q30DP),
	    sprintf("%.2lf",$rawBases),
	    sprintf("%.2lf",$pctOverlapReads),
	    sprintf("%.2lf",$pctOverlapBases),
	    sprintf("%.1lf",$q3IS-$q1IS),
	    sprintf("%.1lf",sqrt($sqSumIS/$cntIS-$sumIS*$sumIS/$cntIS/$cntIS)),
	    sprintf("%.3lf",$gcDepthNum[0]/$gcDepthDen[0]),
	    sprintf("%.3lf",$gcDepthNum[1]/$gcDepthDen[1]),
	    sprintf("%.3lf",$gcDepthNum[2]/$gcDepthDen[2]),
	    sprintf("%.3lf",$gcDepthNum[3]/$gcDepthDen[3]),
	    sprintf("%.3lf",$gcDepthNum[4]/$gcDepthDen[4]),
	    sprintf("%.3lf",$gcDepthNum[5]/$gcDepthDen[5]),
	    sprintf("%.3lf",$gcDepthNum[6]/$gcDepthDen[6]),	    
	    sprintf("%.1lf",$nReadsM/(1e-6 - log(1-$pctDups/100))));
	#my @cols = parseQCColumns();
	#for (my $i=0; $i<=$#cols; $i++) { print " $i  $cols[$i]=$data[$i]\n"; }
	return @data;
}

#==================================================================
# Subroutine:
#   parseRVector - Parse a vector from R
#
# Arguments:
#   line = of the form c(a,b,c)
#
# Returns:
#   An array strings within () or dies
#==================================================================
sub parseRVector {
    my $line = shift;
    chomp($line);
    if ( /c\((\S+)\)/ ) {
	    my @F = split(/,/,$1);
	    return @F;
    }
    die "parseRVector: Cannot parse R vector in '$line'\n";
}

1;

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedqplot.pl - Update the database with qplot values for a NWDID sampleid

=head1 SYNOPSIS

  topmedqplot.pl /net/topmed/incoming/qc.results/broad/2015may01.single.fram NWD321439

  topmedqplot.pl -finddups  ignore ignore       # Debugging, show number qc_results duplicates
  topmedqplot.pl -finddups -verbose ignore ignore 
  topmedqplot.pl -remove=please ignore ignore   # Remove show qc_results duplicates

=head1 DESCRIPTION

This program extracts values of interest from qplot output and saves them
in the database.

See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-finddups>

Scan the database of qc results and show instances of multiple bam_id 
This is only for admins to fing bugs.
No processing of normal parameters is done (e.g. add/modify database).

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-remove=password>

Remove instances of multiple bam_id records from the database.
This generates the SQL to remove all duplicates (all but the last record).
This is only for admins to fing bugs.
No processing of normal parameters is done.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<dirname>

The directory where the qplot output resides

=item B<nwdid>

The NWDID for the data to extract. This serves as the prefix for each of the qplot files.

=item B<-verbose>

Provided for developers to see additional information.

=back


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

