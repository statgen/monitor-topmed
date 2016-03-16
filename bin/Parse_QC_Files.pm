package Parse_QC_Files;
#==================================================================
# Parse_QC_Files.pm
#==================================================================
use strict;
use warnings;

use vars qw(@EXPORT_OK @EXPORT @ISA);
use Exporter;
use POSIX qw(strftime);

@ISA = qw(Exporter);
@EXPORT_OK = qw(parseRVector parseQCFiles parseQCColumns);
@EXPORT    = qw(parseRVector parseQCFiles parseQCColumns);

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
    push(@metrics,sprintf("%.1lf",$F[6]*100));

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

Parse_QC_Files.pm

=head1 SYNOPSIS

  use Parse_QC_Files;
  
=head1 DESCRIPTION

This provides common functions extract values of interest from qplot output

=head1 Functions

=over 4

=item B<parseQCColumns()>

Returns an array of the names of the data returned.

=item B<parseQCFiles($genocheckfile, $statfile, $rfile)>

Returns an array of the qplot output

=item B<parseRVector($line)>

Internal function to parse vector created by R.

=back

=head1 AUTHOR

Written by HM Kang I<E<lt>hmkang@umich.eduE<gt>> in 2015 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut
