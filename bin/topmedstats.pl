#!/usr/bin/perl
###################################################################
#
# Name: topmedstats.pl
#
# Description:
#   Use this program to collect statistics from the jobids files
#   This updates the stepstats database to create statistics 
#   about each step.
#
# ChangeLog:
#   $Log: topmedstats.pl,v $
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
use Date::Parse;

use POSIX qw(strftime tmpnam);

my %STEPS = (                       # List of steps we collect data for
    verify => 1,
    qplot => 1,
    cram => 1,
    ncbiexpt => 1,
    ncbiorig => 1,
    ncbib37 => 1,
    gcepush => 1,
    gcepull => 1,
    b38 => 1,
    bcf => 1,
    backup => 1,
    gcecopy => 1,
    gcecpbcf => 1,
    awscopy => 1,
    fix => 1,
);
my $NOTSET    = 0;            # Not set
my $REQUESTED = 1;            # Task requested
my $SUBMITTED = 2;            # Task submitted to be run
my $STARTED   = 3;            # Task started
my $DELIVERED = 19;           # Data delivered, but not confirmed
my $COMPLETED = 20;           # Task completed successfully
my $CANCELLED = 89;           # Task cancelled
my $FAILEDCHECKSUM = 98;      # Task failed, because checksum at NCBI bad
my $FAILED    = 99;           # Task failed

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    stats_table => 'stepstats',
    bamfiles_table => 'bamfiles',
    outdir => '/net/topmed/working/topmed-output/',
    summaryfile => 'XMLfiles/latest_loaded_files.txt.gz',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose 
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options] jobid yyyy/mm/dd";
    warn "$Script [options] jobid yyyy/mm/dd\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

my $dbh = DBConnect($opts{realm});

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'jobid')    { Jobids(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   Jobids($yyyymmdd)
#
#   Read the jobids files for this particular day and sets
#   stats for all the steps completed on that day
#==================================================================
sub Jobids {
    my ($yyyymmdd) = @_;
    my %mon2number = (
        'Jan' => '01','Feb' => '02','Mar' => '03','Apr' => '04','May' => '05','Jun' => '06',
        'Jul' => '07','Aug' => '08','Sep' => '09','Oct' => '10','Nov' => '11','Dec' => '12');
    my %daytomon = (
        '01' => 'Jan','02' => 'Feb','03' => 'Mar','04' => 'Apr','05' => 'May','06' => 'Jun',
        '07' => 'Jul','08' => 'Aug','09' => 'Sep','10' => 'Oct','11' => 'Nov','12' => 'Dec');
    my $tmpfile = "/tmp/Jobids.$$";

    if ($yyyymmdd eq 'today') {
        $yyyymmdd = strftime('%Y/%m/%d', localtime);
    }
    if ($yyyymmdd !~ /^(20\d\d)\/(\d\d)\/(\d\d)$/) {
        die "$Script Invalid date '$yyyymmdd' syntax. Try '$Script -help'\n";
    }
    my ($y, $m, $d) = ($1, $2, $3);
    my $mon = $daytomon{$m};                # 01 is Jan
    my $ymd = $y . $m . $d;                 # Used to compare in database
    chdir($opts{outdir}) ||
        die "$Script - Unable to CD to '$opts{outdir}: $!\n";

    #   Get counts, but just once
    my ($sql, $sth, $href, $rowsofdata);
    $sql = "SELECT yyyymmdd FROM $opts{stats_table} WHERE yyyymmdd='$yyyymmdd'";
    $sth = DoSQL($sql, 0);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) {
        my ($count, $errcount) = (0, 0);             # This day's counts
        $sql = "SELECT count(*) FROM $opts{bamfiles_table} WHERE state_verify=$COMPLETED";
        $sth = DoSQL($sql, 0);
        $href = $sth->fetchrow_hashref;
        if ($href->{'count(*)'}) { $count = $href->{'count(*)'}; }
        $sql = "SELECT count(*) FROM $opts{bamfiles_table} WHERE state_ncbiorig=$FAILED OR state_ncbib37=$FAILED";
        $sth = DoSQL($sql, 0);
        $href = $sth->fetchrow_hashref;
        if ($href->{'count(*)'}) { $errcount = $href->{'count(*)'}; }
        #   INSERT this row, from now on we can do UPDATE
        $sql = "INSERT INTO $opts{stats_table} " .
            "(yyyymmdd,bamcount,errcount) VALUES ('$yyyymmdd',$count,$errcount)";
        $sth = DoSQL($sql);
        print "Created record for '$yyyymmdd' [$count,$errcount]\n";
    }
    

    #   Get all fields from jobids files for the date of interest
    #   Lines could look like:
    #     10164.jobids:Wed Jan 6 02:03:39 EST 2016 expt 17999839 ok secs
    #     10164.jobids:DOW 01 06 12:12:12 EST 2016 sexpt 17999839 ok 1 secs
    #   Get all jobids because some jobs take more than one day to complete
    #   This scales because old jobids files are eventually deleted
    system("grep -e E.T *.jobids > $tmpfile") &&
        die "$Script - Unable to get a list of *.jobids\n";
    open(IN, $tmpfile) ||
        die "$Script - Unable to read file '$tmpfile'\n";
    unlink($tmpfile);
    my %jobdata = ();                   # Hash of completed jobs by date
    my %startedjobs = ();               # Save start of job date        
    while (my $l = <IN>) {    
        my @c = split(' ',$l);          # Break into columns
        if ($c[6] eq 'redo') { next; }
        if ($c[6] eq 'ncbib37') { next; }
        if ($c[6] eq 'ncbiorig') { next; }
        if (! exists($STEPS{$c[6]})) { warn "$Script - Ignoring step '$c[6]' from line:  $l"; }
        if ($c[6] eq 'sexpt') { $c[6] = 'expt'; }   # Correct dumb names
        if ($c[6] eq 'ncbiorig') { $c[6] = 'orig'; }
        my $d = $c[5] . $mon2number{$c[1]} . sprintf('%02d',$c[2]);
        if (exists($startedjobs{$c[7]})) {   # Already found this jobid, calculate time
            my $t = str2time("$c[1] $c[2] $c[3] $c[5]");
            #print "finish $c[7] => $c[1] $c[2] $c[3] $c[5]  [$t]\n";
            $t = $t - $startedjobs{$c[7]};
            if ($t < 0) { $t = -$t; }
            $jobdata{$c[7]} = "$d $c[6] $t";  # Save end date, step and duration for job
            next;
        }
        #   New jobid, save start time
        $startedjobs{$c[7]} = str2time("$c[1] $c[2] $c[3] $c[5]");
        #print "initialize $c[7] => $c[1] $c[2] $c[3] $c[5]  [$startedjobs{$c[7]}]\n";
    }
    close(IN);
    if (! %jobdata) { die "No step data found for '$yyyymmdd'\n"; }

    #   Now calculate stats for each kind of step
    my %data = ();
    foreach my $j (keys %jobdata) {
        my ($d, $step, $t) = split(' ', $jobdata{$j});
        if (! exists($data{$d}{$step})) {
            $data{$d}{$step}{count} = 0;
            $data{$d}{$step}{totaltime} = 0;
        }
        $data{$d}{$step}{count}++;
        $data{$d}{$step}{totaltime} += $t;
    }

    #   Maybe show what we figured out
    my @keydates = sort keys %data;
    if ($opts{verbose}) {
        foreach my $d (@keydates) {
            if ($d ne $ymd) { next; }      # Only show date of interest
            foreach my $step (keys %{$data{$d}}) {
                print "    $step  $data{$d}{$step}{count} $data{$d}{$step}{totaltime} ave=" .
                    int($data{$d}{$step}{totaltime}/$data{$d}{$step}{count}) . " secs\n";
            }
        }
        die "Verbose specified, no database update made\n";
    }

    my $steps = '';
    $sql = '';
    foreach my $d (@keydates) {
        if ($d ne $ymd) { next; }      # Only show date of interest
        #if ($d lt '20170805') { next; }        # Hack to update a 
        my @k  = keys %{$data{$d}};
        $steps = join(',', @k);
        $sql = '';
        foreach my $step (@k) {
            if ($data{$d}{$step}) {
                $sql .= "count_$step=$data{$d}{$step}{count},";
                $sql .= "avetime_$step=" . int($data{$d}{$step}{totaltime}/$data{$d}{$step}{count}) . ',';
            }
        }
        if ($sql) {
            chop($sql);                     # Remove ,
            $sql = "UPDATE $opts{stats_table} SET $sql WHERE yyyymmdd='" . 
                substr($d,0,4) . '/' . substr($d,4,2) . '/' . substr($d,6,2) . "'";
            $sth = DoSQL($sql);
            print "Updated data for '$d' [$steps]\n";
        }
    }
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedstats.pl - Collect statistics from the jobids files

=head1 SYNOPSIS

  topmedstats.pl jobid 2015/03/02       # Collect stats for a particular day

=head1 DESCRIPTION

This program is used to collect statistics from the jobids files
and update the statistics database used to generate plots.

See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<jobid yyyy/mm/dd>

Collect statistics from the jobids files for a particular date.

=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

