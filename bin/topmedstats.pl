#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
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
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/perl5";
use TopMed_Get;
use My_DB;
use Getopt::Long;

use POSIX qw(strftime tmpnam);

my %STEPS = (                       # List of steps we collect data for
    verify => 1,
    bai => 1,
    qplot => 1,
    cram => 1,
    expt => 1,
    orig => 1,
    b37 => 1,
    b38 => 1,
);
my $NOTSET    = 0;            # Not set
my $REQUESTED = 1;            # Task requested
my $SUBMITTED = 2;            # Task submitted to be run
my $STARTED   = 3;            # Task started
my $DELIVERED = 19;           # Data delivered, but not confirmed
my $COMPLETED = 20;           # Task completed successfully
my $CANCELLED = 89;           # Task cancelled
my $FAILED    = 99;           # Task failed

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
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
        "  or\n" .
        "$Script [options] summary yyyy/mm/dd\n" .
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
if ($fcn eq 'summary')  {
    Summary(@ARGV, "$opts{outdir}/$opts{summaryfile}");
    exit;
}

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   Summary($yyyymmdd, $file)
#
#   Read the NCBI summary file, identifying all BAMs loaded on a date
#   The loadedbamcount column is the total loaded counts each day
#   Once the database is initialized nicely with SummaryHistory
#   we only need get the count for this day and then add it to
#   the previous database entry - avoiding re-processing everything
#   everyday -- and relying the the NCBI summary file is consistent
#==================================================================
sub Summary {
    my ($yyyymmdd, $file) = @_;

    my $origsum = 0;
    my $b37sum = 0;
    my $b38sum = 0;

    #   Read the summary file in reverse order. This allows us to get the LAST state for a file
    my $cmd = "tac $file";
    if ($file =~ /\.gz$/) { $cmd = "gunzip -c $file | tac - 2> /dev/null |"; }
    my $rc = open(IN, $cmd);
    if (! $rc) { die "$Script - Unable to read data from command: $cmd\n"; }
    while (<IN>) {
        chomp();
        #   [0]=protected [2]=date [3]=NWDnnnnnn [5]=checksum
        #   [9]=error|loaded|received [10]=BAM|CRAM|UKNOWN [12]=message
        my @cols = split(' ',$_);
        if ($cols[0] ne 'protected') { next; }
        if ($cols[3] !~ /(NWD\d+)\.(.+)/) { next; } # Break into NWDnnnnnn and suffix
        my ($nwd, $suffix) = ($1, $2);
        if ($suffix eq 'expt.xml') { next; }
        my $date = substr($cols[2],0,10);           # Date of this action
        $date =~ s/-/\//g;
        if ($date ne $yyyymmdd) { next; }
        if ($opts{verbose}) { print "Processing data for date=$date\n"; }

        if (($suffix eq 'src.bam' && $cols[10] eq 'BAM') ||
            ($suffix eq 'src.cram' && ($cols[10] eq 'CRAM' || $cols[10] eq 'UNKNOWN'))) {
            $origsum++;             
            next;
        }

        if (($suffix eq 'recal.bam' && $cols[10] eq 'BAM') ||   # From a bug of mine
            ($suffix eq 'recal.cram' && ($cols[10] eq 'CRAM' || $cols[10] eq 'UNKNOWN')) || # My bug
            ($suffix eq 'recal.37.bam' && ($cols[10] eq 'BAM' || $cols[10] eq 'UNKNOWN')) || # My bug
            ($suffix eq 'recal.37.cram' && ($cols[10] eq 'CRAM' || $cols[10] eq 'UNKNOWN'))) {
            $b37sum++;             
            next;
        }

        if (($suffix eq 'recal.38.cram' && ($cols[10] eq 'CRAM' || $cols[10] eq 'UNKNOWN'))) {
            $b38sum++;           
            next;
        }
    }
    close(IN);
    print "Loaded $origsum,$b37sum,$b38sum BAMs on $yyyymmdd\n";

    my $errorigcount = 0;
    my $errb37count = 0;
    my $errb38count = 0;
    my $sql = "SELECT count(*) FROM $opts{bamfiles_table} WHERE state_ncbiorig=$FAILED";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if ($rowsofdata > 0) {
        my $href = $sth->fetchrow_hashref;
        $errorigcount  = $href->{'count(*)'};
    }
    $sql = "SELECT count(*) FROM $opts{bamfiles_table} WHERE state_ncbib37=$FAILED";
    $sth = DoSQL($sql);
    $rowsofdata = $sth->rows();
    if ($rowsofdata > 0) {
        my $href = $sth->fetchrow_hashref;
        $errb37count  = $href->{'count(*)'};
    }
    $sql = "SELECT count(*) FROM $opts{bamfiles_table} WHERE state_ncbib38=$FAILED";
    $sth = DoSQL($sql);
    $rowsofdata = $sth->rows();
    if ($rowsofdata > 0) {
        my $href = $sth->fetchrow_hashref;
        $errb38count  = $href->{'count(*)'};
    }

    DoSQL("UPDATE $opts{stats_table} SET " .
        "loadedorigbamcount=$origsum," .
        "loadedb37bamcount=$b37sum," .
        "loadedb38bamcount=$b38sum," .
        "errorigcount=$errorigcount," .
        "errb37count=$errb37count," .
        "errb38count=$errb38count" .
        " WHERE yyyymmdd='$yyyymmdd'");
}

#==================================================================
# Subroutine:
#   SummaryHistory($file)
#
#   Read the NCBI summary file, identifying all BAMs loaded to a date
#   and generate the loadedbamcount in stepstats for the past
#
#   This was used just once to create history for stepstats
#==================================================================
sub SummaryHistory {
    my ($yyyymmdd, $file) = @_;

    #   Find all NWDIDnnnnn.src.bam in
    #   protected 3563129 2014-09-06T09:24:12 NWDnnnnn.src.bam 10844214270 93ed94e9918b868d0ecd7009c3e427e8 = = = loaded BAM etc
    #     or
    #   protected NWD792235-remap.37.run.xml ... error RUN_XML - some_error_message
    my %loadedbydate = ();
    my ($loaddate, $rc);
    if ($file =~ /\.gz$/) { $rc = open(IN, "gunzip -c $file |"); }
    else { $rc = open(IN, $file); }
    if (! $rc) { die "$Script - Unable to read file '$file': $!\n"; }
    while (<IN>) {
        if (! /protected\s+\d+\s+(20\S+)T.+\s+NWD/) { next; }
        $loaddate = $1;
        $loaddate =~ s/-/\//g;
        if ($loaddate lt $yyyymmdd) { next; }
        if (/protected\s+\d+\s+20.+\s+NWD\S+bam\s+.+=\s+=\s+=\s+loaded\sBAM/) {
            if (! exists($loadedbydate{$loaddate})) { $loadedbydate{$loaddate} = 0; }
            $loadedbydate{$loaddate}++;
            next;
        }
    }
    close(IN);
    my @ldates = sort keys %loadedbydate;   # List of dates when bams were loaded

    my $sql = "SELECT yyyymmdd FROM $opts{stats_table}";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "No rows?\n"; }
    my @sqlupdates = ();
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        my $yyyymmdd = $href->{yyyymmdd};
        #   Get count of all loaded bams until this date
        my $sum = 0;
        foreach my $d (@ldates) {
            if ($d gt $yyyymmdd) { last; }
            $sum += $loadedbydate{$d};
        }
        push @sqlupdates,"UPDATE $opts{stats_table} SET loadedbamcount=$sum WHERE yyyymmdd='$yyyymmdd';";
    }
    foreach (@sqlupdates) {
        print "$_\n";
    }
}

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
    my $tmpfile = "/tmp/Jobids.$$";

    if ($yyyymmdd eq 'today') {
        $yyyymmdd = strftime('%Y/%m/%d', localtime);
    }

    if ($yyyymmdd !~ /^(20\d\d)\/(\d\d)\/(\d\d)$/) {
        die "$Script Invalid date '$yyyymmdd' syntax. Try '$Script -help'\n";
    }
    my ($y, $m, $d) = ($1, $2, $3);
    chdir($opts{outdir}) ||
        die "$Script - Unable to CD to '$opts{outdir}: $!\n";

    #   Get counts of bams and errors sending bams, but just once
    my ($sql, $sth, $href, $rowsofdata);
    $sql = "SELECT yyyymmdd FROM $opts{stats_table} WHERE yyyymmdd='$yyyymmdd'";
    $sth = DoSQL($sql, 0);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) {
        my ($bamcount, $errcount) = (0, 0);             # Today's counts
        $sql = "SELECT count(*) FROM $opts{bamfiles_table} WHERE state_md5ver=$COMPLETED";
        $sth = DoSQL($sql, 0);
        $href = $sth->fetchrow_hashref;
        if ($href->{'count(*)'}) { $bamcount = $href->{'count(*)'}; }
        $sql = "SELECT count(*) FROM $opts{bamfiles_table} WHERE state_ncbiorig=$FAILED OR " .
            "state_ncbib37=$FAILED OR state_ncbib38=$FAILED";
        $sth = DoSQL($sql, 0);
        $href = $sth->fetchrow_hashref;
        if ($href->{'count(*)'}) { $errcount = $href->{'count(*)'}; }
        #   INSERT this row
        $sql = "INSERT INTO $opts{stats_table} " .
            "(yyyymmdd,bamcount,errcount) VALUES ('$yyyymmdd',$bamcount,$errcount)";
        $sth = DoSQL($sql);
        print "Created record for '$yyyymmdd' [$bamcount,$errcount]\n";
    }

    #   Get all OK files from jobids files and extract dates we want
    #   Lines could look like:
    #     10164.jobids:Wed Jan 6 02:03:39 EST 2016 expt 17999839 ok secs
    #     10164.jobids:DOW 01 06 12:12:12 EST 2016 sexpt 17999839 ok 1 secs
    system("grep -w ok *.jobids> $tmpfile") &&
        die "$Script - Unable to get a list of *.jobids\n";
    open(IN, $tmpfile) ||
        die "$Script - Unable to read file '$tmpfile'\n";
    unlink($tmpfile);
    my %data = ();                      # Hash of step->[count][times]
    while (my $l = <IN>) {    
        my @c = split(' ',$l);          # Break into columns
        if ($#c == 9 && $c[9] eq 'secs') {    # Sometimes # secs missing
            $c[10] = $c[9];
            $c[9] = '1';
        }
        if ($c[1] =~ /\D/) {
            if (exists($mon2number{$c[1]})) { $c[1] = $mon2number{$c[1]}; }
            else { die "$Script - Unable to get month '$c[1]' in line:\n  $l"; }
        }
    
        #   If this is the date we want, collect details for this step
        if (! ($c[1] == $m && $c[2] == $d && $c[5] == $y)) { next; }

        #   Correct step names from old/bad data
        if ($c[6] eq 'sexpt') { $c[6] = 'expt'; }
        if ($c[6] eq 'ncbiorig') { $c[6] = 'orig'; }
        if (! exists($STEPS{$c[6]})) { die "$Script - unknown step '$c[6]' from line:\n  $l"; }
        if (! exists($data{$c[6]})) { $data{$c[6]}{count} = 0; $data{$c[6]}{totaltime} = 0; }
        $data{$c[6]}{count}++;
        $data{$c[6]}{totaltime} += $c[9];
    }
    close(IN);
    if (! %data) { die "No data found for '$yyyymmdd'\n"; }
    my @keys = sort keys %data;
    if ($opts{verbose}) {
        foreach my $s (@keys) {
            print "    $s  $data{$s}{count} $data{$s}{totaltime} ave=" . int($data{$s}{totaltime}/$data{$s}{count}) . " secs\n";
        }
    }

    $sql = "UPDATE $opts{stats_table} SET ";
    foreach my $s (@keys) {
        $sql .= "count_$s=$data{$s}{count},";
        $sql .= "avetime_$s=" . int($data{$s}{totaltime}/$data{$s}{count}) . ',';
    }
    chop($sql);
    $sql .= " WHERE yyyymmdd='$yyyymmdd'";
    $sth = DoSQL($sql);
    print "Updated data for '$yyyymmdd' [" . join(',', @keys) . "]\n";
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedcmd.pl - Update the database for NHLBI TopMed

=head1 SYNOPSIS

  topmedcmd.pl mark 33 arrived completed   # BAM has arrived
  topmedcmd.pl mark NWD00234  arrived completed   # Same BAM has arrived
  topmedcmd.pl unmark 33 arrived           # Reset BAM has arrived

  topmedcmd.pl set 33 jobidqplot 123445    # Set jobidqplot in bamfiles
 
  topmedcmd.pl where 2199                  # Returns path to bam, path to backup, bamname

  topmedcmd.pl permit add bai braod 2015oct18   # Stop bai job submissions for a run
  topmedcmd.pl permit remove 12             # Remove a permit control
  topmedcmd.pl permit test bai 4567         # Test if we should submit a bai job for one bam

=head1 DESCRIPTION

This program supports simple commands to set key elements of the NHLBI database.
The queue of tasks is kept in a MySQL database.
See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-runs NAME[,NAME,...]>

Specifies a specific set of runs on which to run the action,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program are grouped into several groups which are used
to deal with specific sets of information in the monitor databases.

B<mark bamid dirname  [verb] [state]>
Use this to set the state for a particular BAM file.
You may specify the bamid or the NWDID.
Mark will set a date for the process (e.g. arrived sets state_arrive)
and unmark will set that entry to NULL.
The list of verbs and states can be seen by B<perldoc topmedcmd.pl>.

B<permit enable/disable operation center run>
Use this to control the database which allows one enable or disable topmed operations
(e.g. backup, verify etc) for a center or run.
Use B<all> for all centers or all runs or all operations.

B<permit test operation bamid>
Use this to test if an operation (e.g. backup, verify etc) may be submitted 
for a particular bam.

B<set bamid columnname value>
Use this to set the value for a column for a particular BAM file.

B<send2ncbi filelist>
Use this to copy data to NCBI with ascp.

B<show arrived>
Use this to show the bamids for all BAMs that are marked arrived.

B<step yyyy/mm/dd stepname slurmid seconds>
Update the step database with time for various steps so they can be easily plotted.
This only applies to successful steps.
You can extract the data from steptime with a command like:

B<mysql --user=sqlnhlbiro --host=f-db --password=PASSWORD --batch \>
  B<-e 'select step,stepdate,seconds from steptime' nhlbi | \>
  B<grep -v stepdate | sed "s/\t/,/g"> /tmp/step.csv>

B<unmark bamid [verb]>
Use this to reset the state for a particular BAM file to the default
database value.

B<whatnwdid bamid|NWDnnnnnn>
Use this to get some details for a particular bam.

<where bamid>
Use this to display the directory of the bam file, 
the path to the backup direcotry,
the name of the bam without any path,
the real host where the file leaves (e.g. B<topmed3>),
and the index of the hostname (e.g. B<2> for topmed2).


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

