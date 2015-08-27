#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name: topmed_status.pl
#
# Description:
#   Use this program to automatically update the states for runs
#   from NHLBI TopMed centers.
#   This program is expected to run as a crontab job.
#
# ChangeLog:
#   $Log: topmed_status.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use File::Basename;
use Cwd;
use Getopt::Long;
use DBIx::Connector;

my($me, $mepath, $mesuffix) = fileparse($0, '\.pl');
(my $version = '$Revision: 1.0 $ ') =~ tr/[0-9].//cd;

#   These represent states the run can be in. column name to letter
my %attributes2letter = (
    datearrived => 'A',
    datebai => 'I',
    datemd5ver => '5',
    datemapping => 'M',
    dateqplot => 'Q',
    datebackup => 'B',
    datecram => 'C',
    datecp2ncbi => 'N',
);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    topdir => '/net/topmed/incoming/topmed',
    verbose => 0,
);
my ($dbc);                      # For access to DB
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose topdir=s center=s runs=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$me$mesuffix [options] runstatus\n" .
        "\nVersion $version\n" .
        "Update the status fields in the database.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift(@ARGV);

my $dbh = DBConnect($opts{realm});
if ($opts{verbose}) { print "$me$mesuffix Version $version, realm=$opts{realm}\n"; }

#--------------------------------------------------------------
#   dateXXXX columns are either a time when something was done
#   or a flag indicating varying states. For example with verify:
#   datemd5ver not defined   - nothing every happened
#   datemd5ver > 10    verify successfully
#   datemd5ver < 0     started to verify
#   datemd5ver = -1    verify failed
#   datemd5ver = 0     verify requested
#   datemd5ver = 2     verify submitted to be done (not done for arrived)
#   datemd5ver = 1     verify cancelled
#--------------------------------------------------------------

#--------------------------------------------------------------
#   Calculate overall status for each run 
#--------------------------------------------------------------
if ($fcn eq 'runstatus') {
    #   Get all the known centers in the database
    my @attributes = keys %attributes2letter;
    my $attributes_str = join(',', @attributes);
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        #   For each run, calculate a status for all the bam files in the run
        my $runsref = GetRuns($cid);
        if (! %$runsref) { next; }          # No directories we want
        foreach my $runid (keys %{$runsref}) {
            my $numberbams = $runsref->{$runid};
            #   Get list of all bams for this run
            my $sql = "SELECT bamid,$attributes_str FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            my %counts = ();
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;      # Data for this BAM
                #   Get counts of states we are interested for each attribute
                foreach my $a (@attributes) {
                    if (! defined($href->{$a})) { next; }
                    if ($href->{$a} eq '') { next; }
                    if ($href->{$a} == 1)  { next; }    # cancelled
                    #   These go into done, processing or neither states
                    if ($href->{$a} > 10)  { $counts{$a}{done}++;   next; } # yea!
                    if ($href->{$a} == -1) { $counts{$a}{failed}++; next; } # boo!
                    $counts{$a}{processing}++;
                }
            }
            #   Figure out the status for this run
            #   Status consists of a letter and a summary value
            #   The letter represents the step (e.g. A for arrived)
            #   the summary value is one of done processing or unknown
            my $s = '';
            my $n;
            foreach my $a (@attributes) {
                if (! exists($counts{$a}{failed}))     { $counts{$a}{failed} = 0; }
                if (! exists($counts{$a}{done}))       { $counts{$a}{done} = 0; }
                if (! exists($counts{$a}{processing})) { $counts{$a}{processing} = 0; }
                #   If we have done all these, it's done
                if ($counts{$a}{failed} > 0) {
                    $s .= "$attributes2letter{$a}=failed,";
                    next;
                }
                if ($counts{$a}{processing} > 0) {
                    $s .= "$attributes2letter{$a}=processing,";
                    next;
                }
                if ($counts{$a}{done} > 0 && $counts{$a}{done} < $numberbams) {
                    $s .= "$attributes2letter{$a}=processing,";
                    next;
                }
                if ($counts{$a}{done} >= $numberbams) {
                    $s .= "$attributes2letter{$a}=done,";
                    next;
                }
                $s .= "$attributes2letter{$a}=unknown,";
            }
            if ($opts{verbose}) { warn "Status runid $runid=$s\n"; }
            #   Update status for this run 
            $sql = "UPDATE $opts{runs_table} SET status='$s' WHERE runid=$runid";
            $sth = DoSQL($sql);
        }
    }
    exit;
}

die "Invalid request '$fcn'. Try '$me$mesuffix --help'\n";

#==================================================================
# Subroutine:
#   GetCenters - Get list of centers
#
# Arguments:
#   none
#
# Returns:
#   Reference to hash of center ids  to center names
#==================================================================
sub GetCenters {
#    my ($d) = @_;
    my %center2name = ();

    #   Get all the known centers in the database
    my $sql = "SELECT centerid,centername FROM $opts{centers_table}";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { warn "$me$mesuffix No centers found\n"; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{center} && $href->{centername} ne $opts{center}) { next; }
        $center2name{$href->{centerid}} = $href->{centername};
    }
    if ($opts{center}) {            # Only do one center
        if (! %center2name) { die "$me$mesuffix Center '$opts{center}' unknown\n"; }
        print "$me$mesuffix - Running on center '$opts{center}' only\n";
    }
    return \%center2name;
}

#==================================================================
# Subroutine:
#   GetRuns - Get list of runs for a center and the BAM count
#
# Arguments:
#   cid = center id
#
# Returns:
#   Reference to hash of run ids  to run dirnames
#==================================================================
sub GetRuns {
    my ($cid) = @_;
    my %runs2count = ();

    my $sql = "SELECT runid,dirname,bamcount FROM $opts{runs_table} WHERE centerid=$cid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) {
        if ($opts{verbose}) { warn "$me$mesuffix Found no runs for center '$cid'\n"; }
        return \%runs2count;
    }
    my %theseruns = ();
    if ($opts{runs}) {              # We only want some runs
        $opts{runs} =~ s/,/ /g;
        foreach my $r (split(' ', $opts{runs})) {
            $theseruns{$r} = 1;
            if (! $opts{onlyrun}) { 
                print "$me$mesuffix - Using only run '$r'\n";
                $opts{onlyrun}++;
            }
        }
    }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{runs} && (! exists($theseruns{$href->{dirname}}))) { next; }
        $runs2count{$href->{runid}} = $href->{bamcount};
    }
    return \%runs2count;
}

#==================================================================
# Subroutine:
#   DBConnect($realm)
#
#   Connect to our database using realm '$realm'. Return a DB handle.
#   Get the connection information from DBIx::Connector
#   Fully qualified realm file may be provided
#==================================================================
sub DBConnect {
    my ($realm) = @_;
    if (! $realm) { return 0; }
    #   Get the connection information FROM DBIx::Connector
    #   Fully qualified realm file may be provided
    if ($realm =~ /^(\/.+)\/([^\/]+)$/) {
        my $d = $1;
        $realm = $2;
        $dbc = new DBIx::Connector(-realm => $realm, -connection_dir => $d,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    else {
        $dbc = new DBIx::Connector(-realm => $realm,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    return $dbc->connect();
}

#==================================================================
# Subroutine:
#   DoSQL - Execute an SQL command
#
# Arguments:
#   sql - String of SQL to run
#   die - boolean if we should die on error
#
# Returns:
#   SQL handle for subsequent MySQL actions
#   Does not return if error detected
#==================================================================
sub DoSQL {
    my ($sql, $die) = @_;
    if (! defined($die)) { $die = 1; }
    if ($opts{verbose}) { warn "DEBUG: SQL=$sql\n"; }
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    if ($DBI::err) {
        if (! $die) { return 0; }
        die "$me$mesuffix SQL failure: $DBI::errstr\n  SQL=$sql\n";
    }
    return $sth;
}
#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_monitor.pl - Find runs that need some action

=head1 SYNOPSIS

  topmed_monitor.pl verify

=head1 DESCRIPTION

Use program as a crontab job to find runs that have not had certain
steps completed. When one is found a request is queued to
complete the step.

This process will be most successful if this program is run
after one might expect the previous step has completed.
For instance for verifying a run, run this in the early morning when you
might expect all the data has arrived.

=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-help>

Generates this output.

=item B<-runs NAME[,NAME,...]>

Specifies a specific set of runs,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<arrive | verify | backup | qplot | report | cp2ncbi>

Directs this program to look for runs that have not been through the process name
you provided and to queue a request they be verified.

=back


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

