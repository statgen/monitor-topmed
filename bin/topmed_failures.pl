#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name: topmed_failures.pl
#
# Description:
#   Use this program to automatically find jobs which have
#   failed in some way that the master database has not
#   been updated. This will appear as a job that has
#   started, but never finishes.
#
# ChangeLog:
#   $Log: topmed_failures.pl,v $
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
use POSIX qw(strftime);

my($me, $mepath, $mesuffix) = fileparse($0, '\.pl');
(my $version = '$Revision: 1.0 $ ') =~ tr/[0-9].//cd;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my $topmedbin = '/usr/cluster/monitor/bin';
my %opts = (
    topmedcmd => "$topmedbin/topmedcmd.pl",
    sacct => '/usr/cluster/bin/sacct',
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    topdir => '/net/topmed/incoming/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    outdir => '/working/topmed-output',
    verbose => 0,
    maxjobs => 50,
);

#   This is a map of part of a column name to the name in the log file
my %col2logname = (
    md5ver => 'verify',
    backup => 'backup',
    bai => 'bai',
    qplot => 'qplot',
);
my %col2verb = (                    # Mark job as failed with topmedcmd.pl
    md5ver => 'md5verified',
    backup => 'backedup',
    bai => 'baid',
    qplot => 'qploted',
);
my @colnames = keys %col2logname;   # Used in *_Failure functions
my $sqlors = '';                    # Build OR string for SQL
foreach my $c (@colnames) { $sqlors .= 'date' . $c . '=-1 OR '; }
$sqlors = substr($sqlors,0,-4);     # Drop ' OR '

my ($dbc);                      # For access to DB
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose=n center=s runs=s resubmit maxjobs=n
    memory=s partition=s qos=s noprompt
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$me$mesuffix [options] look4failures\n" .
        "or\n" .
        "$me$mesuffix [options] resubmit\n" .
        "\nVersion $version\n" .
        "Find runs which were cancelled and we do not know about them.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift(@ARGV);

my $dbh = DBConnect($opts{realm});
if ($opts{verbose}) { print "$me$mesuffix Version $version, realm=$opts{realm}\n"; }

my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

#   Set environment variables for shell scripts that do -submit
if ($opts{memory})    { $ENV{TOPMED_MEMORY} = $opts{memory}; }
if ($opts{partition}) { $ENV{TOPMED_PARTITION} = $opts{partition}; }
if ($opts{qos})       { $ENV{TOPMED_QOS} = $opts{qos}; }

#--------------------------------------------------------------
#   Get a list of jobs that have been submitted, see if any failed
#--------------------------------------------------------------
if ($fcn =~ /^look/) {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        if ($opts{verbose}) { print $centername . "\n"; }
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles that arrived
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            if ($opts{verbose}) { print "  $dirname\n"; }
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                if (! $opts{maxjobs}) { last; }
                my $href = $sth->fetchrow_hashref;
                Find_Failure($href);
            }
        }
    }
    if ($opts{jobcount}) { print "$opts{jobcount} job failures detected\n"; }
    exit;
}

#--------------------------------------------------------------
#   Find all failed jobs and attempt to resubmit them
#--------------------------------------------------------------
if ($fcn eq 'resubmit') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        if ($opts{verbose}) { print $centername . "\n"; }
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles that arrived
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            if ($opts{verbose}) { print "  $dirname\n"; }
            #   Get list of all bams that have arrived
            my $sql = "SELECT * FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid' AND ($sqlors)";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                if (! $opts{maxjobs}) { last; }
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                if (! -f $f) { next; }          # If BAM not there, do not submit
                Resubmit_Failure($f, $href);    # Something here failed
            }
        }
    }
    if ($opts{jobcount}) { print "$opts{jobcount} failed jobs have been resubmitted\n"; }
    exit;
}

die "Invalid request '$fcn'. Try '$me$mesuffix --help'\n";

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

#==================================================================
# Subroutine:
#   Find_Failure - For any tasks which have jobs submitted
#       see if the job has failed in SLURM
#       Mark any failed step in the database
#
#   sacct might return something like:
#       7819459      2674-qplot topmed-in+  csg  1 CANCELLED+
#   and the reason why it failed is in the output (e.g. output/2721-qplot.out)
#   like this:
#       slurmstepd: *** JOB 7819459 CANCELLED AT 2015-06-22T11:37:17 
#
# Arguments:
#   href - hash reference to columns in a database row
#
# Returns:
#   nothing
#==================================================================
sub Find_Failure {
    my ($href) = @_;

    foreach my $date (@colnames) {
        my $col = 'date' . $date;
        if (! defined($href->{$col})) { next; }
        if ($href->{$col} eq '') { next; }
        if (! ($href->{$col} eq '2' || $href->{$col} < -10)) { next; }
        $col = 'jobid' . $date;
        if (! defined($href->{$col})) { next; }
        if ($href->{$col} eq '') { next; }
        #if ($href->{$col} =~ /\D/) { next; }
        #   This jobid has been submitted, but not completed
        my $cmd = $opts{sacct} . ' -j ' . $href->{$col};
        my $s = `$cmd 2>&1`;
        if ($opts{verbose} > 1) { print "Checking job $href->{$col}\n$s"; }
        if ($s !~ /CANCELLED/ && $s !~ /NODE_FAIL/  && $s !~ /FAILED/) { next; }

        #   Job failed 
        print "Jobid='$href->{$col}' $date FAILED because of SLURM\n" .
            "  BAMID=$href->{bamid} BAMNAME=$href->{bamid}\n$s\n";
        # See if we can find something in the output that tells us why
        my $logfile = $opts{outdir} . '/' . $href->{bamid} . '-' .
            $col2logname{$date} . '.out';
        if (open(IN,$logfile)) {
            while(<IN>) {
                if (/slurmstepd/) { print '  ' . $_; }
            }
            close(IN);
        }
        #   Mark this as a failed job $topmedcmd mark $bamid backedup failed
        $cmd = "$opts{topmedcmd} mark $href->{bamid} $col2verb{$date} failed";
        system($cmd) && print "Failed to mark job as failed\n  CMD=$cmd\n";
        $opts{jobcount}++;
        return;
    }
    return;
}

#==================================================================
# Subroutine:
#   Resubmit_Failure - Resubmit a failed task
#
# Arguments:
#   file - fully qualified path to BAM
#   href - hash reference to columns in a database row
#
# Returns:
#   nothing
#==================================================================
sub Resubmit_Failure {
    my ($file, $href) = @_;

    #   One or more of these tasks need to be resubmitted
    foreach my $date (@colnames) {
        my $col = 'date' . $date;
        if (! defined($href->{$col})) { next; }
        if ($href->{$col} eq '') { next; }
        if ($href->{$col} ne '-1') { next; }    # Not failed
        #   This task failed, resubmit it
        my $cmd = "$topmedbin/topmed_$col2logname{$date}.sh -submit $href->{bamid} $file";
        if ($date eq 'md5ver') {
            $cmd = "$topmedbin/topmed_$col2logname{$date}.sh -submit $href->{bamid} $href->{checksum} $file";
        }
        if (! $opts{noprompt}) {
            print "Submit this job?\n  $cmd\n  Enter 'y' or 'n' or 'q': ";
            $_ = <STDIN>;
            if ($_ eq "q\n") { exit; }
            if ($_ ne "y\n") { print "Nothing submitted\n"; return; }
        }
        system($cmd) && die "Failed to submit job\n  cmd=$cmd\n";
        print "Submitted $date job for bamid=$href->{bamid}\n";
        $opts{jobcount}++;
        $opts{maxjobs}--;
        if ($opts{maxjobs} == 0) { print "Maximum limit of jobs that can be submitted has been reached\n"; }
    }
}

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
#   GetRuns - Get list of runs for a center
#
# Arguments:
#   cid = center id
#
# Returns:
#   Reference to hash of run ids  to run dirnames
#==================================================================
sub GetRuns {
    my ($cid) = @_;
    my %run2dir = ();

    my $sql = "SELECT runid,dirname FROM $opts{runs_table} WHERE centerid=$cid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if ((! $rowsofdata) && $opts{verbose}) {
        warn "$me$mesuffix Found no runs for center '$cid'\n";
        return \%run2dir;
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
        $run2dir{$href->{runid}} = $href->{dirname};
    }
    return \%run2dir;
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
    if ($opts{verbose} > 1) { warn "DEBUG: SQL=$sql\n"; }
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

topmed_failures.pl - Find runs that had a failure

=head1 SYNOPSIS

  topmed_failures.pl look

=head1 DESCRIPTION

Use program as a crontab job to find runs that have failed
(usually cancelled by SLURM).
When one is found the database will be correct to show failure
and if we are lucky the task can be requeued to be tried
again (e.g. with a larger memory)

=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-help>

Generates this output.

=item B<-memory nG>

Force the sbatch --mem setting when submitting a job.

=item B<-noprompt>

Do not prompt to submit a job.

=item B<-partition name>

Force the sbatch --partition setting when submitting a job.

=item B<-qos name>

Force the sbatch --qos setting when submitting a job.

=item B<-realm NAME>

Specifies the database realm to read data from. This defaults to B<topmed>;

=item B<-runs NAME[,NAME,...]>

Specifies a specific set of runs on which to run the action,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-verbose N>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<look4failures>

Directs this program to look for runs that have failed in some step.
These are detected, summary information about it is generated
and the task is marked as failed in the database.

=item B<resubmit>

Directs this program to attempt to resubmit failed jobs.
You are prompted for each submit.

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

