#!/usr/bin/perl
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
use Cwd qw(realpath abs_path);
use POSIX qw(strftime);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
#   Pre-check options for project
for (my $i=0; $i<=$#ARGV; $i++) {
    if (($ARGV[$i] eq '-p' || $ARGV[$i] eq '-project') && defined($ARGV[$i+1])) {
        $ENV{PROJECT} = $ARGV[$i+1];
        last;
    }
}
if (! -d "/usr/cluster/$ENV{PROJECT}") { die "$Script - Environment variable PROJECT '$ENV{PROJECT}' incorrect\n"; }

our %opts = (
    topmedcmd => "/usr/cluster/$ENV{PROJECT}/bin/topmedcmd.pl",
    sacct => '/usr/cluster/bin/sacct',
    realm => "/usr/cluster/$ENV{PROJECT}/etc/.db_connections/$ENV{PROJECT}",
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    outdir => "/net/$ENV{PROJECT}/working/$ENV{PROJECT}-output",
    verbose => 0,                       # Must be defined for My_DB
);

my $STARTED   = 3;                      # Task started
my @statecols =                         #   States to be checked
    qw/ state_verify
    state_backup
    state_cram
    state_qplot
    state_gce38push
    state_gce38pull
    state_gce38bcf
    state_gce38copy
    state_gce38cpbcf
    state_gcecleanup
    state_aws38copy
/;

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s fail requeue project=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] show\n" .
        "Find tasks which we think were started, but do not seem to be running.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift(@ARGV);
my $dbh = DBConnect($opts{realm});
my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

#   User might provide runid rather than name of run
if (exists($opts{runs}) &&  $opts{runs} =~ /^\d+$/) {
    my $sql = "SELECT dirname from $opts{runs_table} WHERE runid=$opts{runs}";
    my $sth = DoSQL($sql);
    if ($sth) {
        my $href = $sth->fetchrow_hashref;
        $opts{runs} = $href->{dirname};
    }
}

#--------------------------------------------------------------
#   Get a list of jobs that have been submitted, see if any failed
#--------------------------------------------------------------
if ($fcn =~ /^show/) {
    foreach my $statecol (@statecols) {
        my $sql = BuildSQL("SELECT bamid,bamname,$statecol FROM $opts{bamfiles_table}", 
            "WHERE $statecol=$STARTED");
        my $sth = DoSQL($sql);
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) {
            if ($opts{verbose}) { print "No jobs running for '$statecol'\n"; }
            next;
        }
        #   For each bamid with something started, check each task state
        for (my $i=1; $i<=$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref;
            $opts{"job_$statecol"}++;
            #   Find failure for this sample
            my $s = $statecol;
            $s =~ s/state_//;               # state_cram becomes cram
            if ($s eq 'gce38push') { $s = 'gcepush'; }
            if ($s eq 'gce38pull') { $s = 'gcepull'; }
            if ($s eq 'gce38bcf')  { $s = 'bcf'; }
            if ($s eq 'gce38copy') { $s = 'gcecopy'; }
            Find_Failure($href, $s);
        }
    }
    #   Give summary of what we found out   Stats start with job_
    print $nowdate . ' ';
    foreach my $k (sort keys %opts) {
        if ($k !~ /^job_/) { next; }
        my $kmsg = $k;
        $kmsg =~ s/_/ /g;
        print "  $k = $opts{$k}";
    }
    print "\n";
    exit;
}

die "Invalid request '$fcn'. Try '$Script --help'\n";

#==================================================================
# Subroutine:
#   Find_Failure - Check state of a job that has started
#
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
#   task - name of task that is running, e.g. cram
#
# Returns:
#   nothing
#==================================================================
sub Find_Failure {
    my ($h, $task) = @_;

    #   See if there is a log file laying around for this task. If so, the job is
    #   recent and we can check it. No log file and and started is a bad thing.    
    my $f = $opts{outdir} . "/$h->{bamid}-$task.out";
    if (! open(IN,$f)) {
        $opts{job_no_log_file}++;
        if ($opts{verbose}) { print "Job started but no log found for $task bamid=$h->{bamid}\n"; }
        return;
    }
    #   Read log file, searching for line with slurm Jobid
    my $slurmid;
    while (<IN>) {
        #   Maybe the job was cancelled by SLURM right away
        if (/CANCELLED AT \S+\s+(.+)/) {
            my $msg = "$task job (bamid=$h->{bamid}) never got started, cancelled by SLURM $1\n";
            print $msg;
            $opts{job_started_and_failed}++;
            if ($opts{fail}) {
                my $cmd = "$opts{topmedcmd} mark $h->{bamid} $task failed";
                my $rc = system($cmd);
                if (! $rc) { print "$h->{bamid} $task marked as failed\n"; }
                else { print "Unable to mark $h->{bamid} $task as failed: $cmd\n"; }
                $opts{job_started_and_failed}++;
            }
            else { $opts{job_started_and_failed}++; }   # Just count it as failed
            close(IN);
            return;
        }
        #   Look through log file for SLURMID printed by script
        my @words = split(' ', $_);
        if ((! @words) || $words[0] !~ /^#==/) { next; }
        $slurmid = $words[3];
        last;
    }
    close(IN);
    if (! $slurmid) {
        warn "$Script - Unable to get SLURMID for '$f' - status unknown\n";
        return;
    }

    #   See what SLURM thinks of job $slurmid
    my $cmd = $opts{sacct} . ' -j ' . $slurmid;
    my @lines = split("\n", `$cmd 2>&1`);
    if (! defined($lines[2])) { return; }
    my @words = split(' ', $lines[2]);
    if ($words[0] ne $slurmid) { return; }
    my $msg = "We think $task job is running, but it was $words[5]: " .
        "bamid=$h->{bamid} slurmid=$slurmid status=$words[6]\n";
    if ($words[5] eq 'COMPLETED') {         # Running job succeeded
        print $msg;
        $opts{job_started_and_completed}++;
        return;
    }
    if ($words[5] =~ /CANCELLED|NODE_FAIL|FAILED|TIMEOUT/) {    # Running job failed
        print $msg;
        $opts{job_started_and_failed}++;
        if ($opts{fail}) {
            $cmd = "$opts{topmedcmd} mark $h->{bamid} $task failed";
            my $rc = system($cmd);
            if (! $rc) { print "$h->{bamid} $task marked as failed\n"; }
            else { print "Unable to mark $h->{bamid} $task as failed: $cmd\n"; }
        }
       if ($opts{requeue}) {
            my $cmd = "$opts{topmedcmd} mark $h->{bamid} $task requested";
            my $rc = system($cmd);
            if (! $rc) { print "$h->{bamid} $task marked as requested\n"; }
            else { print "Unable to mark $h->{bamid} $task as requested: $cmd\n"; }
            $opts{job_started_and_requested}++;
        }
        return;
    }
    if ($opts{verbose}) { print "Job $h->{bamid} $h->{bamname} $f is still running\n"; }
    $opts{jobs_still_running}++;
    return;
}

#==================================================================
# Subroutine:
#   BuildSQL - Complete SQL statement based on options
#
# Arguments:
#   sql - initial sql
#   where - optional WHERE clause, without WHERE
#
# Returns
#   Completed SQL statement
#==================================================================
sub BuildSQL {
    my ($sql, $where) = @_;
    my $s = $sql;

    #   For a specific center
    if ($opts{center}) {
        $s = $sql . " as b JOIN $opts{runs_table} AS r on b.runid=r.runid " .
            "JOIN $opts{centers_table} AS c on r.centerid=c.centerid " .
            "WHERE c.centername='$opts{center}'";
        if ($opts{datayear}) { $s .= " AND b.datayear=$opts{datayear}"; }
        $opts{datayear} = '';
        $where =~ s/where //i;          # Remove WHERE from caller
    }
    #   For a specific run (overrides center)
    if ($opts{runs}) {
        $s = $sql . " as b JOIN $opts{runs_table} AS r on b.runid=r.runid " .
            "WHERE r.dirname='$opts{runs}'";
        if ($opts{datayear}) { $s .= " AND b.datayear=$opts{datayear}"; }
        $opts{datayear} = '';
        $where =~ s/where //i;          # Remove WHERE from caller
    }

    #   Add in caller's WHERE
    if ($where =~ /\s*WHERE\s/) { $s .= ' ' . $where; }
    else { $s .= " AND $where"; }

    #   Add support for datayear
    if ($opts{datayear}) { $s .= " AND datayear=$opts{datayear}"; }

    #   Support randomization
    if ($opts{random}) { $s .= ' ORDER BY RAND()'; }
    return $s;
}


#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_failures.pl - Find runs that had a failure

=head1 SYNOPSIS

  topmed_failures.pl show
  topmed_failures.pl -fail show
  topmed_failures.pl -requeue show

=head1 DESCRIPTION

Use program as a crontab job to find runs that have failed
(sometimes cancelled by SLURM).
When one is found the database will be corrected to show a failure
or perhaps even be requeued to try again.


=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-fail>

If a job has failed and we thought it was running, mark it as failed in the database.

=item B<-help>

Generates this output.

=item B<-project PROJECT>

Specifies these commands are to be used for a specific project.
Warning, this can only be abbreviated as B<-p> or <-project>.
The default is to use the environment variable PROJECT.

=item B<-realm NAME>

Specifies the database realm to read data from. This defaults to B<topmed>;

=item B<-requeue>

If a job has failed and we thought it was running, mark it as requested in the database.

=item B<-runs NAME[,NAME,...]>

Specifies a specific set of runs on which to run the action,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.


=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<show>

Directs this program to look for runs that have failed in a non-standard way
(e.g. typically cancelled by SLURM).
You may want to specify B<-mark> so that the database record is updated as failed.

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

