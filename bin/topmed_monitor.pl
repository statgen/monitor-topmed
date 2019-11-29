#!/usr/bin/perl
###################################################################
#
# Name: topmed_monitor.pl
#GetRuns
# Description:
#   Use this program to automatically request actions on data
#   from NHLBI TopMed centers.
#   This program is expected to run as a crontab job.
#
# ChangeLog:
#   $Log: topmed_monitor.pl,v $
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
);
use My_DB;
use Getopt::Long;

use Fcntl ':flock';
use POSIX qw(strftime);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my $NOTSET    = 0;            # Not set
my $REQUESTED = 1;            # Task requested
my $SUBMITTED = 2;            # Task submitted to be run
my $STARTED   = 3;            # Task started
my $DELIVERED = 19;           # Data delivered, but not confirmed
my $COMPLETED = 20;           # Task completed successfully
my $CANCELLED = 89;           # Task cancelled
my $FAILEDCHECKSUM = 98;      # Task failed, because checksum at NCBI bad
my $FAILED    = 99;           # Task failed

#   Pre-check options for project
for (my $i=0; $i<=$#ARGV; $i++) {
    if (($ARGV[$i] eq '-p' || $ARGV[$i] eq '-project') && defined($ARGV[$i+1])) {
        $ENV{PROJECT} = $ARGV[$i+1];
        last;
    }
}
if (! -d "/usr/cluster/$ENV{PROJECT}") { die "$Script - Environment variable PROJECT '$ENV{PROJECT}' incorrect\n"; }

our %opts = (
	datatype => 'genome',		# What kind of data we are handling
    realm => "/usr/cluster/$ENV{PROJECT}/etc/.db_connections/$ENV{PROJECT}",
    topmedcmd => $Bin . "/topmedcmd.pl",
    topmedarrive => $Bin . "/topmed_arrive.sh",
    topmedverify => $Bin . "/topmed_verify.sh",
    topmedbackup => $Bin . "/topmed_backup.sh",
    topmedcram   => $Bin . "/topmed_cram.sh",
    topmedqplot  => $Bin . "/topmed_qplot.sh",
    topmedexpt  => $Bin . "/topmed_ncbiexpt.sh",
    #topmedncbiorig => $Bin . "/topmed_ncbiorig.sh",
    #topmedncbib37 => $Bin . "/topmed_ncbib37.sh",
    topmedncbib38 => $Bin . "/topmed_ncbib38.sh",
    topmedgce38push => $Bin . "/topmed_gcepush.sh",
    topmedgce38pull => $Bin . "/topmed_gcepull.sh",
    topmedgce38post => $Bin . "/topmed_gcepost.sh",
    topmedgcecopy => $Bin . "/topmed_gcecopy.sh",
    topmedgcecpbcf => $Bin . "/topmed_gcecpbcf.sh",
    topmedgcecleanup => $Bin . "/topmed_gcecleanup.sh",
    topmedawscopy => $Bin . "/topmed_awscopy.sh",
    topmedfix => $Bin . "/topmed_fix.sh",
    topmedbcf  => $Bin . "/topmed_bcf.sh",
    topmedxml    => $Bin . "/topmed_xml.pl",
    centers_table => 'centers',    
    runs_table => 'runs',
    runs_pkey => 'runid',
    samples_table => 'bamfiles',
    samples_pkey => 'bamid',
    topdir => "/net/$ENV{PROJECT}/incoming/$ENV{PROJECT}",
    studies_table => 'studies',
    submitlog => '/run/shm/batchsubmit.log',
    consoledir => "/net/$ENV{PROJECT}/working/$ENV{PROJECT}-output",
    dryrun => 0,
    verbose => 0,
    maxjobs => 100,
    jobcount => 0,              # Not actually an option, but stats
    jobsnotpermitted => 0,
    jobsfailedsubmission => 0,
);
my %validverbs = (              # List of valid directives
    arrive => 1,
    verify => 1,
    qplot => 1,
    cram => 1,
    backup => 1,
    qplot => 1,
    gcepush => 1,
    gcepull => 1,
    bcf => 1,
    gcecopy => 1,
    gcecpbcf => 1,
    gcecleanup => 1,
    awscopy => 1,
    fix => 1,
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose topdir=s center=s runs=s piname=s studyname=s maxjobs=i random
    dryrun suberr datayear=i build=i project=s descending datatype=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $verbs = (join '|', sort keys %validverbs);
    warn "$Script [options] [-datatype genome|rnaseq] $verbs\n" .
        "Find runs which need some action and queue a request to do it.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
if ($opts{datatype} eq 'rnaseq') {
	$opts{samples_table} = 'tx_samples';
	$opts{samples_pkey} = 'txseqid';
	$opts{runs_table} = 'tx_projects';
	$opts{runs_pkey} = 'rnaprojectid';
	$opts{topmedarrive} = $Bin . "/rnaseq_arrive.sh";
    $opts{topmedverify} = $Bin . "/rnaseq_verify.sh";
    $opts{topmedbackup} = $Bin . "/rnaseq_backup.sh";
    $opts{topmedawscopy} = $Bin . "/rnaseq_awscopy.sh";
    $opts{topmedfix} = $Bin . "/rnaseq_fix.sh";
	%validverbs = (
		arrive => 1,
		verify => 1,
		backup => 1,
		qplot => 1,
		awscopy => 1,
		fix => 1,
	);
}

my $dbh = DBConnect($opts{realm});
my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

#   User might provide dirname rather than runid
if (exists($opts{runs}) &&  $opts{runs} =~ /[^0-9,]/) {
    my @r = split(',',$opts{runs});
    my $s = "'" . join("','",@r) . "'";
    my $sql = "SELECT $opts{runs_pkey} from $opts{runs_table} WHERE dirname IN ($s)";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) {
        die "$Script -  $nowdate Unknown run in '$opts{runs}'\n";
    }
    @r = ();
    for (1 .. $rowsofdata) {
        my $href = $sth->fetchrow_hashref;
        push @r, $href->{$opts{runs_pkey}};
    }
    if (! @r) { die "$Script - Unknown run in '$opts{runs}'\n"; }
    $s = join(',', @r);
    print "Getting data for runs: $s\n";
    $opts{runs} = $s;
}

#--------------------------------------------------------------
#   Ready to submit jobs for whatever actions are needed
#--------------------------------------------------------------
my $fcn;
foreach $fcn (@ARGV) {
    #   Make sure only one of this subcommand is running
    #   This lock is released when the program ends by whatever means
    my $f = "/run/lock/topmed.$fcn.lock";
    my $fh;
    if (! open($fh, '>' . $f)) {
        warn "$Script - Unable to create file '$f': $!\n";
        next;
    }
    if (! flock($fh, LOCK_EX|LOCK_NB)) { 
        warn "Skipping - another instance of '$Script($fcn)' is running\n";
        next;
    }

    #   Do whatever is asked of us
    if (! exists($validverbs{$fcn})) { die "Invalid request '$fcn'. Try '$Script --help'\n"; }
    my $callsubmit = "Submit_$fcn()";
    eval $callsubmit;

    flock($fh, LOCK_UN);
}
exit;

#==================================================================
# Subroutine:
#   Submit_arrive() - Mark new file as arrived
#==================================================================
sub Submit_arrive {
    my $runsref = GetUnArrivedRuns() || next;
    #   For each unarrived run, see if there are bamfiles that arrived
    foreach my $runid (keys %{$runsref}) {
        #   Get list of all samples
        my $s = '';
        my $bs = '';
        if ($opts{datatype} eq 'genome') { $bs = ',b.'; $s = $bs . 'bamname'; }
        my $sql = "SELECT r.dirname$s,c.centername,b.$opts{samples_pkey},b.expt_sampleid,b.state_arrive " .
            "FROM $opts{samples_table} AS b " .
            "JOIN $opts{runs_table} AS r ON b.$opts{runs_pkey}=r.$opts{runs_pkey} " .
            "JOIN $opts{centers_table} AS c ON r.centerid = c.centerid " .
            "WHERE r.$opts{runs_pkey}=$runid AND r.arrived!='Y'";
        my $sth = DoSQL($sql);
        my $rowsofdata = $sth->rows() || next;
        for (my $i=1; $i<=$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref;
            #   Build path to sample manually
            if ($opts{datatype} eq 'genome') {	# Check if file exists
            	my $f = $opts{topdir} . "/$href->{centername}/$href->{dirname}/" . $href->{bamname};
            	my @stats = stat($f);
            	if (! @stats) { next; }     # No real data
            	#   Check ownership of run. Things break if not owned by real owner
            	my $owner = getpwuid($stats[4]);
            	if ($owner ne $ENV{PROJECT}) {
                	print "$nowdate Ignoring run '$runsref->{$runid}' owned by $owner\n";
            	}
            	#   If the mtime on the file is very recent, it might still be coming
            	if ((time() - $stats[9]) < 3600) { next; }  # Catch it next time
            }
            #   See if we should mark this as arrived
            if ($href->{state_arrive} == $COMPLETED) { next; }
            #   Run the command
            my $rc = BatchSubmit("$opts{topmedarrive} $href->{$opts{samples_pkey}}");  # Not run in SLURM
            if ($rc > 0) {          # Unable to submit, capture output
                rename($opts{submitlog}, "$opts{consoledir}/$href->{$opts{samples_pkey}}-arrive.out")
            }
        }
    }
    ShowSummary('Samples arrived');
}

#==================================================================
# Subroutine:
#   Submit_generic(name, state1, state2)
#		Submit jobs for steps whose checking all looks alike
#
#  	Parms:
#		name - name of step to run for msgs
#		prevstate - database column for previous step that must be complete
#		state2run - database column for step to be submitted
#==================================================================
sub Submit_generic {
    my ($name, $prevstate, $state2run) = @_;

	my $k = 'topmed' . $name;		# Index to shell step bash to run
    my $sql = BuildSQL("SELECT $opts{samples_pkey},$prevstate,$state2run",
        " WHERE $prevstate=$COMPLETED AND $state2run!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { return; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #	If submitting failed step and step failed, force to requested
        if ($opts{suberr} && $href->{$state2run} >= $FAILEDCHECKSUM) {
            $href->{$state2run} = $REQUESTED;
        }
		#	Deal with special cases that don't use normal process
		if ($name eq 'awscopy') {
        	if ($href->{state_aws38copy} != $NOTSET && 
        		$href->{state_aws38copy} != $REQUESTED) { next; }
		}
        #	Ignore anything not requested
        if ($href->{$state2run} != $NOTSET && $href->{$state2run} != $REQUESTED) { next; }

		#	Execute bash script step
        if (! BatchSubmit("$opts{$k} -submit $href->{$opts{samples_pkey}}")) { last; }
    }
    ShowSummary($name);
}

#==================================================================
# Subroutine:
#   Submit_verify() - Submit jobs to verify input sample
#==================================================================
sub Submit_verify {
	Submit_generic('verify', 'state_arrive', 'state_verify');
}

#==================================================================
# Subroutine:
#   Submit_cram() - Submit jobs to create a cram if necessary
#==================================================================
sub Submit_cram {
	Submit_generic('cram', 'state_verify', 'state_cram');
}

#==================================================================
# Subroutine:
#   Submit_backup() - Submit jobs to backup files locally
#==================================================================
sub Submit_backup {
	my $dependency = 'state_cram';
	if ($opts{datatype} eq 'rnaseq') { $dependency = 'state_verify' };
	Submit_generic('backup', $dependency, 'state_backup');
}

#==================================================================
# Subroutine:
#   Submit_qplot() - Submit jobs to run qplot
#==================================================================
sub Submit_qplot {
	Submit_generic('qplot', 'state_backup', 'state_qplot');
}

#==================================================================
# Subroutine:
#   Submit_gcepush() - Submit jobs to push data to GCE for processing
#==================================================================
sub Submit_gcepush {
	Submit_generic('gcepush', 'state_cram', 'state_gce38push');
}

#==================================================================
# Subroutine:
#   Submit_gcepull() - Submit jobs to pull processed data from GCE
#==================================================================
sub Submit_gcepull {
	Submit_generic('gcepull', 'state_gce38push', 'state_gce38pull');
}

#==================================================================
# Subroutine:
#   Submit_bcf() - Submit jobs to create BCF file locally
#==================================================================
sub Submit_bcf {
	Submit_generic('bcf', 'state_b38', 'state_gce38bcf');
}

#==================================================================
# Subroutine:
#   Submit_gcecopy() - Submit jobs to copy local cram data to GCE storage
#==================================================================
sub Submit_gcecopy {
	Submit_generic('gcecopy', 'state_gce38bcf', 'state_gce38copy');
}

#==================================================================
# Subroutine:
#   Submit_gcecpbcf() - Submit jobs to copy local bcf data to GCE storage
#==================================================================
sub Submit_gcecpbc {
	Submit_generic('gcecpbc', 'state_gce38bcf', 'state_gce38cpbcf');
}

#==================================================================
# Subroutine:
#   Submit_gcecleanup() - Submit jobs to cleanup unnecessary backup files
#==================================================================
sub Submit_gcecleanup {
	Submit_generic('gcecleanup', 'state_gce38copy', 'state_gcecleanup');
}

#==================================================================
# Subroutine:
#   Submit_awscopy() - Submit jobs to copy local data to AWS storage
#==================================================================
sub Submit_awscopy {
	Submit_generic('awscopy', 'state_gce38cpbcf', 'state_aws38copy');
}

#==================================================================
# Subroutine:
#   Submit_fix() - Submit jobs for fix
#==================================================================
sub Submit_fix {
	Submit_generic('fix', 'state_backup', 'state_fix');
}

#==================================================================
# Subroutine:
#   BatchSubmit - Run a command which submits the command to batch
#
# Arguments:
#   cmd - command
#
# Returns:
#   0  - no more jobs can be submitted
#   <0 - submit successful, more jobs can be submitted
#   >0 - submit failed, more jobs can be submitted
#==================================================================
sub BatchSubmit {
    my ($cmd) = @_;
    $opts{maxjobs}--;
    if ($opts{maxjobs} < 0) { return 0; }
    if ($opts{maxjobs} == 0 && $opts{verbose}) { print "Limit of jobs to be submitted has been reached\n"; }
    if ($opts{dryrun}) { print "dryrun => $cmd\n"; return -1; }
    my $rc = system("$cmd 2>&1 >> $opts{submitlog}");
    $rc = $rc >> 8;
    if ($rc == 0) {
        $opts{jobcount}++;
        if ($opts{verbose}) { print "submitted => $cmd\n"; }
        #   Capture the sampleid for ShowSummary
        if ($cmd =~ /\s(\d+)$/) {
            $opts{jobbamid} .= $1 . ' ';
        }
        return -1;
    }
    if ($rc == 4) {
        $opts{jobsnotpermitted}++;
        if ($opts{verbose}) { print "Too many on host => $cmd\n"; }
        return 1;
    }
    if ($rc == 5) {
        $opts{jobsnotpermitted}++;
        $opts{maxjobs} = 0;
        if ($opts{verbose}) { print "Too many on system => $cmd\n"; }
        return 0;
    }
    $opts{jobsfailedsubmission}++;
    print "Submit failed => $cmd\n";
    return 1;
}

#==================================================================
# Subroutine:
#   ShowSummary - Print summary of jobs activity
#
# Arguments:
#   type - type of job, verify, backup etc
#==================================================================
sub ShowSummary {
    my ($type) = @_;

    my $s = '';
    if ($opts{jobcount})             { $s .= "$opts{jobcount} jobs submitted"; }
    if ($opts{jobsnotpermitted})     { $s .=  "  $opts{jobsnotpermitted} job submissions not permitted"; }
    if ($opts{jobsfailedsubmission}) { $s .=  "  $opts{jobsfailedsubmission} job submissions failed"; }
    if (! $s) { return; }
    if ($opts{jobbamid}) { $s .= "  bamids=$opts{jobbamid}"; }
    print "$nowdate $type: $s\n";
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
    $sql .= " FROM $opts{samples_table} AS b ";
    my $s = $sql;

    #   For a specific center
    if ($opts{center}) {
        $s = $sql . " JOIN $opts{runs_table} AS r on b.$opts{runs_pkey}=r.$opts{runs_pkey} " .
            "JOIN $opts{centers_table} AS c on r.centerid=c.centerid " .
            "WHERE c.centername='$opts{center}'";
        if ($opts{datayear}) { $s .= " AND b.datayear=$opts{datayear}"; }
        $opts{datayear} = '';
        if ($opts{build}) { $s .= " AND b.build=$opts{build}"; }
        $opts{build} = '';
        $where =~ s/where //i;          # Remove WHERE from caller
    }
    #   For a specific run (overrides center)
    if ($opts{runs}) {
        $s = $sql . " JOIN $opts{runs_table} AS r on b.$opts{runs_pkey}=r.$opts{runs_pkey} " .
            "WHERE r.$opts{runs_pkey} IN ($opts{runs})";
        if ($opts{datayear}) { $s .= " AND b.datayear=$opts{datayear}"; }
        $opts{datayear} = '';
        if ($opts{build}) { $s .= " AND b.build=$opts{build}"; }
        $opts{build} = '';
        $where =~ s/where //i;          # Remove WHERE from caller
    }

    #   Add in caller's WHERE
    if ($where =~ /\s*WHERE\s/) { $s .= ' ' . $where; }
    else { $s .= " AND $where"; }

    #   Add support for datayear
    if ($opts{datayear}) { $s .= " AND b.datayear=$opts{datayear}"; }

    #   Add support for build
    if ($opts{build}) { $s .= " AND b.build=$opts{build}"; }

    #   Add support for piname
    if ($opts{piname}) { $s .= " AND b.piname='$opts{piname}'"; }

    #   Add support for piname
    if ($opts{studyname}) { $s .= " AND b.studyname='$opts{studyname}'"; }

    #   Support randomization
    if ($opts{random}) { $s .= ' ORDER BY RAND()'; }
    if ($opts{descending}) { $s .= ' ORDER BY b.$opts{samples_pkey} DESC'; }
     if ($opts{verbose}) { print "SQL=$s\n"; }
    return $s;
}

#==================================================================
# Subroutine:
#   GetUnArrivedRuns - Get list of all runs that have not arrived.
#       Uses $opts{runs}
#
# Returns:
#   Reference to hash of run ids to run dirnames
#==================================================================
sub GetUnArrivedRuns {
    my %run2dir = ();

    my $sql = "SELECT $opts{runs_pkey},dirname FROM $opts{runs_table}";
    my $where = " WHERE arrived!='Y'";
    #   Maybe want some runs
    if ($opts{runs}) { 
        $where .= " AND $opts{runs_pkey} IN ($opts{runs})";
    }
    $sql .= $where;
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Run does not exist $sql\n"; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        $run2dir{$href->{$opts{runs_pkey}}} = $href->{dirname};
    }
    return \%run2dir;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_monitor.pl - Find runs that need some action

=head1 SYNOPSIS

  topmed_monitor.pl verify
  topmed_monitor.pl -run 20150604 qplot    # Select only samples from one run
  topmed_monitor.pl -center nygc verify    # Select only samples from a center
  topmed_monitor.pl -maxjobs 5 cram        # Only submit a few jobs
  topmed_monitor.pl -datayear 2 qplot      # Do qplot on year 2 samples
  topmed_monitor.pl -random bcf            # Randomly find samples for bcf

  topmed_monitor.pl -project inpsyght -random bcf   # Force same cmd for a project

  topmed_monitor.pl -center broad arrive cram verify qplot   # Submit in this order

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

=item B<-build N>

Submit only jobs for samples from a specific build.

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-datatype genome|rnaseq>

Specifies this is a particular datatype of data. The default is 'genome'.

=item B<-datayear N>

Submit only jobs for samples in a specific year.

=item B<-descending>

Select the rows based on the sampleid in B<descending> order.
The MySQL default is B<ascending> order.

=item B<-dryrun>

Do not submit any jobs, just show the command to be executed.

=item B<-help>

Generates this output.

=item B<-maxjobs N>

Do not submit more than N jobs for this invocation.
The default for B<-maxjobs> is B<100>.

=item B<-piname NAME>

Specifies a piname for runs on which to run the action,
e.g. B<Ellinor>.
The default is to run against all pinames.

=item B<-project PROJECT>

Specifies these commands are to be used for a specific project.
Warning, this can only be abbreviated as B<-p> or <-project>.
The default is to use the environment variable PROJECT.

=item B<-random>

Randomly select data to be processed. This may not be used with B<-center> or B<-runs>. 
This is intended for cases where a large set of data is to be selected
and we want it to run over a wide set of hosts.

=item B<-runs NAME>

Specifies a run on which to run the action,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-studyname NAME>

Specifies a study name for runs on which to run the action,
e.g. B<MESA>.
The default is to run against all studies.

=item B<-suberr>

Submit the job if the state is B<error>. Normally tasks in this state
are not submitted to be run.

=item B<-topdir PATH>

Specifies the path to where the tree of BAMs exists. This defaults to  B</incoming/topmed>;

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<arrive | verify | qplot | cram | backup | qplot | gcepush | gcepull | bcf | gcecopy | gcecpbcf | gcecleanup | fix\n" .
y>

Directs this program to look for runs that have not been through the process name
you provided and to queue a request they be verified.

=back


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2019 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

