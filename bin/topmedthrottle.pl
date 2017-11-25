#!/usr/bin/env perl
###################################################################
#
# Name: topmedthrottle.pl
#
# Description:
#   Use this program to control the rate at which topmed jobs are submitted to SLURM
#   This is necessary because SLURM cannot properly manage multiple QOS
#   and using only one QOS, we can easily run too many of one job type,
#   overwhelming the filesystem, NFS or the network
#
# ChangeLog:
#   $Log: topmedthrottle.pl,v $
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

use Getopt::Long;
use My_DB;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    runs_table => 'runs',
    conf => "$Bin/../etc/topmedthrottle.conf",
    topmedsummary => '/run/shm/Slurm.summary',
    topmedcluster => '/usr/cluster/topmed/bin/topmedcluster.pl summary',
    topmedmonitor => '/usr/cluster/topmed/bin/topmed_monitor.pl',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose conf=s dry-run action=s
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] submit\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
my %JOBS = ();                  # All the type of jobs we know how to run
if ($fcn eq 'submit') { SubmitJobs(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   SubmitJobs()
#
#   Print summary of topmed-related jobs
#==================================================================
sub SubmitJobs {
    #my ($p1, $p2) = @_;
    my %jobsqueued = ();
    my %countbyhost = ();

    ParseConf($opts{conf});         # Sets JOBS
    if (! %JOBS) { die "$Script - No jobtypes found in '$opts{conf}'\n"; }

    #   Get details about jobs that are running and how many per host
    #   This should already be waiting, but if not get it dynamically
    #   $opts{topmedsummary} might be created about the same time as we start.
    #   If so, wait a little so the file might be updated.
    my @s = stat($opts{topmedsummary});
    my $fileage = time();               # How old is file to read
    if (@s) { $fileage = $fileage - $s[9]; }
    if ($fileage > 295 && $fileage < 303) { sleep(10); }  # Wait for update of file
    if ($fileage > 320) {               # Older than ~5 min
        $opts{topmedsummary} = "$opts{topmedcluster} |";    # Invoke pgm to get info
    }
    if ($opts{verbose}) { print "Reading SLURM summary from $opts{topmedsummary}\n"; }
    my $in;
    open($in, $opts{topmedsummary}) ||
        die "$Script - Unable to get details about jobs from '$opts{topmedsummary}': $!\n";
    while (my $l = <$in>) {         # jobtype queued=N running=N qhosts= h n... rhosts= h n ...
        if ($l =~ /\s*(\S+) queued=(\d+) running=(\d+)/) {
            $jobsqueued{$1} = $2+$3;
        }
        if ($l =~ /qhosts=(.*)\s+rhosts=(.*)/) {
            my $q = $1;
            my $r = $2;
            my @w = split(' ',$q);          # Get all queued 'host count' pairs
            for (my $i=0; $i<$#w; $i=$i+2) { $countbyhost{$w[$i]} += $w[$i+1]; }
            @w = split(' ',$r);          # Get all running 'host count' pairs
            for (my $i=0; $i<$#w; $i=$i+2) { $countbyhost{$w[$i]} += $w[$i+1]; }
        }
    }
    close($in);

    #   %JOBS has list of jobtypes and total number of running and queued
    #   Now go through the limits specified in the conf file and
    #   generate a call to submit jobs if there is the need
    foreach my $j (keys %JOBS) {    # Calculate jobs to be submitted
        if ($opts{action} && $opts{action} ne $j) { next; }
        if (! exists($jobsqueued{$j})) { $jobsqueued{$j} = 0; }   # None of these running
        my $max = 0;                    # Maybe submit this many
        if ($jobsqueued{$j} < $JOBS{$j}{max}) { $max = $JOBS{$j}{max} - $jobsqueued{$j}; }
        if ($max > 0) {                 # Submit more jobs 
            my $cmd = $opts{topmedmonitor} . " -max $max";
            if ($opts{verbose}) { $cmd .= ' -v'; }
            if (exists($JOBS{$j}{'mon-options'})) { $cmd .= ' ' . $JOBS{$j}{'mon-options'}; }
            $cmd .= ' ' . $j;
            if (exists($JOBS{$j}{'mon-log'})) { $cmd .= ' | tee -a ' . $JOBS{$j}{'mon-log'} . ' 2>&1'; }
            if ($opts{verbose}) { print "==>$cmd\n"; }
            if ($opts{'dry-run'}) { print "DRYRUN: $cmd\n"; }               
            else {
                system($cmd) &&
                    die "$Script - Submit failed: $cmd\n";
            }
        }
        else {
            if ($opts{verbose}) { print "May not submit job for '$j'\n"; }
        }
    }

    return;
}

#==================================================================
# Subroutine:
#   ParseConf(file)
#
#   Parse config file setting JOBS
#==================================================================
sub ParseConf {
    my ($file) = @_;
    my $in;
    my $OPENAMP = '{';
    my $CLOSEAMP = '}';
    open($in, $file) ||
        die "$Script - Unable to read file '$file': $!\n";
    while (<$in>) {
        chomp();
        if (/^#/) { next; }
        if (! /\S/) { next; }
        if (/^job\s+(\S+)\s+$OPENAMP/i) {
            my $job = $1;
            while (<$in>) {
                if (/\s*$CLOSEAMP/) { last; }       # End of job
                if (/\s*(\S+)\s*=\s*(.+)\s*$/) {    # Capture key = value
                    $JOBS{$job}{$1} = $2;
                }
            }
            next;
        }
        print "Did not recognize: $_\n";
    }
    close($in);
    if ($opts{verbose}) {
        my @k = keys %JOBS;
        print "Defined " . scalar(@k) . " jobs:\n";
        foreach my $j (@k) {
            print "  $j ";
            foreach my $key (keys %{$JOBS{$j}}) { print "$key=$JOBS{$j}{$key} "; }
            print "\n";
        }
        print "\n";
    }
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedthrottle.pl - Control the rate at which topmed jobs are submitted


=head1 SYNOPSIS

  topmedcluster.pl submit


=head1 DESCRIPTION

Use this program to control the rate at which topmed jobs are submitted to SLURM

This is necessary because SLURM cannot properly manage multiple QOS
and using only one QOS, we can easily run too many of one job type,
overwhelming the filesystem, NFS or the network


=head1 OPTIONS

=over 4

=item B<-action NAME>

Use this to restrict any jobs to be submitted to only one action, e.g. cram, gcecopy etc.

=item B<-conf PATH>

Specifies the path to a configuration file

=item B<-dry-run>

Prevents jobs from actually being submitted, but shows what would be done.

=item B<-help>

Generates this output.

=item B<-verbose>

Provided for developers to see additional information.

=back


=head1 PARAMETERS

Parameters to this program are used
to deal with specific sets of information in the monitor databases.

B<submit>
Use this to submit jobs to SLURM.


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

