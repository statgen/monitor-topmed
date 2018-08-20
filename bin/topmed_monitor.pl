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

if (! -d "/usr/cluster/$ENV{PROJECT}") { die "$Script - Environment variable PROJECT '$ENV{PROJECT}' incorrect\n"; }
our %opts = (
    realm => "/usr/cluster/$ENV{PROJECT}/etc/.db_connections/$ENV{PROJECT}",
    topmedcmd => $Bin . "/topmedcmd.pl",
    topmedarrive => $Bin . "/topmed_arrive.sh",
    topmedverify => $Bin . "/topmed_verify.sh",
    topmedbackup => $Bin . "/topmed_backup.sh",
    topmedcram   => $Bin . "/topmed_cram.sh",
    topmedqplot  => $Bin . "/topmed_qplot.sh",
    topmedexpt  => $Bin . "/topmed_ncbiexpt.sh",
    topmedncbiorig => $Bin . "/topmed_ncbiorig.sh",
    topmedncbib37 => $Bin . "/topmed_ncbib37.sh",
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
    studies_table => 'studies',
    bamfiles_table => 'bamfiles',
    topdir => "/net/$ENV{PROJECT}/incoming/$ENV{PROJECT}",
    submitlog => '/tmp/batchsubmit.log',
    consoledir => "/net/$ENV{PROJECT}/working/$ENV{PROJECT}-output",
    dryrun => 0,
    verbose => 0,
    maxjobs => 100,
    jobcount => 0,              # Not actually an option, but stats
    jobsnotpermitted => 0,
    jobsfailedsubmission => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help verbose topdir=s center=s runs=s piname=s studyname=s maxjobs=i random
    dryrun suberr datayear=i build=i nopermit descending
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] arrive|verify|qplot|cram|backup|qplot|gcepush|gcepull|bcf|gcecopy|gcecpbcf|gcecleanup|awscopy|fix\n" .
        "Find runs which need some action and queue a request to do it.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift(@ARGV);

my $dbh = DBConnect($opts{realm});
my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

if ($opts{nopermit}) { $ENV{IGNORE_PERMIT} = 1; }   # Stop topmedpermit.pl

#   User might provide dirname rather than runid
if (exists($opts{runs}) &&  $opts{runs} =~ /[^0-9,]/) {
    my @r = split(',',$opts{runs});
    my $s = "'" . join("','",@r) . "'";
    my $sql = "SELECT runid from $opts{runs_table} WHERE dirname IN ($s)";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) {
        die "$Script -  $nowdate Unknown run in '$opts{runs}'\n";
    }
    @r = ();
    for (1 .. $rowsofdata) {
        my $href = $sth->fetchrow_hashref;
        push @r, $href->{runid};
    }
    if (! @r) { die "$Script - Unknown run in '$opts{runs}'\n"; }
    $s = join(',', @r);
    print "Getting data for runs: $s\n";
    $opts{runs} = $s;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have arrived, but not processed yet
#--------------------------------------------------------------
if ($fcn eq 'arrive') {
    my $runsref = GetUnArrivedRuns() || next;
    #   For each unarrived run, see if there are bamfiles that arrived
    foreach my $runid (keys %{$runsref}) {
        #   Get list of all bams
        my $sql = "SELECT r.dirname,b.bamname,c.centername,b.bamid,b.expt_sampleid,b.state_arrive " .
            "FROM $opts{bamfiles_table} AS b " .
            "JOIN $opts{runs_table} AS r ON b.runid=r.runid " .
            "JOIN $opts{centers_table} AS c ON r.centerid = c.centerid " .
            "WHERE r.runid=$runid AND r.arrived!='Y'";
        my $sth = DoSQL($sql);
        my $rowsofdata = $sth->rows() || next;
        for (my $i=1; $i<=$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref;
            #   Build path to sample manually
            my $f = $opts{topdir} . "/$href->{centername}/$href->{dirname}/" . $href->{bamname};
            my @stats = stat($f);
            if (! @stats) { next; }     # No real data
            #   Check ownership of run. Things break if not owned by real owner
            my $owner = getpwuid($stats[4]);
            if ($owner ne $ENV{PROJECT}) {
                print "$nowdate Ignoring run '$runsref->{$runid}' owned by $owner\n";
            }
            #   See if we should mark this as arrived
            if ($href->{state_arrive} == $COMPLETED) { next; }
            #   If the mtime on the file is very recent, it might still be coming
            if ((time() - $stats[9]) < 3600) { next; }  # Catch it next time
            #   Run the command
            my $rc = BatchSubmit("$opts{topmedarrive} $href->{bamid}");  # Not run in SLURM
            if ($rc > 0) {          # Unable to submit, capture output
                rename($opts{submitlog}, "$opts{consoledir}/$href->{bamid}-arrive.out")
            }
        }
    }
    ShowSummary('Samples arrived');
    exit;
}

#==================================================================
#   Make sure only one of this subcommand is running
#   This lock is released when the program ends by whatever means
#==================================================================
my $f = "/run/lock/topmed.$fcn.lock";
open(my $fh, '>' . $f) || die "$Script - Unable to create file '$f': $!\n";
if (! flock($fh, LOCK_EX|LOCK_NB)) { die "Stopping - another instance of '$Script($fcn)' is running\n"; }

#--------------------------------------------------------------
#   Get a list of BAMs that have not been verified
#--------------------------------------------------------------
if ($fcn eq 'verify') {
    my $sql = BuildSQL("SELECT bamid,bamname,state_arrive,state_verify,checksum",
        "WHERE b.state_arrive=$COMPLETED AND b.state_verify!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_verify} >= $FAILEDCHECKSUM) {
            $href->{state_verify} = $REQUESTED;
        }
        if ($href->{state_verify} != $NOTSET && $href->{state_verify} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedverify} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have not been converted to CRAMs
#--------------------------------------------------------------
if ($fcn eq 'cram') {
    my $sql = BuildSQL("SELECT bamid,state_verify,state_cram",
        "WHERE b.state_verify=$COMPLETED AND b.state_cram!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_cram} >= $FAILEDCHECKSUM) {
            $href->{state_cram} = $REQUESTED;
        }
        if ($href->{state_cram} != $NOTSET && $href->{state_cram} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedcram} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Backup some files offsite
#--------------------------------------------------------------
if ($fcn eq 'backup') {
    my $sql = BuildSQL("SELECT bamid,bamname,state_cram,state_backup",
        "WHERE b.state_cram=$COMPLETED AND b.state_backup!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_backup} >= $FAILEDCHECKSUM) {
            $href->{state_backup} = $REQUESTED;
        }
        if ($href->{state_backup} != $NOTSET && $href->{state_backup} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedbackup} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Run QPLOT on BAMs
#   Special hook using state_fix so we know when a sample has run the new aplot for Tom
#--------------------------------------------------------------
if ($fcn eq 'qplot') {
    my $sql = BuildSQL("SELECT bamid,state_backup,state_qplot",
        "WHERE b.state_backup=$COMPLETED AND b.state_qplot!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_qplot} >= $FAILEDCHECKSUM) {
            $href->{state_qplot} = $REQUESTED;
        }
        if ($href->{state_qplot} != $NOTSET && $href->{state_qplot} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedqplot} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Push data to Google Cloud for processing
#--------------------------------------------------------------
if ($fcn eq 'gcepush' || $fcn eq 'push') {
    my $sql = BuildSQL("SELECT bamid,state_cram,state_gce38push",
        "WHERE b.state_cram=$COMPLETED AND b.state_gce38push!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_gce38push} >= $FAILED) {
            $href->{state_gce38push} = $REQUESTED;
        }
        if ($href->{state_gce38push} != $NOTSET &&
            $href->{state_gce38push} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedgce38push} -submit $href->{bamid}")) { last ; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Pull processed data from Google Cloud
#--------------------------------------------------------------
if ($fcn eq 'gcepull' || $fcn eq 'pull') {
    my $sql = BuildSQL("SELECT bamid,state_gce38push,state_gce38pull", 
        "WHERE b.state_gce38push=$COMPLETED AND b.state_gce38pull!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_gce38pull} >= $FAILED) {
            $href->{state_gce38pull} = $REQUESTED;
        }
        if ($href->{state_gce38pull} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedgce38pull} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Create BCF file locally at the moment
#--------------------------------------------------------------
if ($fcn eq 'bcf') {
    my $sql = BuildSQL("SELECT bamid,state_b38,state_gce38bcf",
        "WHERE b.state_b38=$COMPLETED AND b.state_gce38bcf!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_gce38bcf} >= $FAILEDCHECKSUM) {
            $href->{state_gce38bcf} = $REQUESTED;
        }
        if ($href->{state_gce38bcf} != $NOTSET && $href->{state_gce38bcf} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedbcf} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Copy local cram data to GCE storage
#--------------------------------------------------------------
if ($fcn eq 'gcecopy') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_b38,state_gce38bcf,state_gce38copy",
        "WHERE b.state_b38=$COMPLETED AND b.state_gce38copy!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_gce38copy} >= $FAILEDCHECKSUM) {
            $href->{state_gce38copy} = $REQUESTED;
        }
        if ($href->{state_gce38copy} != $NOTSET && $href->{state_gce38copy} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedgcecopy} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Copy local bcf data to GCE storage
#--------------------------------------------------------------
if ($fcn eq 'gcecpbcf') {
    my $sql = BuildSQL("SELECT bamid,state_b38,state_gce38bcf,state_gce38cpbcf",
        "WHERE b.state_gce38bcf=$COMPLETED AND b.state_gce38cpbcf!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_gce38cpbcf} >= $FAILEDCHECKSUM) {
            $href->{state_gce38cpbcf} = $REQUESTED;
        }
        if ($href->{state_gce38cpbcf} != $NOTSET && $href->{state_gce38cpbcf} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedgcecpbcf} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Cleanup unnecessary backup files
#--------------------------------------------------------------
if ($fcn eq 'gcecleanup') {
    my $sql = BuildSQL("SELECT bamid,state_gcecleanup",
        "WHERE b.state_gce38cpbcf=$COMPLETED AND b.state_gcecleanup!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_gcecleanup} >= $FAILEDCHECKSUM) {
            $href->{state_gcecleanup} = $REQUESTED;
        }
        if ($href->{state_gcecleanup} != $NOTSET && $href->{state_gcecleanup} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedgcecleanup} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Copy local data to AWS storage
#--------------------------------------------------------------
if ($fcn eq 'awscopy') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_b38,state_gce38bcf,state_aws38copy",
        "WHERE b.state_b38=$COMPLETED AND b.state_aws38copy!=$COMPLETED AND " .
        "b.datayear!=3 AND b.send2aws='Y'");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_aws38copy} >= $FAILEDCHECKSUM) {
            $href->{state_aws38copy} = $REQUESTED;
        }
        if ($href->{state_aws38copy} != $NOTSET && $href->{state_aws38copy} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedawscopy} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Fix some screwup
#--------------------------------------------------------------
if ($fcn eq 'fix') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_backup,state_qplot",
        "WHERE b.state_backup=$COMPLETED AND b.state_fix!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_qplot} >= $FAILEDCHECKSUM) {
            $href->{state_qplot} = $REQUESTED;
        }
        if ($href->{state_qplot} == $STARTED || $href->{state_qplot} == $SUBMITTED) { next; }
        if (! BatchSubmit("$opts{topmedqplot} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Create experiment XML and sent it to NCBI for this NWDID
#--------------------------------------------------------------
if ($fcn eq 'sexpt') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT *", 
        "WHERE b.datayear=3 AND b.nwdid_known='Y' AND b.poorquality='N' " .
        "b.state_cram=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_ncbiexpt} >= $FAILEDCHECKSUM) { $href->{state_ncbiexpt} = $REQUESTED; }
        if ($href->{state_ncbiexpt} != $NOTSET && $href->{state_ncbiexpt} != $REQUESTED) { next; }
        #   Tell NCBI about this NWDID experiment
        if (! BatchSubmit("$opts{topmedexpt} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of secondary BAMs to be sent to NCBI
#--------------------------------------------------------------
if ($fcn eq 'sorig') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT *",
        "WHERE b.datayear=1 AND b.state_ncbiexpt=$COMPLETED AND " .
        "b.poorquality='N'");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_ncbiorig} >= $FAILEDCHECKSUM) { $href->{state_ncbiorig} = $REQUESTED; }
        if ($href->{state_ncbiorig} != $NOTSET && $href->{state_ncbiorig} != $REQUESTED) { next; }
        #   Send the secondary BAM to NCBI
        if (! BatchSubmit("$opts{topmedncbiorig} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of remapped primary BAMs to be sent to NCBI
#--------------------------------------------------------------
if ($fcn eq 'sb37') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT *",
        "WHERE b.datayear=1 AND b.state_ncbiexpt=$COMPLETED AND" .
        "b.poorquality='N'");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{suberr} && $href->{state_ncbib37} >= $FAILEDCHECKSUM) { $href->{state_ncbib37} = $REQUESTED; }
        if ($href->{state_ncbib37} != $NOTSET && $href->{state_ncbib37} != $REQUESTED) { next; }
        #   Send the remapped CRAM to NCBI
        if (! BatchSubmit("$opts{topmedncbib37} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

die "Invalid request '$fcn'. Try '$Script --help'\n";

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
        #   Capture the bamid for ShowSummary
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
    if ($opts{jobcount})            { $s .= "$opts{jobcount} jobs submitted"; }
    if ($opts{jobsnotpermitted})    { $s .=  "  $opts{jobsnotpermitted} job submissions not permitted"; }
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
    $sql .= " FROM $opts{bamfiles_table} AS b ";
    my $s = $sql;

    #   For a specific center
    if ($opts{center}) {
        $s = $sql . " JOIN $opts{runs_table} AS r on b.runid=r.runid " .
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
        $s = $sql . " JOIN $opts{runs_table} AS r on b.runid=r.runid " .
            "WHERE r.runid IN ($opts{runs})";
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
    if ($opts{descending}) { $s .= ' ORDER BY b.bamid DESC'; }
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

    my $sql = "SELECT runid,dirname FROM $opts{runs_table}";
    my $where = " WHERE arrived!='Y'";
    #   Maybe want some runs
    if ($opts{runs}) { 
        $where .= " AND runid IN ($opts{runs})";
    }
    $sql .= $where;
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Run does not exist $sql\n"; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        $run2dir{$href->{runid}} = $href->{dirname};
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

=item B<-datayear N>

Submit only jobs for samples in a specific year.

=item B<-descending>

Select the rows based on the bamid in B<descending> order.
The MySQL default is B<ascending> order.

=item B<-dryrun>

Do not submit any jobs, just show the command to be executed.

=item B<-help>

Generates this output.

=item B<-maxjobs N>

Do not submit more than N jobs for this invocation.
The default for B<-maxjobs> is B<100>.

=item B<-nopermit>

Disable topmedpermit.pl check.  Useful to overwhelm SLURM some times.

=item B<-piname NAME>

Specifies a piname for runs on which to run the action,
e.g. B<Ellinor>.
The default is to run against all pinames.

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

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

