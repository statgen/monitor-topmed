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
  qq(/usr/cluster/topmed/lib/perl5),
  qq(/usr/cluster/topmed/local/lib/perl5),
);
use My_DB;
use Getopt::Long;

use Cwd qw(realpath abs_path);
use Fcntl ':flock';
use POSIX qw(strftime tmpnam);

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

my $topmedbin = '/usr/cluster/topmed/bin';
our %opts = (
    topmedcmd => "$topmedbin/topmedcmd.pl",
    topmedarrive => "$topmedbin/topmed_arrive.sh",
    topmedverify => "$topmedbin/topmed_verify.sh",
    topmedbackup => "$topmedbin/topmed_gcebackup.sh",
    topmedcram   => "$topmedbin/topmed_cram.sh",
    topmedqplot  => "$topmedbin/topmed_qplot.sh",
    topmedexpt  => "$topmedbin/topmed_ncbiexpt.sh",
    topmedncbiorig => "$topmedbin/topmed_ncbiorig.sh",
    topmedncbib37 => "$topmedbin/topmed_ncbib37.sh",
    topmedncbib38 => "$topmedbin/topmed_ncbib38.sh",
    topmedgce38push => "$topmedbin/topmed_gcepush.sh",
    topmedgce38pull => "$topmedbin/topmed_gcepull.sh",
    topmedgce38post => "$topmedbin/topmed_gcepost.sh",
    topmedgce38bcfpush => "$topmedbin/topmed_gcebcfpush.sh",
    topmedgce38bcfpull => "$topmedbin/topmed_gcebcfpull.sh",
    topmedgcecopy => "$topmedbin/topmed_gcecopy.sh",
    topmedbcf  => "$topmedbin/topmed_bcf.sh",
    topmedxml    => "$topmedbin/topmed_xml.pl",
    netdir => '/net/topmed',
    incomingdir => 'incoming/topmed',
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    studies_table => 'studies',
    bamfiles_table => 'bamfiles',
    topdir => '/net/topmed/incoming/topmed',
    backupdir => '/working/backups/incoming/topmed',
    resultsdir => '/incoming/qc.results',    
    dryrun => 0,
    verbose => 0,
    maxjobs => 100,
    jobcount => 0,              # Not actually an option, but stats
    jobsnotpermitted => 0,
    jobsfailedsubmission => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose topdir=s center=s runs=s maxjobs=i random
    dryrun suberr datayear=i
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] arrive|verify|qplot|cram|push|pull|post|pushbcf|pullbcf|gcecopy\n" .
        "Find runs which need some action and queue a request to do it.\n" .
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
            #   See if we should mark this as arrived
            #   All BAMs must start with NWD. Only skip this if the BAM name
            #   is NWD and it is marked as completed
            if ($href->{bamname} =~ /^NWD/ && $href->{state_arrive} == $COMPLETED) { next; }
            #   If the mtime on the file is very recent, it might still be coming
            if ((time() - $stats[9]) < 3600) { next; }  # Catch it next time
            #   Run the command
            BatchSubmit("$opts{topmedarrive} $href->{bamid}");  # Note, not run in SLURM
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
open(my $fh, '>' . $f) || die "Unable to create file '$f': $!\n";
if (! flock($fh, LOCK_EX|LOCK_NB)) { die "Stopping - another instance of '$Script($fcn)' is running\n"; }

#--------------------------------------------------------------
#   Get a list of BAMs that have not been verified
#--------------------------------------------------------------
if ($fcn eq 'verify') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,bamname,state_arrive,state_verify,checksum FROM $opts{bamfiles_table}",
        "WHERE state_verify!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only do verify if file has arrived and ready to be run
        if ($href->{state_arrive} != $COMPLETED) { next; }
        if ($opts{suberr} && $href->{state_verify} >= $FAILEDCHECKSUM) {
            #$href->{state_verify} = $REQUESTED;
            $href->{state_verify} = $REQUESTED;
        }
        if ($href->{state_verify} != $NOTSET && $href->{state_verify} != $REQUESTED) { next; }
        #if ($href->{state_verify} != $NOTSET && $href->{state_verify} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedverify} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have not been converted to CRAMs
#--------------------------------------------------------------
if ($fcn eq 'cram') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_verify,state_cram FROM $opts{bamfiles_table}",
        "WHERE state_cram!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only create cram if file has been verified
        if ($href->{state_verify} != $COMPLETED) { next; }
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
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,bamname,state_cram,state_gcebackup FROM $opts{bamfiles_table}",
        "WHERE state_gcebackup!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only do backup if file has been verified
        if ($href->{state_cram} != $COMPLETED) { next; }
        if ($opts{suberr} && $href->{state_gcebackup} >= $FAILEDCHECKSUM) {
            $href->{state_gcebackup} = $REQUESTED;
        }
        if ($href->{state_gcebackup} != $NOTSET && $href->{state_gcebackup} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedbackup} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Run QPLOT on BAMs
#--------------------------------------------------------------
if ($fcn eq 'qplot') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_verify,state_qplot FROM $opts{bamfiles_table}",
        "WHERE state_qplot!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only do qplot if verify finished
        if ($href->{state_verify} != $COMPLETED) { next; }
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
if ($fcn eq 'push') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_cram,state_gce38push,poorquality FROM $opts{bamfiles_table}",
        "WHERE state_gce38push!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only send data if cram was done
        if ($href->{poorquality} ne 'N') { next; }
        if ($href->{state_cram} != $COMPLETED) { next; }
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
if ($fcn eq 'pull') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_gce38push,state_gce38pull,poorquality FROM $opts{bamfiles_table}", 
        "WHERE state_gce38pull!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only get data if remap was done and requested
        if ($href->{poorquality} ne 'N') { next; }
        if ($href->{state_gce38push} != $COMPLETED) { next; }
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
#   Post process data we fetched from Google Cloud
#--------------------------------------------------------------
if ($fcn eq 'post') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_gce38pull,state_gce38post,poorquality FROM $opts{bamfiles_table}",
        "WHERE state_gce38post!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only post process data that we already fetched
        if ($href->{poorquality} ne 'N') { next; }
        if ($href->{state_gce38pull} != $COMPLETED) { next; }
        if ($opts{suberr} && $href->{state_gce38post} >= $FAILED) {
            $href->{state_gce38post} = $REQUESTED;
        }
        if ($href->{state_gce38post} != $NOTSET &&
            $href->{state_gce38post} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedgce38post} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Push BCF data to Google Cloud for processing.  Maybe never used?
#--------------------------------------------------------------
if ($fcn eq 'pushbcf') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_cram,state_gce38bcf_push,poorquality FROM $opts{bamfiles_table}",
        "WHERE state_gce38bcf_push!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only send data if cram was done
        if ($href->{poorquality} ne 'N') { next; }
        if ($href->{state_cram} != $COMPLETED) { next; }
        if ($opts{suberr} && $href->{state_gce38push} >= $FAILED) {
            $href->{state_gce38bcf_push} = $REQUESTED;
        }
        if ($href->{state_gce38bcf_push} != $NOTSET &&
            $href->{state_gce38bcf_push} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedgce38bcfpush} -submit $href->{bamid} b38")) { last ; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Pull BCF processed data from Google Cloud
#--------------------------------------------------------------
if ($fcn eq 'pullbcf') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_gce38bcf_push,state_gce38bcf_pull,poorquality FROM $opts{bamfiles_table}",
        "WHERE state_gce38bcf_pull!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only get data if bcf was done and requested
        if ($href->{poorquality} ne 'N') { next; }
        if ($href->{state_gce38bcf_push} != $COMPLETED) { next; }
        if ($href->{state_gce38bcf_pull} != $COMPLETED) { next; }
        if ($opts{suberr} && $href->{state_gce38bcf_pull} >= $FAILED) {
            $href->{state_gce38bcf_pull} = $REQUESTED;
        }
        if ($href->{state_gce38pull} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{state_gce38bcf_pull} -submit $href->{bamid} b38")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Create BCF file locally at the moment
#--------------------------------------------------------------
if ($fcn eq 'bcf') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_b38,state_gce38bcf FROM $opts{bamfiles_table}",
        "WHERE state_gce38bcf!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($href->{state_b38} != $COMPLETED) { next; }
        if ($href->{state_gce38bcf} == $COMPLETED) { next; }
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
#   Copy local data to GCE storage   NOT CORRECT
#--------------------------------------------------------------
if ($fcn eq 'gcecopy') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_b38,state_gce38bcf,state_38cp2gce,poorquality FROM $opts{bamfiles_table}",
        "WHERE state_38cp2gce!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($href->{poorquality} ne 'N') { next; }
        if (! $href->{state_b38}) { next; }
        if (! $href->{state_gce38bcf}) { next; }
        if ($href->{state_b38} != $COMPLETED) { next; }
        if ($href->{state_gce38bcf} != $COMPLETED) { next; }
        if ($href->{state_38cp2gce} != $COMPLETED) { next; }
        if ($opts{suberr} && $href->{state_38cp2gce} >= $FAILEDCHECKSUM) {
            $href->{state_38cp2gce} = $REQUESTED;
        }
        if ($href->{state_38cp2gce} != $NOTSET && $href->{state_38cp2gce} != $REQUESTED) { next; }
        if (! BatchSubmit("$opts{topmedgcecopy} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Create experiment XML and sent it to NCBI for this NWDID
#--------------------------------------------------------------
if ($fcn eq 'sexpt') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT * FROM $opts{bamfiles_table}", "WHERE datayear=1");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only do if this NWDID if all the local steps have completed
        if ($href->{nwdid_known} ne 'Y') { next; }
        if ($href->{poorquality} ne 'N') { next; }
        if ($href->{state_qplot} != $COMPLETED) { next; }
        if ($href->{state_cram} != $COMPLETED) { next; }

        #   Check important fields for this BAM are possibly correct
        my $skip = '';
        foreach my $col (qw(expt_sampleid phs library_name nominal_length nominal_sdev base_coord)) {
            if (exists($href->{$col}) && $href->{$col}) { next; }
            $skip .= "$col ";
        }
        #   Do better checking than just is the field non-blank
        if ($href->{expt_sampleid} !~ /^NWD/) {
            $skip .= ' Invalid_NWDID_' . $href->{expt_sampleid};
        }
        if ($skip) {
            print "  BAM '$href->{bamname}' [$href->{bamid}] is ignored because of incomplete data for: $skip\n";
            next;
        }

        if ($opts{suberr} && $href->{state_ncbiexpt} >= $FAILEDCHECKSUM) { $href->{state_ncbiexpt} = $REQUESTED; }
        if ($href->{state_ncbiexpt} != $NOTSET && $href->{state_ncbiexpt} != $REQUESTED) { next; }
        #   Tell NCBI about this NWDID experiment
        BatchSubmit("$opts{topmedexpt} -submit $href->{bamid}");
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of secondary BAMs to be sent to NCBI
#--------------------------------------------------------------
if ($fcn eq 'sorig') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT * FROM $opts{bamfiles_table}",
        "WHERE datayear=1 AND state_ncbiexpt=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only send the original BAM if the experiment was accepted at NCBI
        if ($href->{poorquality} ne 'N') { next; }
        if ($href->{state_ncbiexpt} != $COMPLETED) { next; }
        #   Check important fields for this BAM are possibly correct
        my $skip = '';
        foreach my $col (qw(checksum expt_sampleid)) {
            if (exists($href->{$col}) && $href->{$col}) { next; }
            if ($opts{verbose}) { print "  No value for '$col'\n"; }
            $skip .= "$col ";
        }
        if ($skip) {
            print "  BAM '$href->{bamname}' [$href->{bamid}] is ignored because of incomplete data for: $skip\n";
            next;
        }
        if ($opts{suberr} && $href->{state_ncbiorig} >= $FAILEDCHECKSUM) { $href->{state_ncbiorig} = $REQUESTED; }
        if ($href->{state_ncbiorig} != $NOTSET && $href->{state_ncbiorig} != $REQUESTED) { next; }
        #   Send the secondary BAM to NCBI
        BatchSubmit("$opts{topmedncbiorig} -submit $href->{bamid}");
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of remapped primary BAMs to be sent to NCBI
#--------------------------------------------------------------
if ($fcn eq 'sb37') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT * FROM $opts{bamfiles_table}",
        "WHERE datayear=1 AND state_ncbiexpt=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($href->{poorquality} ne 'N') { next; }
        my $skip = '';
        foreach my $col (qw(checksum expt_sampleid)) {
            if (exists($href->{$col}) && $href->{$col}) { next; }
            if ($opts{verbose}) { print "  No value for '$col'\n"; }
            $skip .= "$col ";
        }
        if ($skip) {
            print "  BAM '$href->{bamname}' [$href->{bamid}] is ignored because of incomplete data for: $skip\n";
            next;
        }
        if ($opts{suberr} && $href->{state_ncbib37} >= $FAILEDCHECKSUM) { $href->{state_ncbib37} = $REQUESTED; }
        if ($href->{state_ncbib37} != $NOTSET && $href->{state_ncbib37} != $REQUESTED) { next; }
        #   Send the remapped CRAM to NCBI
        BatchSubmit("$opts{topmedncbib37} -submit $href->{bamid}");
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
#   boolean if jobs may still be submitted
#==================================================================
sub BatchSubmit {
    my ($cmd) = @_;
    $opts{maxjobs}--;
    if ($opts{maxjobs} < 0) { return 0; }
    if ($opts{maxjobs} == 0 && $opts{verbose}) { print "Limit of jobs to be submitted has been reached\n"; }
    if ($opts{dryrun}) { print "dryrun => $cmd\n"; return 1; }
    my $rc = system("$cmd 2>&1");
    $rc = $rc >> 8;
    if ($rc == 0) {
        $opts{jobcount}++;
        if ($opts{verbose}) { print "submitted => $cmd\n"; }
        return 1;
    }
    $opts{maxjobs}++;                   # Submit failed, keep trying
    if ($rc == 4) { $opts{jobsnotpermitted}++; }
    else { $opts{jobsfailedsubmission}++; }
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
    print "$nowdate $type: $s\n";
}

#==================================================================
# Subroutine:
#   ExecSQL($sql, $die)
#
#   Execute SQL.  Keep trying if connection lost.
#
#   Returns handle for SQL
#==================================================================
sub ExecSQL {
    my ($sql, $die) = @_;
    if ($opts{persist}) { return PersistDoSQL($opts{realm}, $sql); }
    return DoSQL($sql, $die);
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
        $where =~ s/where //i;
    }
    #   For a specific run (overrides center)
    if ($opts{runs}) {
        $s = $sql . " as b JOIN $opts{runs_table} AS r on b.runid=r.runid " .
            "WHERE r.dirname='$opts{runs}'";
        $where =~ s/where //i;
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
        $where .= " AND dirname='$opts{runs}'";
    }
    $sql .= $where;
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Run does not exist   $where\n"; }
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

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-datayear N>

Submit only jobs for samples in a specific year.

=item B<-dryrun>

Do not submit any jobs, just show the command to be executed.

=item B<-help>

Generates this output.

=item B<-maxjobs N>

Do not submit more than N jobs for this invocation.
The default for B<-maxjobs> is B<100>.

=item B<-random>

Randomly select data to be processed. This may not be used with B<-center> or B<-runs>. 
This is intended for cases where a large set of data is to be selected
and we want it to run over a wide set of hosts.
In practice this is only useful for B<push>, B<pull>, B<post>, B<pushbcf> and B<pullbcf>.

=item B<-realm NAME>

Specifies the database realm to read data from. This defaults to B<topmed>;

=item B<-runs NAME[,NAME,...]>

Specifies a specific set of runs on which to run the action,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.

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

=item B<arrive | verify | qplot | cram | push | pull | post | pushbcf | pullbcf | gcecopy>

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

