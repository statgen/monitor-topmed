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
    topmedgcecpbcf => "$topmedbin/topmed_gcecpbcf.sh",
    topmedawscopy => "$topmedbin/topmed_awscopy.sh",
    topmedfix => "$topmedbin/topmed_fix.sh",
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
    submitlog => '/tmp/batchsubmit.log',
    consoledir => '/net/topmed/working/topmed-output',
    dryrun => 0,
    verbose => 0,
    maxjobs => 100,
    jobcount => 0,              # Not actually an option, but stats
    jobsnotpermitted => 0,
    jobsfailedsubmission => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help verbose topdir=s center=s runs=s piname=s studyname=s maxjobs=i random
    dryrun suberr datayear=i build=i nopermit
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] arrive|verify|qplot|cram|gcebackup|qplot|gcepush|gcepull|bcf|gcecopy|gcecpbcf|awscopy|fix\n" .
        "Find runs which need some action and queue a request to do it.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift(@ARGV);

my $dbh = DBConnect($opts{realm});
my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

if ($opts{nopermit}) { $ENV{IGNORE_PERMIT} = 1; }   # Stop topmedpermit.pl

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
            #   Check ownership of run. Things break if now owned by topmed
            my $owner = getpwuid($stats[4]);
            my $gid = $stats[5];
            if ($owner ne 'topmed' || $gid ne '2307982') {
                print "$nowdate Ignoring run '$runsref->{$runid}' owned by $owner/$gid\n";
            }
            #   See if we should mark this as arrived
            #   All BAMs must start with NWD. Only skip this if the BAM name
            #   is NWD and it is marked as completed
            if ($href->{bamname} =~ /^NWD/ && $href->{state_arrive} == $COMPLETED) { next; }
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
if ($fcn eq 'gcebackup' || $fcn eq 'backup') {
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
#   Special hook using state_fix so we know when a sample has run the new aplot for Tom
#--------------------------------------------------------------
if ($fcn eq 'qplot') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_gcebackup,state_qplot FROM $opts{bamfiles_table}",
        "WHERE state_qplot!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only do qplot if backup finished
        if ($href->{state_gcebackup} != $COMPLETED) { next; }
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
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_cram,state_gce38push FROM $opts{bamfiles_table}",
        "WHERE state_gce38push!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only send data if cram was done
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
if ($fcn eq 'gcepull' || $fcn eq 'pull') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_gce38push,state_gce38pull FROM $opts{bamfiles_table}", 
        "WHERE state_gce38pull!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only get data if remap was done and requested
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
#   Copy local cram data to GCE storage
#--------------------------------------------------------------
if ($fcn eq 'gcecopy') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_b38,state_gce38bcf,state_gce38copy FROM $opts{bamfiles_table}",
        "WHERE state_gce38copy!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($href->{state_b38} != $COMPLETED) { next; }
        if ($href->{state_gce38bcf} != $COMPLETED) { next; }
        if ($href->{state_gce38copy} == $COMPLETED) { next; }
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
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_b38,state_gce38bcf,state_gce38cpbcf FROM $opts{bamfiles_table}",
        "WHERE state_gce38cpbcf!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($href->{state_b38} != $COMPLETED) { next; }
        if ($href->{state_gce38bcf} != $COMPLETED) { next; }
        if ($href->{state_gce38cpbcf} == $COMPLETED) { next; }
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
#   Copy local data to AWS storage
#--------------------------------------------------------------
if ($fcn eq 'awscopy') {
    #   Get list of all samples yet to process
    my $sql = BuildSQL("SELECT bamid,state_b38,state_gce38bcf,state_aws38copy FROM $opts{bamfiles_table}",
        "WHERE state_aws38copy!=$COMPLETED AND datayear!=3 AND send2aws='Y'");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($href->{state_b38} != $COMPLETED) { next; }
        if ($href->{state_gce38bcf} != $COMPLETED) { next; }
        if ($href->{state_aws38copy} == $COMPLETED) { next; }
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
    my $sql = BuildSQL("SELECT bamid,state_gcebackup,state_qplot FROM $opts{bamfiles_table}",
        "WHERE state_fix!=$COMPLETED");
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { exit; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        #   Only do qplot if backup finished
        if ($href->{state_gcebackup} != $COMPLETED) { next; }
        if ($opts{suberr} && $href->{state_qplot} >= $FAILEDCHECKSUM) {
            $href->{state_qplot} = $REQUESTED;
        }
        if ($href->{state_qplot} == $STARTED || $href->{state_qplot} == $SUBMITTED) { next; }
        if (! BatchSubmit("$opts{topmedqplot} -submit $href->{bamid}")) { last; }
    }
    ShowSummary($fcn);
    exit;

    #   Get list of all samples yet to process
    #my $s = '';
    #if ($opts{center}) { $s = 'b.'; }   # Deal with quirk BuildSQL
    #my $sql = BuildSQL("SELECT bamid,state_fix FROM $opts{bamfiles_table}",
    #    "WHERE ${s}state_fix!=$COMPLETED");
    #my $sth = DoSQL($sql);
    #my $rowsofdata = $sth->rows();
    #if (! $rowsofdata) { exit; }
    #for (my $i=1; $i<=$rowsofdata; $i++) {
    #    my $href = $sth->fetchrow_hashref;
    #    if ($href->{state_fix} == $COMPLETED) { next; }
    #    if ($opts{suberr} && $href->{state_fix} >= $FAILEDCHECKSUM) {
    #        $href->{state_fix} = $REQUESTED;
    #    }
    #    if ($href->{state_fix} != $REQUESTED) { next; }
    #    if (! BatchSubmit("$opts{topmedfix} -submit $href->{bamid}")) { last; }
    #}
    #ShowSummary($fcn);
    #exit;
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
    my $s = $sql;

    #   For a specific center
    if ($opts{center}) {
        $s = $sql . " as b JOIN $opts{runs_table} AS r on b.runid=r.runid " .
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
        $s = $sql . " as b JOIN $opts{runs_table} AS r on b.runid=r.runid " .
            "WHERE r.dirname='$opts{runs}'";
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
    if ($opts{datayear}) { $s .= " AND datayear=$opts{datayear}"; }

    #   Add support for build
    if ($opts{build}) { $s .= " AND build=$opts{build}"; }

    #   Add support for piname
    if ($opts{piname}) { $s .= " AND piname='$opts{piname}'"; }

    #   Add support for piname
    if ($opts{studyname}) { $s .= " AND studyname='$opts{studyname}'"; }

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

=item B<arrive | verify | qplot | cram | gcebackup | qplot | gcepush | gcepull | bcf | gcecopy | gcecpbcf | fix\n" .
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

