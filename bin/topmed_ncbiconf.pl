#!/usr/bin/perl
###################################################################
#
# Name: topmed_ncbiconf.pl
#
# Description:
#   Use this program to confirm what data was sent to NCBI
#   and that they are happy with it.
#
# ChangeLog:
#   $Log: topmed_ncbiconf.pl,v $
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
use TopMed_Get;
use Getopt::Long;

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

#   Ignore files with these suffixes as they were sent incorrectly
my %IGNORESUFFIXES = (
    'recal.bam' => 1,
    'recal.cram' => 1,
    'recal.37.bam' => 1,
    'final.sqz.bam' => 1,
    'squeezed.bam' => 1,
    'hg19.bam' => 1,
);

my $topmedbin = '/usr/cluster/topmed/bin';
my $ascphost = 'asp-um-sph@gap-submit.ncbi.nlm.nih.gov';
our %opts = (
    topmedcmd => "$topmedbin/topmedcmd.pl",
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    ascpcmd => '/usr/cluster/bin/ascp -i /net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri -Q -l 200m -k 1 -q',
    centers_table => 'centers',
    runs_table => 'runs',
    studies_table => 'studies',
    bamfiles_table => 'bamfiles',
    summary_table => 'ncbi_summary',
    summarydir => 'ncbisummaries',
    topdir  => '/net/topmed/incoming/topmed',
    studystatusurl => 'http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetSampleStatus.cgi?study_id=phs000954.v1.p1&rettype=xml',
    xmlfilesdir => '/net/topmed/working/topmed-output/XMLfiles',
    ascpfiles => '/net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri ' .
        '-Q -l 200m -k 1 $ascphost:outgoing/Files',
    ascpinfiles => 'outgoing/Files',
    studystatus => 'latest_samples.xml',
    bamsstatus => 'latest_loaded_files.txt.gz',
    summarymysqlcmd => '/home/topmed/.mysql/topmed_ncbiconf.cmd',      # Used to load summary file into database
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s fetchfiles loadsummary xmlfilesdir=s all
)) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] [-fetchfiles -loadsummary] updatedb\n" .
        "$Script [options] [-fetchfiles] ignore     # Fetch files only\n" .
        "Confirm files sent to NCBI were error free.\n";
    exit 1;
}
my $fcn = shift(@ARGV);
my $dbh = DBConnect($opts{realm});
my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

#--------------------------------------------------------------
#   Get files of interest
#--------------------------------------------------------------
if ($opts{fetchfiles}) {
    chdir $opts{xmlfilesdir} ||
        die "$Script Unable to CD to '$opts{xmlfilesdir}'\n";
    #   Get XML of all known experiments
    my $cmd = "wget -O $opts{studystatus} -o /dev/null '$opts{studystatusurl}'";
    system($cmd) &&
        die "$Script Failed to fetch sample list to '$opts{studystatus}'\n";
    print "$Script Fetched sample list to '$opts{xmlfilesdir}/$opts{studystatus}'\n";

    #   Get today's list of files sent to NCBI and their status
    #   Unfortunately, if there is no activity yesterday, there is none for today
    #   so we must do some guessing as to the last report generated
    my $rc = 1;
    if (! $opts{all}) {
        my $ym = strftime('%Y%m', localtime);
        my $today = strftime('%d', localtime);
        my $d = sprintf('%02d', $today);
        my $report = 'NCBI_SRA_Files_Full_UM-SPH_';
        #   Just for consistency's sake, the first of the month is different
        if ($d eq '01') {
            my $f = $report . $ym . $d . '.gz';
            $cmd = "$opts{ascpcmd} $ascphost:$opts{ascpinfiles}/$f .";
            $rc = system($cmd);
            if ($rc == 0) {                     # This day's file exists
                SummaryTableMerge($f);
                rename($f, $opts{bamsstatus});  # Save as lastest file
                print "Fetched file $f, saved as '$opts{xmlfilesdir}/$opts{bamsstatus}'\n";
            }
        }
        if ($rc) {                              # Not the first, search for the last file
            $report = 'NCBI_SRA_Files_UM-SPH_';
            foreach (0 .. 4) {                  # Try to get one of the last few days
                my $d = sprintf('%02d', $today - $_);
                my $f = $report . $ym . $d . '.gz';
                $cmd = "$opts{ascpcmd} $ascphost:$opts{ascpinfiles}/$f .";
                $rc = system($cmd);                 # Failure is OK, up to a point
                if ($rc == 0) {                     # This day's file exists
                    SummaryTableMerge($f);
                    rename($f, $opts{bamsstatus});  # Save as lastest file
                    print "Fetched file $f, saved as '$opts{xmlfilesdir}/$opts{bamsstatus}'\n";
                    last;
                }
            }
        }
    }
    if ($opts{all} || $rc) {                    # Pretty lost, get all files 
        # Here is command:  /usr/cluster/bin/ascp -i /net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri \
        #   -Q -l 200m -k 1 -q asp-um-sph@gap-submit.ncbi.nlm.nih.gov:outgoing/Files .
        $cmd = "$opts{ascpcmd} $ascphost:$opts{ascpinfiles} .";   # Creates directory Files
        $rc = system($cmd);
        if ($rc != 0) { 
            die "$Script Unable to find any listing of loaded files. Totally lost: $!\n  CMD=$cmd\n"
        }
        # Got lots of files, find latest of loaded files
        opendir(my $dh, 'Files') ||
            die "$Script Unable to read directory 'Files': $!\n";
        my @files = grep { /\.gz/ } readdir($dh);
        close($dh);

        #   The monthly summaries are named special, rename them to look like the others
        my @files2 = ();
        foreach my $f (@files) {
            $f = 'Files/' . $f;
            if ($f =~ /_Full/) {
                my $f2 = $f;
                $f2 =~ s/_Full//;
                rename($f, $f2) ||
                    die "$Script - cmd failed: 'rename $f $f2': $!\n";
                push @files2,$f2;
            }
            else { push @files2,$f; }
        }

        #   Walk through each file and update the database
        foreach my $f (sort @files2) {
            SummaryTableMerge($f);
        }
        if ($opts{all}) { print "All files were downloaded to output/XMLfiles/Files\n"; }
    }
    if ($rc) { die "$Script Unable to fetch files of loaded files at NCBI\n"; }
}

#--------------------------------------------------------------
#   Load summary file into database
#--------------------------------------------------------------
if ($opts{loadsummary}) {
    chdir $opts{xmlfilesdir} ||
        die "$Script Unable to CD to '$opts{xmlfilesdir}'\n";
    #   Get XML of all known experiments
    my $cmd = "gunzip -c $opts{bamsstatus} > $opts{summary_table}";
    system($cmd) &&
        die "$Script Failed to unzip '$opts{bamsstatus}': CMD=$cmd\n";
    $cmd = "$opts{summarymysqlcmd}";
    my $rc = system($cmd);
    unlink($opts{summary_table});
    if ($rc) { die "$Script Failed to load summary file: CMD=$cmd\n"; }
    print "$Script Loaded summary file into database table '$opts{summary_table}'\n";
}

#--------------------------------------------------------------
#   Update state_expt fields in database when NCBI is happy with what we sent
#--------------------------------------------------------------
if ($fcn eq 'ignore') { exit; }     # Convenient way to just load files

if ($fcn eq 'updatedb') { CheckSummary(); exit; }

die "Invalid request '$fcn'. Try '$Script --help'\n";

#==================================================================
# Subroutine:
#   SummaryTableMerge - Read a file and merge it into ncbi_summary
#==================================================================
sub SummaryTableMerge {
    my ($f) = @_;
    #   Walk through each file and update the database
    #   3 Might started with NWDxxxxxx
    #   9 loaded, error etc
    #   12 Possibly a useful error msg
    my @xdbcols = qw(
        realm
        upload_id
        upload_date
        file_name
        file_size
        file_md5sum
        upload_name
        upload_size
        upload_md5sum
        file_status
        file_type
        load_date
        file_error
        submissions
        loaded_runs
        unloaded_runs
        suppressed_runs
        loaded_analyses
        unloaded_analyses
        suppressed_analyses
    );
    open (INSQL, "gunzip -c $f |") ||
        die "$Script - Unable to open file '$f': $!\n";
    my @dbcols = split("\t", <INSQL>);
    print "Merging data from '$f' ";
    my $upd = 0;
    my $ins = 0;
    my %ERRORS = ();                    # NWDID to hash of path to error msg
    my %LOADED = ();                    # NWDID loaded files of interest to NWDID
    while (<INSQL>) {
        chomp();
        my @cols = split("\t",$_);
        if ($cols[3] !~ /^NWD/) { next; }   # Ignore cruft from other projects
        if ($#cols != $#dbcols) {
            die "$Script - Number columns in $. of '$f' was incorrect " .
                "($#cols != $#dbcols)\n";
        }
        if ($cols[3] !~ /^(NWD\d{6})/) {
            die "$Script - error field detected, but unable to get NWDid, line $. of '$f'. Line=$_\n";
        }
        my $nwd = $1;

        #   Clean up some fields for easier post processing or to fix errors in names
        if ($cols[12] eq '-') { $cols[12] = ''; }
        if ($cols[3] =~ /\.recal\./)  { $cols[3] =~ s/\.recal\./\.remap\./; }
        if ($cols[3] =~ /\.remap\.bam/)  { $cols[3] =~ s/\.remap\.bam/\.remap\.37\.bam/; }
        if ($cols[3] =~ /\.remap\.cram/) { $cols[3] =~ s/\.remap\.cram/\.remap\.37\.cram/; }
        if ($cols[3] =~ /\.sqz.bam/) { next; }  # Ignore old bug in my code

        my $setsql = 'SET ';
        for (my $i=0; $i<=$#dbcols; $i++) {
            $setsql .= $dbcols[$i] . "='" . $cols[$i] . "',";
        }
        chop($setsql);
        my $sql = "UPDATE $opts{summary_table} $setsql WHERE file_name='$cols[3]'";
        my $sth = DoSQL($sql,0);        # This can fail
        if (! $sth->rows()) {           # If failed, then do insert
            $sql = "INSERT INTO $opts{summary_table} $setsql";
            DoSQL($sql);                # This may not fail
            $ins++;
        }
        else { $upd++; }

        #   Keep track of all BAMs or CRAMs that were loaded so we can delete the errors
        if ($cols[9] eq 'loaded' && $cols[3] =~ /am$/) { $LOADED{$cols[3]} = $nwd; }

        #   Catch errors here in %ERRORS    cols[3]=file_name cols[12]=file_error
        if ($cols[9] ne 'error') { next; }
        if ($cols[3] =~ /submit.xml/) { next; }     # Ignore errors on submits
        my %E = ();
        my $eref = $ERRORS{$nwd} || \%E;        # Either previous errors or new one
        if ($cols[12] =~ /^Run_\w+_is_rej(.+)/) { $cols[12] = 'Rej' . $1; }
        $eref->{$cols[3]} .= $cols[12] . "\n";
        $ERRORS{$nwd} = $eref;
    }
    close(INSQL);

    #   Remove any error messages for loaded Files from msgs like this
    #     protected 5399247 2016-02-07T14:27:12 NWD433184-remap.37.run.xml ... error
    #     protected 5399248 2016-02-08T11:56:21 NWD433184.remap.37.cram ... loaded   (cram or bam)
    #     protected 5401844 2016-02-08T12:03:34	NWD433184-secondary.37.run.xml ... error
    #     protected 5468239 2016-02-25T15:43:55	NWD433184.src.bam ... loaded   (cram or bam)
    my $removederrs = 0;
    foreach my $file (keys %LOADED) {
        my $nwd = $LOADED{$file};
        if ($file =~ /$nwd\.src\.\w+am$/) {
            if (delete $ERRORS{$nwd}{"$nwd-secondary.37.run.xml"}) { 
                if ($opts{verbose}) { print "Removed errors for $file\n"; }
                $removederrs++;
            }
            next;
        }
        if ($file =~ /$nwd\.remap\.(3\d)\.\w+am$/) {
            if (delete $ERRORS{$nwd}{"$nwd-remap.$1.run.xml"}) { 
                if ($opts{verbose}) { print "Removed errors for $file\n"; }
                $removederrs++;
            }
            next;
        }
        print "Did not know what to do with '$file'\n";
    }

    #   Update database with all the error messages
    my $nwderrs = 0;
    foreach my $nwdid (keys %ERRORS) {
        if ($opts{verbose}) { print $nwdid . ":\n"; }
        my $eref = $ERRORS{$nwdid};
        foreach my $k (keys %$eref) {
            if ($opts{verbose}) { print "  $k => $eref->{$k}"; }
            $eref->{$k} =~ s/_/ /g;         # Spaces make the msg more readable
            my $sql =  "UPDATE $opts{bamfiles_table} SET emsg='$k => $eref->{$k}' WHERE expt_sampleid='$nwdid'";
            ##DoSQL($sql);      # Emsg is used for all errors
            $nwderrs++;
        }
    }
    print "  ($upd updates, $ins inserts, found $nwderrs NWDID errors, removed $removederrs errors)\n";
}

#==================================================================
# Subroutine:
#   CheckSummary - Check the summary database for changes in status
#==================================================================
sub CheckSummary {
    #my ($cref, $table) = @_;
    my %stats = ();                     # Collect stats on attempts

    #   Check status for each expt we sent
    print "$nowdate Checking for loaded experiments:   "; 
    %stats = ();
    my $sql = "SELECT bamid,bamname,expt_sampleid from $opts{bamfiles_table} " .
        "WHERE state_ncbiexpt=$DELIVERED";
    LookFor($sql, 'expt', 'expt.xml', \%stats);
    print SummarizeStats(\%stats) . "\n";

    #   Check status for each original BAM/CRAM that we sent
    print "$nowdate Checking for loaded original files:   "; 
    %stats = ();
    $sql = "SELECT bamid,bamname,expt_sampleid,cramchecksum,checksum from $opts{bamfiles_table} " .
        "WHERE state_ncbiorig=$DELIVERED";
    LookFor($sql, 'orig', 'src.bam src.cram final.sqz.bam', \%stats);
    print SummarizeStats(\%stats) . "\n";

    #   Check status for each b37 BAM/CRAM that we sent
    print "$nowdate Checking for remapped build 37 files:   "; 
    %stats = ();
    $sql = "SELECT bamid,bamname,expt_sampleid,cramchecksum,checksum from $opts{bamfiles_table} " .
        "WHERE state_ncbib37=$DELIVERED";
    #   Last three of extensions are from my bugs
    LookFor($sql, 'b37', 'recal.37.cram recal.37.bam recal.cram recal.bam', \%stats);
    print SummarizeStats(\%stats) . "\n";

    return;
}

#==================================================================
# Subroutine:
#   LookFor - Search database for files of a certain name
#
# Arguments:
#   initialsql - SQL select in bamfiles for entries of interest
#   type - expt orig b37
#   extensions - string of file extensions to the nwdid
#   statsref - hash reference to array to collect stats in
#
#  Returns:
#   Boolean if successful or not
#==================================================================
sub LookFor {
    my ($initialsql, $type, $extensions, $statsref) = @_;

    #   Get list of all nwdids we are interested in
    my $sth = DoSQL($initialsql);
    my $rowsofdata = $sth->rows();
    if ($rowsofdata <= 0) { return undef(); }
    print $type . " delivered=$rowsofdata  ";
    my @exts = split(' ', $extensions);

    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        my $bamid = $href->{bamid};
        my $nwdid = $href->{expt_sampleid};
        my $cramchecksum = $href->{cramchecksum};
        my $checksum = $href->{checksum} || 'notset';
        my @searchfiles = ();
        foreach (@exts) { push @searchfiles,$nwdid . '.' . $_; }
        $href = GetSummaryStatus(@searchfiles);
        if (! $href) {
            #if ($opts{verbose}) { print "$nwdid not found in summary database\n"; }
            next;
        }
        my $upload_date = $href->{upload_date} || '';
        my $loaded_runs = $href->{loaded_runs} || '';
        #   BAM/CRAM is now loaded - this is algorithm Tom likes (Dec 2016) 
        if ($href->{file_status} eq 'loaded' && 
            substr($upload_date,0,2) eq '20' &&     # Year is 20xx
            ($type eq 'expt' || substr($loaded_runs,0,3) eq 'SRR')) {    # Have SRR
            $statsref->{$type . 'loaded'}++;
            if ($opts{verbose}) { print "$nwdid $type loaded $upload_date\n"; }
            my $sql = "UPDATE $opts{bamfiles_table} SET state_ncbi$type=$COMPLETED " .
#                "time_ncbi$type='" . substr($upload_date,0,19) ."'," .
#                "emsg='' " .
                "WHERE bamid=$bamid";
            my $sth = DoSQL($sql);
            next;
        }
        #   BAM/CRAM is received (no state change)
        if ($href->{file_status} eq 'received') {
            $statsref->{$type . 'received'}++;
            if ($opts{verbose}) { print "$nwdid $type received $href->{upload_date}\n"; }
            #   Received could have an incorrect checksum
            if ($href->{file_name} =~ /cram/) { $checksum = $cramchecksum; }
            if ($type ne 'expt' && $href->{file_md5sum} ne $checksum) {
                $statsref->{$type . 'checksumerror'}++;
                my $statecol = 'state_ncbi' . $type;         # Mark this column in error
                my $sql = "UPDATE $opts{bamfiles_table} SET $statecol=$FAILEDCHECKSUM WHERE bamid=$bamid";
                if ($opts{verbose}) { print "$nwdid: NCBI checksum incorrect, $bamid $statecol forced to FAILEDCHECKSUM\n"; }
                DoSQL($sql);
            }
            next;
        }
        #   BAM/CRAM had an error
        if ($href->{file_status} eq 'error') {
            $statsref->{$type . 'error'}++;
            print "$nwdid $type in error $href->{file_error}\n";
            my $sql = "UPDATE $opts{bamfiles_table} SET state_ncbi$type=$FAILED " .
                "WHERE bamid=$bamid";
            my $sth = DoSQL($sql);
            next;
        }
        #   BAM/CRAM was replaced, not a real error
        if ($href->{file_status} =~ /replaced_by/) {
            next;
        }
        $statsref->{$type . 'unknown'}++;
        print "$nwdid $type unknown status '$href->{file_status}'\n";
    }
    return 1;
}

#==================================================================
# Subroutine:
#   SummarizeStats
#
# Arguments:
#   statsref - hash reference to array to collect stats in
#
#   Returns:
#     Summary string
#==================================================================
sub SummarizeStats {
    my ($statsref) = @_;
    my $s = '';
    foreach my $k (sort keys %$statsref) {
        $s .= "$k=$statsref->{$k} ";
    }
    if (! $s) { return ''; }
    return $s;
}

#==================================================================
# Subroutine:
#   GetSummaryStatus - Get status for NWDID
#
# Arguments:
#   filearray - array of files to search for
#
#   Returns:
#     hash reference to database query results
#==================================================================
sub GetSummaryStatus {
    my $sql = "SELECT file_name,file_status,file_error,upload_date,loaded_runs,file_md5sum from $opts{summary_table} WHERE";
    foreach my $f (@_) {
        my $s = $sql . " file_name='$f'";
        my $sth = DoSQL($s);
        if (! $sth) { next; }
        my $rows = $sth->rows();
        if ($rows <= 0) { next; }
        my $href = $sth->fetchrow_hashref;
        return $href;
    }
    return undef();
}

#==================================================================
# Subroutine:
#   CheckEXPT - Check for experiments that have completed
#
# Arguments:
#   cref - reference to hash of center data
#   file - XML summary file from NCBI
#==================================================================
sub CheckEXPT {
    my ($cref, $file) = @_;

    #   Get hash of all experiments delivered to NCBI
    my $nwd2bamid = GetNWDlist($cref, 'state_ncbiexpt');
    if (! %{$nwd2bamid}) { print "$Script - No experiments have been delivered\n"; return; }

    #   Find all submitted_sample_id="NWD560497"
    #     or
    #   protected 5265301 2015-12-31T12:59:26 NWD481739.expt.xml 1720 72b081e2eac659ff8e38b8aa54e360ca
    #     NWD481739-expt.tar 10240 509ed745d862576639f459fe5f115f2d loaded EXPERIMENT_XML - - SRA324210 ...
    my %ncbinwdids = ();
    my $completed = 0;
    my $rc;
    if ($file =~ /\.gz$/) { $rc = open(IN, "gunzip -c $file |"); }
    else { $rc = open(IN, $file); }
    if (! $rc) { die "$Script - Unable to read file '$file': $!\n"; }
    while (<IN>) {
        if (/submitted_sample_id=.(NWD\d+)/) {
            $ncbinwdids{$1} = 1;
            next;
        }
        if (/protected.+\s+(NWD\d+).expt.xml.+loaded\s+EXPERIMENT_XML\s+-\s+-\s+SRA\d+/) {
            $ncbinwdids{$1} = 1;
            next;
        }        
    }
    close(IN);
    if (! %ncbinwdids) { print "  No completed experiments found at NCBI\n"; return; }

    foreach my $nwdid (keys %{$nwd2bamid}) {
        if (! exists($ncbinwdids{$nwdid})) { next; }
        #   This NWDID is now known
        if ($opts{verbose}) { print "  Completed experiment for $nwdid (bamid=$nwd2bamid->{$nwdid})\n"; }
        DoSQL("UPDATE $opts{bamfiles_table} SET state_ncbiexpt=$COMPLETED WHERE bamid=$nwd2bamid->{$nwdid}");
        $completed++;
    }
    print "$nowdate  Marked $completed experiments as completed\n";
    return;
}

#==================================================================
# Subroutine:
#   GetNWDlist - Get hash of bamfiles where something was delivered
#
# Arguments:
#   cref - reference to hash of center data
#   col - database column in $DELIVERED state
#
# Returns:
#   reference to hash of nwdid to bamid
#==================================================================
sub GetNWDlist {
    my ($cref, $col) = @_;

    #   Get hash of all original BAMs delivered to NCBI
    my %nwd2bamid = ();
    my $countexpt = 0;
    foreach my $cid (keys %{$cref}) {
        my $centername = $cref->{$cid};
        my $runsref = GetRuns($cid) || next;
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all data that was delievered
            my $sql = "SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid' AND $col=$DELIVERED";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                $nwd2bamid{$href->{expt_sampleid}} = $href->{bamid};
                $countexpt++;
            }
        }
    }
    return \%nwd2bamid;
}

#==================================================================
# Subroutine:
#   GetStateCol - Given a string, return the database column associated with it
#
# Arguments:
#   str - string
#
# Returns:
#   name of state database column or ''
#==================================================================
sub GetStateCol {
    my ($str) = @_;
    my $statecol = '';
    if ($str =~ /xml/)        { $statecol = 'xmlerror'; }
    if ($str =~ /recal.37/)   { $statecol = 'state_ncbib37'; }    # My bug
    if ($str =~ /recal.bam/)  { $statecol = 'state_ncbib37'; }    # My bug
    if ($str =~ /recal.cram/) { $statecol = 'state_ncbib37'; }    # My bug
    if ($str =~ /remap.37/)   { $statecol = 'state_ncbib37'; }
    if ($str =~ /secondary/)  { $statecol = 'state_ncbiorig'; }
    if ($str =~ /src.bam/)    { $statecol = 'state_ncbiorig'; }
    if ($str =~ /src.cram/)   { $statecol = 'state_ncbiorig'; }
    return $statecol;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_ncbiconf.pl - Confirm that NCBI is happy with data we sent

=head1 SYNOPSIS

  topmed_ncbiconf.pl -fetchfiles updatedb    # Typical case
  topmed_ncbiconf.pl -fetchfiles ignore      # Just force fetch of files

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

=item B<-all>

If B<-fetchfiles> was specified, this option will force all files to
be fetched, not just the single file we think we want.

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-fetchfiles>

Causes the summary log file to be fetched from NCBI.

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the database realm to read data from. This defaults to B<topmed>;

=item B<-runs NAME[,NAME,...]>

Specifies a specific set of runs on which to run the action,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-xmlfilesdir directory>

Specifies a directory where the summary log files will be downloaded.
This defaults to B</net/topmed/working/topmed-output/XMLfiles>.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<updatedb | ignore>

Controls what steps will be run. The usual command is B<updatedb>, but at times
you might want to just fetch the log files, in which case B<ignore> can be
specified to avoid a misleading error message.

=back


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015- and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

