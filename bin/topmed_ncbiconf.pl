#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
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
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/perl5";
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

my $topmedbin = '/usr/cluster/monitor/bin';
my $ascphost = 'asp-um-sph@gap-submit.ncbi.nlm.nih.gov';
our %opts = (
    topmedcmd => "$topmedbin/topmedcmd.pl",
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
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
    days => 0,
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s fetchfiles loadsummary xmlfilesdir=s days=i all ncbistatusfile=s
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
    my $ym = strftime('%Y%m', localtime);
    my $today = strftime('%d', localtime);
    my $d = sprintf('%02d', $today);
    my $report = 'NCBI_SRA_Files_Full_UM-SPH_';
    my $rc = 1;
    #   Just for consistency's sake, the first of the month is different
    if ($d eq '01') {
        my $f = $report . $ym . $d . '.gz';
        $cmd = "$opts{ascpcmd} $ascphost:$opts{ascpinfiles}/$f .";
        $rc = system($cmd);
        if ($rc == 0) {                     # This day's file exists
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
                rename($f, $opts{bamsstatus});  # Save as lastest file
                print "Fetched file $f, saved as '$opts{xmlfilesdir}/$opts{bamsstatus}'\n";
                last;
            }
        }
    }
    if ($opts{all} || $rc) {                              # Pretty lost, get all files and use last loaded file
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
        @files = sort @files;
        my $f = pop(@files);               # Get last entry
        rename("Files/$f", $opts{bamsstatus}) ||
            die "$Script Unable to find any listing of loaded files. Totally lost: $!\n";
        print "Fetched file 'Files/$f', saved as '$opts{xmlfilesdir}/$opts{bamsstatus}'\n";
        if ($opts{all}) { print "All files were downloaded to output/XMLfiles/Files\n"; }
        else { system('rm -rf Files'); }    # Clean up cruft
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

if ($fcn eq 'updatedb_unused') {
    my $centersref = GetCenters();
    CheckEXPT($centersref, "$opts{xmlfilesdir}/$opts{studystatus}");
    CheckEXPT($centersref, "$opts{xmlfilesdir}/$opts{bamsstatus}");
    if ($opts{ncbistatusfile}) { CheckBAMS($centersref, $opts{ncbistatusfile}); }
    else { CheckBAMS($centersref, "$opts{xmlfilesdir}/$opts{bamsstatus}"); }
    exit;
}

die "Invalid request '$fcn'. Try '$Script --help'\n";

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
    $sql = "SELECT bamid,bamname,expt_sampleid from $opts{bamfiles_table} " .
        "WHERE state_ncbiorig=$DELIVERED";
    LookFor($sql, 'orig', 'src.bam src.cram', \%stats);
    print SummarizeStats(\%stats) . "\n";



    #   Check status for each b37 BAM/CRAM that we sent
    print "$nowdate Checking for remapped build 37 files:   "; 
    %stats = ();
    $sql = "SELECT bamid,bamname,expt_sampleid from $opts{bamfiles_table} " .
        "WHERE state_ncbib37=$DELIVERED";
    #   Last three of extensions are from my bugs
    LookFor($sql, 'b37', 'recal.37.cram recal.37.bam recal.cram recal.bam', \%stats);
    print SummarizeStats(\%stats) . "\n";

    #   Check status for each b38 BAM/CRAM that we sent
    print "$nowdate Checking for remapped build 38 files:   "; 
    %stats = ();
    $sql = "SELECT bamid,bamname,expt_sampleid from $opts{bamfiles_table} " .
        "WHERE state_ncbib38=$DELIVERED";
    #   Last three of extensions are from my bugs
    LookFor($sql, 'b38', 'recal.38.cram', \%stats);
    print SummarizeStats(\%stats) . "\n";

    return;
}

#==================================================================
# Subroutine:
#   LookFor - Search database for files of a certain name
#
# Arguments:
#   initialsql - SQL select in bamfiles for entries of interest
#   type - expt orig b37 b38
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
    print $type . "delivered=$rowsofdata  ";
    my @exts = split(' ', $extensions);

    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        my $bamid = $href->{bamid};
        my $nwdid = $href->{expt_sampleid};
        my @searchfiles = ();
        foreach (@exts) { push @searchfiles,$nwdid . '.' . $_; }
        $href = GetSummaryStatus(@searchfiles);
        if (! $href) {
            #if ($opts{verbose}) { print "$nwdid not found in summary database\n"; }
            next;
        }
        #   BAM/CRAM is now loaded
        if ($href->{file_status} eq 'loaded') {
            $statsref->{$type . 'loaded'}++;
            if ($opts{verbose}) { print "$nwdid $type loaded $href->{upload_date}\n"; }
            my $sql = "UPDATE $opts{bamfiles_table} SET state_ncbi$type=$COMPLETED " .
                "WHERE bamid=$bamid";
print $sql . "\n";
            my $sth = DoSQL($sql);
            next;
        }
        #   BAM/CRAM is received (no state change)
        if ($href->{file_status} eq 'received') {
            $statsref->{$type . 'received'}++;
            if ($opts{verbose}) { print "$nwdid $type received $href->{upload_date}\n"; }
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
    my $sql = "SELECT file_name,file_status,file_error,upload_date from $opts{summary_table} WHERE";
    foreach my $f (@_) {
        my $s = $sql .= " file_name='$f'";
        my $sth = DoSQL($sql);
        if (! $sth) { next; }
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
#   CheckBAMS - Check for BAMs that have completed
#
# Arguments:
#   cref - reference to hash of center data
#   file - log summary file from NCBI
#==================================================================
sub CheckBAMS {
    my ($cref, $file) = @_;

    #   Find all NWDIDnnnnn.*.* in our summary file, read backwards (most recent first)
    #   protected 3563129 2014-09-06T09:24:12 NWDnnnnn.src.bam 10844214270 93ed94e9918b868d0ecd7009c3e427e8 ... loaded BAM etc
    #   protected 5370520 2016-01-31T13:17:51 NWDnnnnnn ... error RUN_XML - some_error_message
    #   protected 292535 2016-01-07T09:48:08 NWDnnnnn.src.bam 1610612736 cf4bbff70440c6dc648b886f7b57ebbe ... replaced_by_5302750 BAM

    my $received = 0;                           # Count number of fiels received
    my $errors = 0;                             # Count number of errors
    my $checksumerrors = 0;                     # Count number of mismatched checksums
    my %processedfile = ();                     # Keep track of all files we look seriously at
    my %receivedfile = ();                      # Keep track of all received files
    my %completed = ();                         # Keep track of number of loaded files by type
    $completed{orig} = $completed{b37} = $completed{b38} = 0;
    my $lastdate = '';

    my $rc;
    #   Read the summary file in reverse order. This allows us to get the LAST state for a file
    #   Doing it forward, we detect errors and all sorts of cruft and THEN need to recognize loaded
    my $cmd = "tac $file";
    if ($file =~ /\.gz$/) { $cmd = "gunzip -c $file | tac - 2> /dev/null |"; }
    $rc = open(IN, $cmd);
    if (! $rc) { die "$Script - Unable to read data from command: $cmd\n"; }
    while (<IN>) {
        chomp();
        #   [0]=protected [2]=date [3]=NWDnnnnnn [5]=checksum
        #   [9]=error|loaded|received [10]=BAM|CRAM|UKNOWN [12]=message
        my @cols = split(' ',$_);
        if ($cols[0] ne 'protected') { next; }
        if ($cols[3] !~ /(NWD\d+)\.(.+)/) { next; } # Break into NWDnnnnnn and suffix
        my ($nwd, $suffix) = ($1, $2);
if ($nwd eq 'NWD378142') {
print "$nwd found\n";
}

        if ($suffix eq 'expt.xml') { next; }
        my $date = substr($cols[2],0,10);           # Date of this action
        if ($opts{days} && $lastdate && $date ne $lastdate) {
            $opts{days}--;
            if ($opts{days} == 0) { last; }     # We've processed all we we supposed to
        }
        if ($opts{verbose} && $date ne $lastdate) { print "Processing data for date=$date\n"; }
        $lastdate = $date;

        if ($cols[9] eq 'loaded') {

            my $sql = "UPDATE $opts{bamfiles_table} SET state_ncbi";
            my $sql2 = "=$COMPLETED WHERE expt_sampleid='$nwd'";
            $processedfile{$cols[3]} = 1;               # Note we have something for this file
            #   Note that cols[10] can be UNKNOWN because of a bug at NCBI 

            if (($suffix eq 'src.bam' && $cols[10] eq 'BAM') ||
                ($suffix eq 'src.cram' && ($cols[10] eq 'CRAM' || $cols[10] eq 'UNKNOWN'))) {
                if ($opts{verbose}) { print "  $nwd orig marked as COMPLETED [$suffix] [$date]\n"; }
                $sql  .= 'orig' . $sql2;
                DoSQL($sql);
                $completed{orig}++;               
                next;
            }

            if (($suffix eq 'recal.bam' && $cols[10] eq 'BAM') ||   # From a bug of mine
                ($suffix eq 'recal.cram' && ($cols[10] eq 'CRAM' || $cols[10] eq 'UNKNOWN')) || # My bug
                ($suffix eq 'recal.37.bam' && ($cols[10] eq 'BAM' || $cols[10] eq 'UNKNOWN')) || # My bug
                ($suffix eq 'recal.37.cram' && ($cols[10] eq 'CRAM' || $cols[10] eq 'UNKNOWN'))) {
                if ($opts{verbose}) { print "  $nwd b37 marked as COMPLETED [$suffix] [$date]\n"; }
                $sql  .= 'b37' . $sql2;
                DoSQL($sql);
                $completed{b37}++;               
                next;
            }

            if (($suffix eq 'recal.38.cram' && ($cols[10] eq 'CRAM' || $cols[10] eq 'UNKNOWN'))) {
                if ($opts{verbose}) { print "  $nwd b38 marked as COMPLETED [$suffix] [$date]\n"; }
                $sql  .= 'b38' . $sql2;
                DoSQL($sql);
                $completed{b38}++;               
                next;
            }

            #   We did not know what was loaded
            print $cols[3] . ": Unrecognized type of file was loaded\n";
            next;
        }

        if ($cols[9] eq 'error') {
            if (exists($processedfile{$cols[3]})) {      # If this has been loaded, ignore it
                if ($opts{verbose}) {
                    print "Error ignored because '$cols[3]' was already seen.  MSG=$cols[12]\n";
                    next;
                }
            }
            $errors++;
            $processedfile{$cols[3]} = 1;               # Note we have something for this file
            $cols[12] =~ s/_/ /g;
            print $cols[3] . ': ' . $cols[12] . "\n";
            next;
        }

        #   A file might be received or replaced, but because ASCP can deliver a file
        #   which does not match its input, we must make sure the checksum matches
        #   This happens far more often than you'd expect
        if ($cols[9] eq 'received' || substr($cols[9], 0, 8) eq 'replace') {
            if (exists($processedfile{$cols[3]})) { next; }
            if ($suffix eq 'expt.xml') { next; }        # Already processed these
            $processedfile{$cols[3]} = 1;           # Note we noticed this file

            my $statecol = GetStateCol($suffix);        # Figure out state column name
            #   Now figure out the checksum column
            my $dbcol = 'unknown';              # Figure out which db column to compare 
            if ($suffix eq 'src.bam')       { $dbcol = 'checksum'; }
            if ($suffix eq 'src.cram')      { $dbcol = 'cramchecksum'; }
            if ($suffix eq 'recal.37.cram') { $dbcol = 'b37bamchecksum'; }
            if ($suffix eq 'recal.38.cram') { $dbcol = 'b38bamchecksum'; }
            if ($dbcol eq 'unknown') {          # I've messed up here
                if (exists($IGNORESUFFIXES{$suffix})) { next; }
                print $cols[3] . ": Incorrectly named file was sent\n";
                next;
            }
            #   We know what this is, now get it's checksum value and state 
            my $sql = "SELECT bamid,$dbcol";
            if ($statecol) { $sql .= ',' . $statecol; }
            $sql .= " FROM $opts{bamfiles_table} WHERE expt_sampleid='$nwd'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }        # How can this happen?  Ignore it
            my $href = $sth->fetchrow_hashref;
            #   Checksum matches, so nothing more needs to be considered
            if ($cols[5] eq $href->{$dbcol}) {
                $received++;
                if ($opts{verbose}) { print "  $cols[3] marked as received, but not loaded yet [$date]\n"; }
                next;
            }
            #   NCBI checksum is wrong. If we think it was delivered, mark as failed
            if ($statecol || $href->{$statecol} != $DELIVERED) { next; }
            $checksumerrors++;
            $errors++;
            print $cols[3] . ": NCBI checksum is incorrect, $href->{bamid} $statecol forced to FAILED\n";
            DoSQL("UPDATE $opts{bamfiles_table} SET $statecol=$FAILED WHERE bamid=$href->{bamid}");
            next;
        }
    }
    close(IN);
#
    if ($received) {
        print "$nowdate  $received BAMs/CRAMs received with correct checksum, waiting to be loaded\n";
    }
    if ($errors) {
        print "$nowdate  $errors new BAMs marked as FAILED ($checksumerrors checksum errors)\n";
    }
    print "$nowdate  Completed $completed{orig}/$completed{b37}/$completed{b38}  [orig/b37/b38]\n";
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
    if ($str =~ /remap.38/)   { $statecol = 'state_ncbib38'; }
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

=item B<-batchsize N>

When sending sets of files to NCBI, batch them into sets of this size.
This might possibly make it easier to keep track of what has been sent.
The default for B<-batchsize> is B<120>.

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-days N>

Specifies the number of days of data from the summary log to parse.
The summary log is read backwards, so B<-days 3> would process only
entries for the last three days.
This only applies when checking the status BAMs/CRAMs that have been sent.

=item B<-fetchfiles>

Causes the summary log file to be fetched from NCBI.

=item B<-help>

Generates this output.

=item B<-ncbistatusfile FILE>

Use this to force processing on an older summary file from NCBI.

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

