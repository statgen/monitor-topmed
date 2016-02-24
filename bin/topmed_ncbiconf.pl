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
    summarydir => 'ncbisummaries',
    topdir  => '/net/topmed/incoming/topmed',
    studystatusurl => 'http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetSampleStatus.cgi?study_id=phs000954.v1.p1&rettype=xml',
    xmlfilesdir => '/net/topmed/working/topmed-output/XMLfiles',
    ascpfiles => '/net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri ' .
        '-Q -l 200m -k 1 $ascphost:outgoing/Files',
    ascpinfiles => 'outgoing/Files',
    studystatus => 'latest_samples.xml',
    bamsstatus => 'latest_loaded_files.txt.gz',
    days => 0,
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s fetchfiles xmlfilesdir=s days=i
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] [-fetchfiles] updatedb\n" .
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
    if ($rc) {                              # Pretty lost, get all files and use last loaded file
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
        system('rm -rf Files');         # Clean up cruft
    }
    if ($rc) { die "$Script Unable to fetch files of loaded files at NCBI\n"; }
}

#--------------------------------------------------------------
#   Update state_expt fields in database when NCBI is happy with what we sent
#--------------------------------------------------------------
if ($fcn eq 'ignore') { exit; }     # Convenient way to just load files

if ($fcn eq 'updatedb') {
    my $centersref = GetCenters();
    CheckEXPT($centersref, "$opts{xmlfilesdir}/$opts{studystatus}");
    CheckEXPT($centersref, "$opts{xmlfilesdir}/$opts{bamsstatus}");
    CheckBAMS($centersref, "$opts{xmlfilesdir}/$opts{bamsstatus}");
    exit;
}

die "Invalid request '$fcn'. Try '$Script --help'\n";

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

=item B<-dryrun>

Do not submit any jobs, just show the command to be executed.

=item B<-help>

Generates this output.

=item B<-maxjobs N>

Do not submit more than N jobs for this invocation.
The default for B<-maxjobs> is B<100>.

=item B<-maxsubmits N>

Do not send more than N XML submits to NCBI.
The default for B<-maxjobs> is B<10>.

=item B<-memory nG>

Force the sbatch --mem setting when submitting a job.
This requires the SHELL scripts to handle the environment variable.

=item B<-partition name>

Force the sbatch --partition setting when submitting a job.
This requires the SHELL scripts to handle the environment variable.

=item B<-qos name>

Force the sbatch --qos setting when submitting a job.
This requires the SHELL scripts to handle the environment variable.

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

=item B<arrive | verify | backup | bai | qplot | cram | expt | ncbiorig>

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

