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
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s fetchfiles xmlfilesdir=s
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
    if ($d eq '01') {                       # The first of the month is different
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

    #   Get hash of all original BAMs delivered to NCBI
    my $nwd2bamid = GetNWDlist($cref, 'state_ncbiorig');
    if (! %{$nwd2bamid}) { print "$Script - No original BAMs have been delivered\n"; return; }

    #   Find all NWDIDnnnnn.src.bam in
    #   protected 3563129 2014-09-06T09:24:12 NWDnnnnn.src.bam 10844214270 93ed94e9918b868d0ecd7009c3e427e8 = = = loaded BAM etc
    #     or
    #   protected NWD792235-remap.37.run.xml ... error RUN_XML - some_error_message
    my %loadednwdids = ();
    my %nwdid2errormsg = ();
    my %completed = ();
    $completed{orig} = $completed{b37} = $completed{b38} = 0;
    my $rc;
    if ($file =~ /\.gz$/) { $rc = open(IN, "gunzip -c $file |"); }
    else { $rc = open(IN, $file); }
    if (! $rc) { die "$Script - Unable to read file '$file': $!\n"; }
    while (<IN>) {
        if (/protected\s+.+\s+(NWD\d+\S*.run.xml)\s+.+\s+error\s+RUN_XML\s+-\s+(\S+)/) {
            my ($x, $msg) = ($1, $2);
            $msg =~ s/_/ /g;
            $nwdid2errormsg{$x} .= $msg . "\n";
            next;
        }
        if (/protected\s+.+\s+(NWD\S+).src.bam\s+.+=\s+=\s+=\s+loaded\sBAM/) {
            $loadednwdids{$1} = 'orig';
            next;
        }
        if (/protected\s+.+\s+(NWD\S+).recal.37.bam\s+.+=\s+=\s+=\s+loaded\sBAM/) { 
            $loadednwdids{$1} = 'b37';
            next;
        }
        if (/protected\s+.+\s+(NWD\S+).recal.38.bam\s+.+=\s+=\s+=\s+loaded\sBAM/) { 
            $loadednwdids{$1} = 'b38';
            next;
        }
    }
    close(IN);

    #   Found list of NWDIDs in error
    #   If this changes the state for the NWDID, update the database
    #   Messages we've seen so far:
    #     protected	...	NWD792235-remap.37.run.xml	...	error_msg
    #     protected	...	NWD539447-secondary.37.run.xml	...	Updates_are_prohibited_for_withdrawn_Run(SRR2173553)
    my $errors = 0;
    foreach my $key (keys %nwdid2errormsg) {
        if ($opts{verbose}) { print "  FAILED: $key - $nwdid2errormsg{$key}"; }
        #   Get actual NWDID
        my $nwdid = $key;
        my $statecol = '';
        if ($key =~ /(NWD\d+)/) { $nwdid = $1; }
        else { die "Unable to isolate NWDID:  NWDID=$nwdid Errormsg=$nwdid2errormsg{$key}\n"; }
        #   Figure out state column name
        if ($key =~ /remap.37/) { $statecol = 'state_ncbib37'; }
        if ($key =~ /remap.38/) { $statecol = 'state_ncbib38'; }
        if ($key =~ /secondary/) { $statecol = 'state_ncbiorig'; }
        if (! $statecol) { die "Unable to figure out database column  NWDID=$nwdid Errormsg=$nwdid2errormsg{$key}\n"; }
        #   Now see if we already knew about this error
        delete($loadednwdids{$nwdid});          # Make sure we do not treat this as loaded
        my $sql = "SELECT $statecol FROM $opts{bamfiles_table} " .
                "WHERE expt_sampleid='$nwdid' AND $statecol!=$FAILED";
        my $sth = DoSQL($sql);
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) { next; }            # Already knew about this error
        DoSQL("UPDATE $opts{bamfiles_table} SET $statecol=$FAILED WHERE expt_sampleid='$nwdid'");
        print "$nowdate  New error detected - failure for '$nwdid' $statecol\n";
        $errors++;
    }
    if ($errors) { print "$nowdate  $errors new BAMs marked as FAILED\n"; }

    #   Found list of loaded NWDIDs
    foreach my $nwdid (keys %loadednwdids) {
        if (! exists($nwd2bamid->{$nwdid})) { next; }
        #   This NWDID is known
        my $v = $loadednwdids{$nwdid};
        if ($opts{verbose}) { print "  Completed $v BAM for $nwdid (bamid=$nwd2bamid->{$nwdid})\n"; }
        DoSQL("UPDATE $opts{bamfiles_table} SET state_ncbi$v=$COMPLETED WHERE bamid=$nwd2bamid->{$nwdid}");
        $completed{$v}++;
    }
    foreach my $k (keys %completed) {
        #if ($completed{$k}) {
            print "$nowdate  $completed{$k} $k BAMs marked as COMPLETED\n";
        #}
    }
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

