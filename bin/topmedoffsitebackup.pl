#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
###################################################################
#
# Name: topmedoffsitebackup.pl
#
# Description:
#   Use this program to figure out if the backup directory link
#   for a run that is backed up offsite is correct.
#   If the incoming cram files have been backed up offsite, then
#   the backup directory link should point to the the incoming files.
#   This is needed so that when we send data to NCBI, we can find
#   the crams to be sent.
#
# ChangeLog:
#   $Log: topmedoffsitebackup.pl,v $
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
use TopMed_Get;
use My_DB;
use Getopt::Long;
use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    topmedcmd => '/usr/cluster/monitor/bin/topmedcmd.pl',
    backupprefix => '../../../../../../',
    centers_table => 'centers',
    bamfiles_table => 'bamfiles',
    runs_table => 'runs',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] check\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;
if ($fcn ne 'check') { die "$Script - Unknown directive '$fcn'. Try $Script --help\n"; }

DBConnect($opts{realm});

#--------------------------------------------------------------
#   Look for runs without a proper backup directory
#--------------------------------------------------------------
my %count = ();
#   Get all the known centers in the database
my $centersref = GetCenters();
foreach my $cid (keys %{$centersref}) {
    my $center = $centersref->{$cid};
    my $runsref = GetRuns($cid) || next;
    #   For each run, see if the backup directory is set up
    foreach my $runid (keys %{$runsref}) {
        my $dirname = $runsref->{$runid};
        if ($dirname =~ /slot/) { next; }   # For development
        #   Backed up offsite?
        my $sql = "SELECT offsite FROM " . $opts{runs_table} .
            " WHERE runid=$runid";
        my $sth = DoSQL($sql);
        my $href = $sth->fetchrow_hashref;
        my $offsite = $href->{offsite};
        $sql = "SELECT bamid,bamname FROM " . $opts{bamfiles_table} . 
            " WHERE runid=$runid LIMIT 1";
        $sth = DoSQL($sql);
        if (! $sth) { next; }           # No files for this run yet
        $href = $sth->fetchrow_hashref;
        if (! $href) { next; }          # No files for this run yet

        #   Assumes all files for a run are the same type
        my $original_file = $href->{bamname};
        my $bamid = $href->{bamid};
        my $incomingdir = `$opts{topmedcmd} wherepath $bamid bam`;
        chomp($incomingdir);
        my $backupsdir = `$opts{topmedcmd} wherepath $bamid backup`;
        chomp($backupsdir);

        #   If original file was a bam, just check that backups can be done
        if ($original_file =~ /\.bam$/) {
            if (! -d $backupsdir) {
                print "You must create a local backup directory for '$center/$dirname'\n\n";
                $count{must_create_local_backup}++;
                next;
            }
            $count{local_backup_exists}++;
            next;
        }

        #   If original file was a cram, be sure backup points to incoming
        if ($original_file =~ /\.cram$/) {
            #   Crams backed up to local storage
            if ($offsite eq 'N') {
                if (-d $backupsdir) { next; }   # Local directory exists, we are fine
                my $d = substr($backupsdir,4);          # Remove /net
                $d = $opts{backupprefix} . $d;
                print "Create a local backup directory for '$center/$dirname' like this:\n" .
                    "  cd " . dirname($backupsdir) . "   # As user TOPMED.  Maybe create this somewhere\n" .
                    "  ln -s $d $dirname\n\n";
                next;
            }
            #   Run is a cram to be backed up offsite
            if ($incomingdir eq $backupsdir) { $count{offsite_backup_exists}++; next; }
            if (-d $backupsdir) {                       # Local directory exists, remove this
                print "Local backup directory for '$center/$dirname' should be removed.\n" .
                    "  Correct using something like this -- BE VERY CAREFUL!\n" .
                    "  cd $backupsdir      # As user TOPMED\n" .
                    "  ls -l      # Be sure these should be removed\n" .
                    "  rm -rf $backupsdir\n\n";
                    $count{must_remove_local_backup}++;
            }
            my $d = substr($incomingdir,4);         # Remove /net
            $d = $opts{backupprefix} . $d;
            print "You must create a backup directory for '$center/$dirname'\n" .
                "This backup should point to the incoming directory since backup files are offsite:\n" .                            
                "  cd " . dirname($backupsdir) . "    # As user TOPMED\n" .
                "  ln -s $d $dirname\n\n";
                $count{must_create_offsite_backup}++;
                next;
        }
        print "$Script - Unable to determine type of files in '$center/$dirname'\n";
    }
}

#   Generate summary
print "Summary of check for backups for our runs:\n";
foreach my $k (sort keys %count) {
    $_ = $k;
    $_ =~ s/_/ /g;
    print "  $_ = $count{$k}\n";
}

exit;



#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedcmd.pl - Update the database for NHLBI TopMed

=head1 SYNOPSIS

  topmedcmd.pl mark 33 arrived completed   # BAM has arrived
  topmedcmd.pl mark NWD00234  arrived completed   # Same BAM has arrived
  topmedcmd.pl unmark 33 arrived           # Reset BAM has arrived

  topmedcmd.pl set 33 jobidqplot 123445    # Set jobidqplot in bamfiles
  topmedcmd.pl set NWD123433 jobidqplot 123445    # Set jobidqplot in bamfiles
  topmedcmd.pl set 2016apr20 offsite N     # Set 2016apr20 in runs

  topmedcmd.pl show 2199 state_cram        # Show a column
  topmedcmd.pl show 2199 center            # Show center 
  topmedcmd.pl show NWD00234 run           # Show run
  topmedcmd.pl show 2016apr20 offsite      # Show offsite for run
  topmedcmd.pl show NWD00234 yaml          # Show everything known about a bamid
 
  topmedcmd.pl wherepath 2199 bam          # Returns real path to bam
  topmedcmd.pl wherepath 2199 backup       # Returns path to backups directory
  topmedcmd.pl wherepath 2199 qcresults    # Returns path to directory for qc.results
  topmedcmd.pl wherepath 2199 console      # Returns path to directory for SLURM output
  topmedcmd.pl wherepath NWD00234 b37      # Returns path to remapped b37 directory and to file
  topmedcmd.pl wherepath 2199 b38          # Returns path to remapped b38 directory and to file

  topmedcmd.pl whathost 2199 bam           # Returns host for bam
  topmedcmd.pl whathost 2199 backup        # Returns host for backups directory
  topmedcmd.pl whathost 2199 qcresults     # Returns host for directory for qc.results

  topmedcmd.pl wherefile 2199 bam          # Returns path to bam file (may not exist)
  topmedcmd.pl wherefile 2199 backup       # Returns path for backups file (may not exist)
  topmedcmd.pl wherefile 2199 qcresults    # Returns path for qc.results *.vb.SelfSM file (may not exist)

  topmedcmd.pl permit add bai braod 2015oct18   # Stop bai job submissions for a run
  topmedcmd.pl permit remove 12             # Remove a permit control
  topmedcmd.pl permit test bai 4567         # Test if we should submit a bai job for one bam

  topmedcmd.pl -maxlongjobs 35 squeue       # Show what is running
=head1 DESCRIPTION

This program supports simple commands to set key elements of the NHLBI database.
The queue of tasks is kept in a MySQL database.
See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-maxlongjobs N>

Specifies how many of the longest running jobs to show.
This defaults to B<10>.

=item B<-persist>

Specifies that SQL connection errors will not fail, but will be retried many times
before finally failing.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program are grouped into several groups which are used
to deal with specific sets of information in the monitor databases.

B<mark bamid|nwdid dirname  [verb] [state]>
Use this to set the state for a particular BAM file.
You may specify the bamid or the NWDID.
Mark will set a date for the process (e.g. arrived sets state_arrive)
and unmark will set that entry to NULL.
The list of verbs and states can be seen by B<perldoc topmedcmd.pl>.

B<permit enable/disable operation center run>
Use this to control the database which allows one enable or disable topmed operations
(e.g. backup, verify etc) for a center or run.
Use B<all> for all centers or all runs or all operations.

B<permit test operation bamid>
Use this to test if an operation (e.g. backup, verify etc) may be submitted 
for a particular bam.

B<set bamid|nwdid|dirname columnname value>
Use this to set the value for a column for a particular BAM file
or run.

B<send2ncbi filelist>
Use this to copy data to NCBI with ascp.

B<show bamid|nwdid|dirname colname|center|run|yaml>
Use this to show information about a particular bamid (or expt_sampleid)
or run name.
Use 'yaml' to display everything known about the bam of interest.

B<unmark bamid|nwdid [verb]>
Use this to reset the state for a particular BAM file to the default
database value.

B<whatnwdid bamid|nwdid>
Use this to get some details for a particular bam.

B<wherepath bamid|nwdid bam|backup|qcresults|console|b37|b38>
If B<bam> was specified, display the path to the real bam file.

If B<backup> was specified, display the path to the backup directory.

If B<qcresults> was specified, display the path to the directory where
the qc.results for this bamid will be.

If B<console> was specified, display the path to the directory where
the SLURM console output.

If B<b37> was specified, display the path to the directory of remapped data for build 37 results can be found.

If B<b38> was specified, display the path to the directory of remapped data for build 37 results can be found.

B<whathhost bamid|nwdid bam|backup|qcresults>
returns the host for the bam, backup or qc.results for the bam file.

B<wherefile bamid|nwdid bam|backup|qcresults>
returns the path to the file for the bam, backup or qc.results. This file may not exist


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

