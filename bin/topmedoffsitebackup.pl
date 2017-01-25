#!/usr/bin/perl
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
use lib (
  qq($FindBin::Bin),
  qq($FindBin::Bin/../lib),
  qq($FindBin::Bin/../lib/perl5),
  qq($FindBin::Bin/../local/lib/perl5),
  qq(/usr/cluster/topmed/lib/perl5),
  qq(/usr/cluster/topmed/local/lib/perl5),
);

use TopMed_Get;
use My_DB;
use Getopt::Long;
use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    topmedcmd => '/usr/cluster/topmed/bin/topmedcmd.pl',
    topmedpath => '/usr/cluster/topmed/bin/topmedpath.pl',
    backupprefix => '../../../../../../',
    centers_table => 'centers',
    bamfiles_table => 'bamfiles',
    runs_table => 'runs',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose
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
        my $incomingdir = `$opts{topmedpath} wherepath $bamid bam`;
        chomp($incomingdir);
        my $backupsdir = `$opts{topmedpath} wherepath $bamid backup`;
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
            if (! -l $backupsdir) {                 # Symlink exists, mark to be removed
                my $d = $backupsdir . '.removeme';
                if ((! -e $d) && (-d $backupsdir)) {
                    system("mv $backupsdir $d");
                    print "Local backup directory for '$center/$dirname' renamed to be removed.\n";
                    $count{must_remove_local_backup}++;
                }
            }
            #   Create the symlink for the backup
            my $cmd = "ln -s $opts{backupprefix}" . 'topmed/incoming/topmed/' .
                    "$center/$dirname /net/topmed/working/backups/incoming/topmed/$center/$dirname";
            if (system($cmd)) {
                print "$Script - Unable to create backup symlink:  $cmd\n";
                $count{must_create_offsite_backup}++;
            }
            else {
                print "$Script - Created symlink for backup for '$center/$dirname'\n";
                $count{created_offsite_backup_link}++;
            }
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

topmedoffsitebackup.pl - See if the backup directory link for a run is correct

=head1 SYNOPSIS

  topmedoffsitebackup.pl check

=head1 DESCRIPTION

Use this program to figure out if the backup directory link
for a run that is backed up offsite is correct.
If the incoming cram files have been backed up offsite, then
the backup directory link should point to the the incoming files.
This is needed so that when we send data to NCBI, we can find
the crams to be sent.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS


B<check>
Use this to look for crams that are not in the correct place.

=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

