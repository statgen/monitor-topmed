#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
###################################################################
#
# Name: topmedrename.pl bamid checksum bamfilepath
#
# Description:
#   Use this program to rename a BAM to it's nwdid name
#   and to correct the MD5 file entry for the BAM
#
# ChangeLog:
#   $Log: topmedrename.pl,v $
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
use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my $topmedbin = '/usr/cluster/monitor/bin';
our %opts = (
    topmedcmd => "$topmedbin/topmedcmd.pl",
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 2 || $opts{help}) {
    warn "$Script [options] bamid checksum bamfilepath\n" .
        "Rename the BAM file to its NWDID name, correct the MD5 file entry.\n";
    exit 1;
}
my $bamid = shift(@ARGV);
my $checksum = shift(@ARGV);
my $bamfilepath = shift(@ARGV);
my $dbh = DBConnect($opts{realm});

#--------------------------------------------------------------
#   Rename the BAM file, change the MD5 file
#--------------------------------------------------------------
my $sth = DoSQL("SELECT * FROM $opts{bamfiles_table} WHERE bamid='$bamid'");
my $rowsofdata = $sth->rows();
if (! $rowsofdata) { die "Script BAM '$bamid' does not exist in database [$bamfilepath]\n"; }
my $href = $sth->fetchrow_hashref;
my $nwdid = $href->{expt_sampleid};
if ($nwdid !~ /^NWD/) { die "Script BAM '$bamid' NWDID [$nwdid] was not set [$bamfilepath]\n"; }
#   Checksum hack - allow this to be 'useoldvalue' and then use the checksum value from the db
if ($checksum eq 'useoldvalue' && $href->{checksum}) { $checksum =$href->{checksum}; }

#   CD to place where BAM exists
my $dirname = dirname($bamfilepath);
chdir($dirname) ||
    die "$Script Unable to CD to '$dirname': $!\n";
my $bamfile = basename($bamfilepath);
if ($bamfile =~ /^NWD/) {
    if ($opts{verbose}) { print "$Script Bamid=$bamid already uses NWD format\n"; }
    exit;
}

#   Rename the BAM file
my $newbamfile = $bamfile;
if ($newbamfile !~ /^[^.]+\.(.+)/) { die "$Script Unable to parse '$newbamfile'\n"; }
$newbamfile = $nwdid . '.' . $1;
rename($bamfile, $newbamfile) ||
    die "$Script Unable to rename $bamfile $newbamfile (bamid=$bamid)\n";
if ($opts{verbose}) { print "$Script Renamed $bamfile to $newbamfile\n"; }

#   Save original bamname once in database, change name in database
if ((! defined($href->{bamname_orig})) || $href->{bamname_orig} eq '') {
    DoSQL("UPDATE $opts{bamfiles_table} SET bamname_orig='$bamfile' WHERE bamid=$bamid");
}
DoSQL("UPDATE $opts{bamfiles_table} SET bamname='$newbamfile' WHERE bamid=$bamid");

#   Figure out the MD5 file used and then comment out the proper line
#   Every center has it's own convention, cause that's what standards are for
#   Find all MD5 files and search them, too many ways to fail otherwise
if (opendir(DIR, '.')) {
    my @md5list = ();
    while (readdir(DIR)) {
        if (/\.md5$/) { push @md5list,$_; next; }
        if (/\.MD5$/) { push @md5list,$_; next; }
        if ($_ eq 'Manifest.txt') { push @md5list,$_; next; }
    }
    closedir DIR;
    my $changed = 0;
    my $active = 0;
    foreach my $f (@md5list) {
        # Read this MD5 file and look for checksum line
        open(IN, $f) ||
            die "$Script Unable to open file '$f'\n";
        my $lines = '';
        while(<IN>) {
            if (/$checksum/) { $_ = '# ' . $_; $changed++; }
            if (! /^#/) { $active++; }  # Still have MD5s to do
            $lines .= $_;
        }
        close(IN);
        if ($changed) {
            if ($opts{verbose}) { print "Removing checksum line for $bamfile from '$f'\n"; }
            open(OUT, '>', $f) ||
                die "$Script Unable to overwrite '$f' [$bamid]\n";
            print OUT $lines;
            close(OUT);
            #   If there are no more MD5 sums in this file, rename it
            if (! $active) {
                rename($f, "$f.old") ||
                    print "$Script Unable to rename $f $f.old (bamid=$bamid)\n";
            }
            last;
        }
    }
    if (! $changed) { print "$Script MD5 file not found for $bamfilepath [$bamid] checksum=$checksum\n"; }
}

exit;



#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedrename.pl - Rename a BAM to its nwdid name and update the MD5 file

=head1 SYNOPSIS

  topmedrename.pl  1456 f42dc445bda31bd90e9bace7e2c915aa /incoming/topmed/uw/2015may11.hapmap/89497.bam

=head1 DESCRIPTION

This program renames the original BAM file to it's NWD name.
It will also figure out the MD5 file that contained the checksum
and comment out that line.
If the resulting MD5 file contains no more MD% checksums, then
that file is also renamed.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<bamid>
This is bamid in the database for this BAM.

=item B<checksum>
This is the full checksum value.
This provides us a means to determine where the checksum cam from
and to comment out that line of the file.

=item B<bamfilepath>
This is the fully qualified path to a bamfile.

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

