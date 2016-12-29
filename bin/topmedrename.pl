#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
###################################################################
#
# Name: topmedrename.pl bamid bamfilepath
#
# Description:
#   Use this program to rename a BAM to it's nwdid name
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
if ($#ARGV < 1 || $opts{help}) {
    warn "$Script [options] bamid bamfilepath\n" .
        "Rename the BAM file to its NWDID name, correct the MD5 file entry.\n";
    exit 1;
}
my $bamid = shift(@ARGV);
my $bamfilepath = shift(@ARGV);
my $dbh = DBConnect($opts{realm});

#--------------------------------------------------------------
#   Rename the BAM file and change bamname in database
#--------------------------------------------------------------
my $sth = DoSQL("SELECT * FROM $opts{bamfiles_table} WHERE bamid='$bamid'");
my $rowsofdata = $sth->rows();
if (! $rowsofdata) { die "Script BAM '$bamid' does not exist in database [$bamfilepath]\n"; }
my $href = $sth->fetchrow_hashref;
my $nwdid = $href->{expt_sampleid};
if ($nwdid !~ /^NWD/) { die "Script BAM '$bamid' NWDID [$nwdid] was not set [$bamfilepath]\n"; }


#   CD to place where BAM exists
my $dirname = dirname($bamfilepath);
chdir($dirname) ||
    die "$Script Unable to CD to '$dirname': $!\n";
my $bamfile = basename($bamfilepath);

#   Rename the BAM file if this is bam and the filename does not start with NWDID
if ($bamfile =~ /\.bam$/ && $bamfile !~ /^NWD/) {
    my $newbamfile = $bamfile;
    if ($newbamfile !~ /^[^.]+\.(.+)/) { die "$Script Unable to parse '$newbamfile'\n"; }
    $newbamfile = $nwdid . '.' . $1;
    #   Aspera screws us up by re-transmitting the file if we rename it
    #   so we create a symlink to the original
    if (! -f $newbamfile) {         # Only do this once
        system("ln -s $bamfile $newbamfile") &&
            die "$Script Unable to create symlink to $bamfile for $newbamfile (bamid=$bamid)\n";
        if ($opts{verbose}) { print "$Script Symlink created to $bamfile for $newbamfile\n"; }
        #   Save original bamname once in database, change name in database
        if ((! defined($href->{bamname_orig})) || $href->{bamname_orig} eq '') {
            DoSQL("UPDATE $opts{bamfiles_table} SET bamname_orig='$bamfile' WHERE bamid=$bamid");
        }
        DoSQL("UPDATE $opts{bamfiles_table} SET bamname='$newbamfile' WHERE bamid=$bamid");
    }
}

#   Rename the BAM file if this is a cram and the extension is just .cram
if ($bamfile =~ /NWD\d+\.cram$/) {
    my $newbamfile = $nwdid . '.src.cram';
    #   Aspera screws us up by re-transmitting the file if we rename it
    #   We ignore that in this case, hoping Aspera dies a long slow death :-)
    system("mv $bamfile $newbamfile") &&
        die "$Script Unable to rename $bamfile to $newbamfile (bamid=$bamid)\n";
    if ($opts{verbose}) { print "$Script Renamed $bamfile as $newbamfile\n"; }
    
    #   If there was a crai provided, rename that too
    my $f = $bamfile . '.crai';
    if (-f $f) {
        system("mv $f $newbamfile.crai") &&
            die "$Script Unable to rename $f to $newbamfile.crai (bamid=$bamid)\n";
        if ($opts{verbose}) { print "$Script Renamed $f as $newbamfile.crai\n"; }
    }
    #   Save original bamname once in database, change name in database
    DoSQL("UPDATE $opts{bamfiles_table} SET bamname_orig='$bamfile' WHERE bamid=$bamid");
    DoSQL("UPDATE $opts{bamfiles_table} SET bamname='$newbamfile' WHERE bamid=$bamid");
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

It turns out Aspera can re-transmit the file if we actually rename it,
causing chaos all around as data re-appears after being processed.
So rather than actually renaming files, we create a symlink, making
two versions of the file, under different names.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the database realm to read data from. This defaults to B<topmed>;

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<bamid>
This is bamid in the database for this BAM.

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

