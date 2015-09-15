#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name: topmed_rename.pl
#
# Description:
#   Use this program to rename the BAM files and associated
#   data in any run.
#   This is a one-time correction to all the NHLBI data
#   we will receive in 2015.
#
# ChangeLog:
#   $Log: topmed_rename.pl,v $
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
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    topdir => '/net/topmed/incoming/topmed',
    qcresults => '/net/topmed/incoming/qc.results',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    dryrun => 1,
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose dryrun
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options]";
    warn "$m center [rundir]\n" .
        "Rename BAM files and their associated data\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $centername = shift @ARGV;
my $rundir = shift(@ARGV) || '';            # Run might not be specified

my $dbh = DBConnect($opts{realm});

my $sql = "SELECT centerid FROM $opts{centers_table} WHERE centername='$centername'";
my $sth = DoSQL($sql);
my $rowsofdata = $sth->rows();
if (! $rowsofdata) { die "$Script - Unknown center '$centername'\n"; }
my $href = $sth->fetchrow_hashref;
my $cid =  $href->{centerid};

my $centerdir = $opts{topdir} . '/' . $centername;
if (! -d $centerdir) { die "$Script - Unknown center directory '$centerdir'\n"; }

my $runsref;
if ($rundir) { 
    my %run2dir = ();
    $sql = "SELECT runid,dirname FROM $opts{runs_table} WHERE dirname='$rundir'";
    $sth = DoSQL($sql);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Unknown run '$rundir'\n"; }
    $href = $sth->fetchrow_hashref;
    $run2dir{$href->{runid}} = $rundir;
    $runsref = \%run2dir;
}
else { $runsref = GetRuns($cid); }
if (! $runsref) { die "Unable to get '$rundir' from database\n"; }

#   Get list of all bams to be renamed
foreach my $runid (keys %{$runsref}) {
    my $renamecount = 0;
    my $dirname = $runsref->{$runid};
    $sql = "SELECT bamid,bamname,expt_sampleid FROM $opts{bamfiles_table} WHERE runid='$runid'";
    $sth = DoSQL($sql);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { next; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        $href = $sth->fetchrow_hashref;
        if (Rename($centername, $dirname, $href->{bamname}, $href->{bamid}, $href->{expt_sampleid})) { $renamecount++; }
    }
    print "\nRenamed $renamecount BAMs in '$centername/$dirname'\n";
}

exit;

#==================================================================
# Subroutine:
#   Rename - Rename a BAM file and it's associated files
#
# Arguments:
#   c - center name
#   d - run directory name
#   b - bam name
#   bamid - bamid
#   nwdid - NWDID
#
# Returns:
#   boolean if renames were done
#==================================================================
sub Rename {
    my ($c, $d, $b, $bamid, $nwdid) = @_;

    my $f = $opts{topdir} . "/$c/$d/$b";
    if (! -f $f) {
        warn "$Script - BAM '$f' was not found\n";
        return 0;
    }
    if ($b =~ /nwd/i) {
        if ($opts{verbose}) { warn "$Script - BAM '$f' appears to be renamed already\n"; }
        return 0;
    }
    #   Create list of commands to execute
    my ($firstpart, $newb);
    if ($b =~ /^([^\.]+)\.(\S+)/) { $firstpart = $1; $newb = "$nwdid.$2"; }
    else { die "$Script - Unable to parse BAM '$b'\n"; }

    #   Create small script to do renames
    $f = $opts{topdir} . "/$c/$d/rename_$nwdid.txt";
    open(OUT, '>' . $f) ||
        die "$Script - Unable to create file '$f'\n";
    print OUT "#  These renames were executed for $b\nset -e\n";
    my @filelist = split("\n", `ls $opts{qcresults}/$c/$d/$firstpart* $opts{topdir}/$c/$d/$firstpart*`);
    if (! @filelist) { die "$Script - No files found for BAM '$b' ??\n"; }
    foreach my $l (@filelist) {
        my $newl = $l;
        $newl =~ s/\/$firstpart\./\/$nwdid\./;
        print OUT "mv $l $newl\n";
    }
    #   We now have a list of commands to run, do so
    my $cmd = "bash $f";
    my $sql = "UPDATE $opts{bamfiles_table} SET bamname='$newb' WHERE bamid='$bamid'";
    print "Executing commands: $f ========\n#   SQL=$sql\n";
    if ($opts{dryrun}) { $cmd = "cat $f"; }
    system($cmd) &&
        die "$Script - Rename failed: CMD=$cmd\n";
    if (! $opts{dryrun}) { $sth = DoSQL($sql); }                #   Finally make changes in the database
    return 1;
}


#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_rename.pl - One time script to rename BAM files to use NWDID

=head1 SYNOPSIS

  topmed_rename.pl  uw
  topmed_rename.pl  uw  2015may11.framingham.trio

=head1 DESCRIPTION

This program renames BAMs and their associated files to use NWDID
for the name of the file, rather than whatever convention was used
by the center itself.

This script will only be to clean up files in 2015

=head1 OPTIONS

=over 4

=item B<-bamid id>

If specified the PI name and study will be set in the database for this bam.

=item B<-help>

Generates this output.

=item B<-nonwdid>

If specified, the nwdid file is not created. This can be useful when redoing other parts of the processing.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<bamfilepath>
This is the fully  qualified path to a bamfile.

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

