#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name:	topmed_init.pl
#
# Description:
#   Use this program to check for files that are arriving
#   and initialize the NHLBI TOPMED database
#
# ChangeLog:
#   $Log: nhlbi_init.pl,v $
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
use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    topdir => '/net/topmed/incoming/topmed',
    fakedcomment => 'Created from MD5 data',
    runcount => 0,
    bamcount => 0,
    bamcountruns => '',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose=i center=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] updatedb\n" .
        "Monitor NHLBI data arriving in a directory (default=$opts{topdir}').\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

my $dbh = DBConnect($opts{realm});
my $nowdate = time();
if ($fcn ne 'updatedb') { die "$Script  - Invalid function '$fcn'\n"; }
chdir($opts{topdir}) ||
    die "$Script Unable to CD to '$opts{topdir}': $!\n";

#--------------------------------------------------------------
#   For each center watch for a new run to arrive
#--------------------------------------------------------------
my $centersref = GetCenters();
foreach my $centerid (keys %{$centersref}) {
    my $c = $centersref->{$centerid};
    my $d = $opts{topdir} . '/' . $c;
    if (! chdir($d)) {
        warn "$Script Unable to CD to '$d': $!\n";
        next;
    }
    #   Get all the known batch runs for this center
    my $sql = "SELECT runid,dirname,xmlfound FROM $opts{runs_table} WHERE centerid=$centerid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %knownruns = ();
    my $dir;
    for (my $i=1; $i<=$rowsofdata; $i++) {
        foreach my $href ($sth->fetchrow_hashref) {
            $dir = $href->{dirname};            
            $knownruns{$dir}{runid} = $href->{runid};
            $knownruns{$dir}{xmlfound} = $href->{xmlfound};
        }
    }
    #   Get list of all runs for this center
    #   Find new ones and add details to database
    my $dirsref = GetDirs('.');
    my $runid;
    foreach my $d (@{$dirsref}) {
        $runid = $knownruns{$d}{runid} || CreateRun($centerid, $d);
        if (! defined($runid)) {
            warn "How can runid not be defined?  dir=$d\n";
            next;
        }
        AddBams($runid, $d);
    }
}

$nowdate = strftime('%Y/%m/%d %H:%M', localtime);

if ($opts{runcount}) { print "$nowdate  Added $opts{runcount} runs\n"; }
if ($opts{bamcount}) { print "$nowdate  Added $opts{bamcount} bams from: $opts{bamcountruns}\n"; }
exit;

#==================================================================
# Subroutine:
#   CreateRun - Add details on this run to the database
#   It's not complete, but it can get us going
#
# Arguments:
#   d - directory (e.g. run name)
#   cid - center id
#
# Returns:
#   runid
#==================================================================
sub CreateRun {
    my ($cid, $d) = @_;

    my $sql = "INSERT INTO $opts{runs_table} " .
        "(centerid,dirname,comments,bamcount,dateinit) " .
        "VALUES($cid,'$d','',0,'$nowdate')";
    my $sth = DoSQL($sql);
    my $runid = $sth->{mysql_insertid};
    warn "Added run '$d'\n";
    $opts{runcount}++;
    #   Try to force permissions so things can work later. Can't trust the users
    chmod(0775, $d) || print "Unable to force permissions for '$d'. Too bad for you.\n"; 
    return $runid;
}

#==================================================================
# Subroutine:
#   AddBams - Add details for bams in this directory to the database
#
# Arguments:
#   runid - run id
#   d - directory
#
# Returns:
#   Boolean if any were added or not
#==================================================================
sub AddBams {
    my ($runid, $d) = @_;

    if (! -d $d) {
        warn "Unable to read directory '$d': $!";
        return 0;
    }
    
    #   Get all the known bams for this run
    my $sql = "SELECT bamname FROM $opts{bamfiles_table} WHERE runid=$runid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %knownbams = ();
    for (my $i=1; $i<=$rowsofdata; $i++) {
        foreach my $href ($sth->fetchrow_hashref) {
            $knownbams{$href->{bamname}} = 1;
        }
    }

    #   Merge all the MD5 files together and get list of bams and checksums
    #   There is noconsistency what people do here. Some files are bam checksum,
    #   some are the other way around. Other files have three field and again
    #   the order is anything you can imagine. Grrrr
    my %bams = ();
    my $bamcount = 0;
    if (! open(IN, "cat $d/*.md5 $d/*.MD5 *$d/Manifest.txt 2>/dev/null |")) {
         warn "No MD5 files found in directory '$d'. Maybe later, eh?\n";
         return 0;
    }
    my ($fn, $checksum, $f3);
    while (<IN>) {                  # Sometimes it's file checksum, sometimes not
        ($checksum, $fn, $f3) = split(' ',$_);
        if (! $fn) { warn "Surprising md5 syntax line: $_"; next; }
        if ($f3) { ($fn, $checksum) = ($checksum, $f3); }
        if ($fn !~ /bam$/) { ($checksum, $fn) = ($fn, $checksum); }
        if ($fn =~ /\//) { $fn = basename($fn); }
        if ($fn !~ /bam$/) {
             if ($opts{verbose}) { warn "Ignoring md5 syntax line: $_"; }
             next;
        }
        $bams{$fn} = $checksum;
        $bamcount++;
     }
    close(IN);
    if (! %bams) {
        if ($opts{verbose}) { warn "No bams found in '$d'\n"; }
        return 0;
    }

    #   Generate the bam database entry for each bam
    my $newbams = 0;
    foreach my $f (keys %bams) {
        if (exists($knownbams{$f})) { next; }   # Skip known bams
        $sql = "INSERT INTO $opts{bamfiles_table} " .
            "(runid,bamname,checksum,refname,expt_refname,expt_sampleid,dateinit) " .
            "VALUES($runid,'$f','$bams{$f}','UNKNOWN','UNKNOWN', 'UNKNOWN', $nowdate)";
        $sth = DoSQL($sql);
        $newbams++;
    }

    #   If we added bams, change the bamcount
    if ($newbams) {
        my $n = 
        $sql = "UPDATE $opts{runs_table}  SET bamcount=$bamcount WHERE runid=$runid";
        $sth = DoSQL($sql);
        $opts{bamcount} += $newbams;
        $opts{bamcountruns} .= $d . ' ';
    }
    return 1;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_init.pl - check for files that are arriving and initialize the database

=head1 SYNOPSIS

  topmed_init.pl updatedb

=head1 DESCRIPTION

This program monitors directories for incoming data and then
updates a database with various status values.

=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<topmed>.

=item B<-verbose N>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<updatedb>

This directs the program to monitor the B<-dir> directory for changes
and update the database table specified by B<-realm>.

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

