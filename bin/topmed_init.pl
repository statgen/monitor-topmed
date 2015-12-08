#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
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
    topmedmd5file => 'topmed_md5.txt',          # Consolidate all MD5s in this
    runcount => 0,
    bamcount => 0,
    bamcountruns => '',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s
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
            warn "$Script How can runid not be defined?  dir=$d\n";
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
    chmod(0775, $d) || print "$Script Unable to force permissions for '$d'. Too bad for you.\n"; 
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
        print "$Script - Unable to read directory '$d': $!";
        return 0;
    }

    #   Get all the known bams for this directory/run   Get NWD name and original name
    my $sql = "SELECT bamname,bamname_orig FROM $opts{bamfiles_table} WHERE runid=$runid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %knownbams = ();
    my %knownorigbam = ();
    for (my $i=1; $i<=$rowsofdata; $i++) {
        foreach my $href ($sth->fetchrow_hashref) {
            $knownbams{$href->{bamname}} = 1;
            if (exists($href->{bamname_orig}) && defined($href->{bamname_orig}) &&
                $href->{bamname_orig} ne '' && $href->{bamname_orig} ne ' ') {
                $knownorigbam{$href->{bamname_orig}} = 1;
            }
        }
    }

    #   Get list of all new MD5 files
    #   There is no consistency what people do here.
    my @md5files = ();
    if ( -r "$d/Manifest.txt") {
        push @md5files,'Manifest.txt';
    }
    else {
        opendir(my $dh, $d) ||
            die "$Script - Unable to read directory '$d'\n";
        while (readdir $dh) {
            if (! /\.md5$/) { next; }
            if ($_ eq $opts{topmedmd5file}) { next; }
            push @md5files,$_;
        }
        closedir $dh;
    }
    if (! @md5files) { return 0; }          # No new md5 files found

    #   Foreach md5 file, get the new of the BAM and create the bamfiles record
    #   Append the md5 record to $opts{topmedmd5file} and rename the md5 file
    #   so we do not process it again
    my $newbams = 0;
    my ($fn, $checksum);
    my $md5lines = '';
    foreach (@md5files) {
        my $f = "$d/$_";
        if (! open(IN, $f)) {
            warn "$Script - Unable to read MD5 file '$f': $!\n";
            next;
        }
        my $badmd5 = 0;
        while (my $l = <IN>) {                      # Sometimes it's file checksum, sometimes not
            if ($l =~ /^#/) { next; }
            ($checksum, $fn) = split(' ',$l);
            #   Do sanity check trying to guess the format of their md5 file
            if ($checksum =~ /\./) { ($fn, $checksum) = ($checksum, $fn); }
            if (length($checksum) < 30) {
                print "$Script - Invalid checksum '$checksum' in '$f'. Line: $l";
                $badmd5++;
                next;
            }
            if ($fn =~ /\//) { $fn = basename($fn); }   # Path should not be qualified, but ...
            if ($fn !~ /bam$/) {                    # File better be a bam
                print "$Script - Invalid BAM name '$fn' in '$f'. Line: $l";
                $badmd5++;
                next;
            }
            $md5lines .= $l;                        # Line OK. Save for our own file

            #   Ideally we only create NEW records, but sometimes we might
            #   reprocess an MD5. If so, don't make duplicate records
            if (exists($knownbams{$fn})) { next; }       # Skip known bams
            if (exists($knownorigbam{$fn})) { next; }    # Skip known original bams

            #   New BAM, create database record for it. May not exist yet
            $sql = "INSERT INTO $opts{bamfiles_table} " .
                "(runid,bamname,checksum,dateinit) " .
                "VALUES($runid,'$fn','$checksum', $nowdate)";
            $sth = DoSQL($sql);
            $newbams++;
        }
        close(IN);
        if ($badmd5) { next; }              # This MD5 was in error

        #   Have added new bam. Rename MD5 file and append $md5lines to our own MD5 file
        open(OUT, '>>' . "$d/$opts{topmedmd5file}") ||
            die "$Script - Unable to append MD5 data to '$d/$opts{topmedmd5file}': $!\n";
        print OUT $md5lines;
        close(OUT);
        rename($f, "$f.old") ||
            die "$Script - Unable to rename '$f' to '$f.old': $!\n";
    }

    #   If we added bams, change the bamcount
    if (! $newbams) { return 0; }
    print "$Script $newbams new bams found in '$d'\n";
    $opts{bamcount} += $newbams;            # Stats for ending msg
    $opts{bamcountruns} .= $d . ' ';

    #   Get number of database records
    $sql = "SELECT bamid FROM $opts{bamfiles_table} WHERE runid=$runid";
    $sth = DoSQL($sql);
    my $numbamrecords = $sth->rows();
    $sql = "UPDATE $opts{runs_table}  SET bamcount=$numbamrecords WHERE runid=$runid";
    $sth = DoSQL($sql);

    #   Last sanity check, see if number of BAM files matches records
    #   This might not always be accurate, but ...
    my $n = `ls $d/*.bam | wc -l`;
    chomp($n);
    if ($n eq $numbamrecords) { print "$Script - Congratulations, # bams = # database records\n"; }
    else { print "$Script - Warning, # bams [$n] != # database records [$numbamrecords].  If data is incoming, this might be OK\n"; }
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

=item B<-verbose>

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

