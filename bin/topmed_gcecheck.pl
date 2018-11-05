#!/usr/bin/env perl
###################################################################
#
# Name: topmed_gcecheck.pl
#
# Description:
#   Use this program to update the database to request that remapped
#   data in Google Cloud be pulled to local store.
#
# ChangeLog:
#   $Log: topmed_gcecheck.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;

use FindBin qw($Script);
use IO::File;
use lib (
  qq($FindBin::Bin),
  qq($FindBin::Bin/../lib),
  qq($FindBin::Bin/../lib/perl5),
  qq($FindBin::Bin/../local/lib/perl5),
  qq(/usr/cluster/topmed/lib/perl5),
  qq(/usr/cluster/topmed/local/lib/perl5),
);

use Topmed::Constants qw(:states);
use Getopt::Long;
use My_DB;

use POSIX qw(strftime);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
#   Pre-check options for project
for (my $i=0; $i<=$#ARGV; $i++) {
    if (($ARGV[$i] eq '-p' || $ARGV[$i] eq '-project') && defined($ARGV[$i+1])) {
        $ENV{PROJECT} = $ARGV[$i+1];
        last;
    }
}
if (! -d "/usr/cluster/$ENV{PROJECT}") { die "$Script - Environment variable PROJECT '$ENV{PROJECT}' incorrect\n"; }

our %opts = (
    verbose => 0,
    gcecachefileprefix => "topmed_gcecheck",
    gsuri => 'gs://topmed-recabs/\*/\*.flagstat',
    realm => "/usr/cluster/$ENV{PROJECT}/etc/.db_connections/$ENV{PROJECT}",
    bamfiles_table => 'bamfiles',
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose cache nocache
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] mark|verify\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);
DBConnect($opts{realm});

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'mark')     { MarkState(@ARGV); exit; }
#if ($fcn eq 'verify')   { VerifyState(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   MarkState()
#
#   Find all cases of remapped samples at Google and mark
#   them to be pulled.  This creates a cache file of flagstats.
#==================================================================
sub MarkState {
    #my ($x) = @_;
    print "$nowdate MarkState\n";

    my $nwdref = CacheGSData();
    my $count = 0;
    my $countcannot = 0;
    foreach my $nwdid (keys %{$nwdref}) {
        #   If we have NWDID, this appears to be ready to be pulled
        my $sql = "SELECT state_gce38push,state_gce38pull FROM $opts{bamfiles_table} " .
            "WHERE expt_sampleid='$nwdid'";
        my $sth = DoSQL($sql);
        my $href = $sth->fetchrow_hashref;
        if (! exists($href->{state_gce38push}) ||
            $href->{state_gce38push} ne $COMPLETED) {
            print "Sample $nwdid exists in GCE, but was not pushed\n";
            $countcannot++;
            next;
        }
        if ($href->{state_gce38pull} == $COMPLETED) {
            print "Sample $nwdid has already been pulled, but still exists in GCE\n";
            $countcannot++;
            next;
        }
        if ($href->{state_gce38pull} == $REQUESTED) { next; }
        if ($href->{state_gce38pull} == $SUBMITTED) { next; }
        if ($href->{state_gce38pull} == $STARTED) { next; }

        #   Mark sample to be pulled
        $sql = "UPDATE $opts{bamfiles_table} SET state_gce38pull=$REQUESTED " .
            "WHERE expt_sampleid='$nwdid'";
        $sth = DoSQL($sql);
        if ($opts{verbose}) { print "Requested PULL for $nwdid\n"; }
        $count++;
    }
    print "Requested $count PULLs to be done, $countcannot NWDIDs could not do PULL\n";
}

#==================================================================
# Subroutine:
#   CacheGSData()
#
#   Get all flagstat files at GCE. Cache them in a local file
#   and return a reference to a hash of NWDID and the GCE file
#==================================================================
sub CacheGSData {
    #my ($x) = @_;
    my %gcedata = ();

    #   Look for an existing cache file that isn't too old
    #   and either use it or remove it and create a new one
    my $now = time();
    my ($out, $in);
    opendir($in, '/tmp') ||
        die "$Script - cannot read directory '/tmp': $!\n";
    while (readdir $in) {
        if (! /$opts{gcecachefileprefix}\.(\d+)/) { next; }
        my $t = $1;
        if (($now - $t) < 60*60*24) { $now = $t; last; }
        else {              # Old cache file, remove it
            my $f = "/tmp/$opts{gcecachefileprefix}.$t"; 
            unlink($f);
            print "Removed old cache file '$f'\n";
        }
    }
    closedir($in);

    #   If a cache exists, read it and return the data.  $now is time for file
    my $count = 0;
    my $flagstatscache = "/tmp/$opts{gcecachefileprefix}." . $now;
    if ((! $opts{nocache}) && -r $flagstatscache) {
        print "Reading GCE data from '$flagstatscache'\n";
        open($in, $flagstatscache) ||
            die "$Script - Unable to read file '$flagstatscache'; $!\n";
        while (<$in>) {
            my @c = split(' ', $_);
            if ($c[0] eq 'NWD000') { next; }    # GCE test sample
            if ($c[0] =~ /^NWD/) {
                $gcedata{$c[0]} = $c[1];
                $count++;
            }
        }
        close($in);
        print "Returned $count GCE samples from cache file\n";
        return \%gcedata;
    }

    #   No cache file, get list of files from GCE    
    open($out, '>' . $flagstatscache) ||
        die "$Script - Unable to create file '$flagstatscache'; $!\n";
    my $cmd = "gsutil ls $opts{gsuri}";    
    print "Very slowly creating cache file '$flagstatscache'\n";
    open($in, "$cmd |") ||
        die "$Script - Unable to get list of flagstat files.  CMD=$cmd\n";
    while (my $l = <$in>) {
        chomp($l);
        if ($l !~ /^gs.*\/(NWD\d+)\/NWD\d*\.recab\.cram\.flagstat/) {
            if ($opts{verbose} && $l =~ /NWD/) { print "$Script - cannot parse $l\n"; }
            next;
        }
        if ($1 eq 'NWD000') { next; }
        $gcedata{$1} = $l;
        print $out "$1 $l\n";           # Save data in cache file
        $count++;
    }
    close($in);
    close($out);
    print "Created cache file '$flagstatscache' of $count samples\n";
    return \%gcedata;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_gcecheck.pl - Manage the state_gce38* flags and remapped files

=head1 SYNOPSIS
 
  topmed_gcecheck.pl mark         # Sets state_gce38pull in database

=head1 DESCRIPTION

This program is used to manage the state_gce38* flags in the database.

=head1 OPTIONS

=over 4

=item B<-cache>

A list of GCE files is saved in B</tmp/topmde_gcecheck.timestamp> the first time.
Subsequent calls will load the GCE file list from this.
Old cache files are automatically removed.

=item B<-help>

Generates this output.

=item B<-nocache>

Force removal of the cache file.

=item B<-project PROJECT>

Specifies these commands are to be used for a specific project.
Warning, this can only be abbreviated as B<-p> or <-project>.
The default is to use the environment variable PROJECT.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program can be:

B<mark>
Scans the files at Google and for every remapped sample, sets the database
flag so the sample will be pulled to local storage.

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut
