#!/usr/bin/perl
###################################################################
#
# Name: topmed_phs.pl
#
# Description:
#   Use this program to manage the PHS files and database
#
# ChangeLog:
#   $Log: topmed_phs.pl,v $
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
use My_DB;
use Getopt::Long;
use XML::Simple;

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

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
if (! -d "/usr/cluster/$ENV{PROJECT}") { die "$Script - Environment variable PROJECT '$ENV{PROJECT}' incorrect\n"; }
our %opts = (
    realm => "/usr/cluster/$ENV{PROJECT}/etc/.db_connections/$ENV{PROJECT}",
    bamfiles_table => 'bamfiles',
    phsconfig => '/net/topmed/incoming/study.reference/study.reference/study.phs.numbers.tab',
    phsdir => '/net/topmed/incoming/study.reference/phs',
    phsurl => 'http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetSampleStatus.cgi?study_id=%PHS%&rettype=xml',
    verbose => 0,
    filelist => '',
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose filelist=s dryrun phsdir=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options]";
    warn "$m fetch\n$m update\n$m verify\n" .
        "Manage the PHS files and database\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}

$_ = $opts{filelist};
s/,/ /g;
my @filelist = split(' ', $_);
my $filelist_ref = \@filelist;

my $dbh = DBConnect($opts{realm});
if ($opts{verbose}) { print "$Script realm=$opts{realm}\n"; }

my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

#--------------------------------------------------------------
#   Execute the subcommands requested
#--------------------------------------------------------------
foreach my $fcn (@ARGV) {
    if ($fcn eq 'fetch')  {
        $filelist_ref = Fetch($opts{phsconfig}, $opts{phsdir});
        next;
    }
    if ($fcn eq 'update') {
        Update($filelist_ref);
        next;
    }
    if ($fcn eq 'verify') {
        Verify($filelist_ref);
        next;
    }
    die "$Script - Unknown directive '$fcn'\n";
}

exit;

#==================================================================
# Subroutine:
#   Fetch - Fetch the PHS files 
#
# Arguments:
#   config - Path to the config file
#   dir - Path to directory where files are kept
# Returns:
#   filelist_ref - reference to an array of files
#==================================================================
sub Fetch {
    my ($config, $dir) = @_;
    my $today = strftime('%Y%m%d', localtime);
    my @filelist = ();

    #   Read in config file, fetch PHS file for each PI/study/Center
    #   Unfortunately the tab delimited file has extra blanks scatter through it
    open(IN, $config) ||
        die "Unable to open file '$config': $!\n";
    $_ = <IN>;
    if (! /PI_NAME/) { die "File '$config' does not have the proper header: $_"; }
    my $n = 0;
    while (<IN>) {
        chomp();
        s/\t/ /g;
        my ($piname, $study, $center, $phs) = split(' ', $_);
        if ((! defined($phs)) || $phs eq '.') { next; }
        #   We now have the PHS to fetch
        my $f = $opts{phsdir} . "/$today.$center.$piname.$study.$phs.xml";
        if ($opts{verbose}) { print "Fetch '$phs' to '$f'\n"; }
        my $cmd = $opts{phsurl};
        $cmd =~ s/%PHS%/$phs/;
        $cmd = "wget -o /dev/null -O $f '$cmd'";
        system($cmd) &&
            die "Unable to fetch PHS '$phs'\n";
        $cmd = "gzip -f $f";
        system($cmd) &&
            die "Unable to compress '$f'\n";
        $n++;
        push @filelist, $f . '.gz';
    }
    close(IN);
    print "$nowdate Fetched $n PHS files to '$opts{phsdir}'\n";
    return \@filelist;
} 

#==================================================================
# Subroutine:
#   Update - Fetch the PHS files 
#
# Arguments:
#   filesref - Reference to array of paths to xml files
#==================================================================
sub Update {
    my ($filesref) = @_;

    foreach my $file (@$filesref) {
        my $f = $file;
        if ($f =~ /\.gz$/) {
            $f = "/tmp/$$.tmp";
            my $cmd = "gunzip -c $file > $f";
            system($cmd) &&
                die "Unable to run command: $cmd\n";
        }
        #   Read in each XML             
        my $xml = new XML::Simple();
        my $myxml = eval { $xml->XMLin($f) };
        if ($@) {
            if ($f ne $file) { unlink($f); }
            die "XML parse error in '$file': $@\n";
        }
        #   Walk through array of all NWDID entries
        my $k = $#{$myxml->{Study}->{SampleList}->{Sample}};
        if ($opts{verbose}) { print "  Processing $k NWDID entries\n"; }
        my $changes = 0;
        my $unknowns = 0;
        for (my $i=0; $i<=$k; $i++) {
            my $phs = $myxml->{Study}->{accession};
            if ($phs =~ /^(phs\d+)/) { $phs = $1; }
            my $r = $myxml->{Study}->{SampleList}->{Sample}[$i];    # For convenience
            my $nwdid = $r->{submitted_sample_id};
            my $consent = $r->{consent_short_name} || '';
            my $sra_sample = $r->{sra_sample_id} || '';
            my $sra_details = $r->{sra_data_details} || '';

            #   Update this NWDID in database, else ignore 
            my $sth = DoSQL("SELECT bamid from $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'");
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { 
                if ($opts{verbose}) { print "  No BAMID found for NWDID='$nwdid'\n"; }
                $unknowns++;
                next;
            }
            if ($rowsofdata > 1) { die "Yikes! NWDID '$nwdid' not unique in $opts{bamfiles_table}\n"; }
            my $href = $sth->fetchrow_hashref;

            #   Update with new information
            my $updsql = "UPDATE $opts{bamfiles_table} SET phs='$phs'," .
                "phs_consent_short_name='$consent'," .
                "phs_sra_sample_id='$sra_sample'," .
                "phs_sra_data_details='$sra_details'," .
                "nwdid_known='Y' " .
                "WHERE bamid=$href->{bamid}";
            if ($opts{dryrun}) { print "Did not update with: $updsql\n"; }
            else { DoSQL($updsql); }
            $changes++;
        }
        if ($f ne $file) { unlink($f); }
        print "$nowdate Updated $changes entries from $file, $unknowns unknown bamids\n";  
    }
} 

#==================================================================
# Subroutine:
#   Verify - Compare the state for an NWDID in the PHS files
#       to that in the bamfiles database
#       These might not match because the NCBI summary File
#       is not always correct :-( 
#
# Arguments:
#   filesref - Reference to array of paths to xml files
#==================================================================
sub Verify {
    my ($filesref) = @_;
    my %status = ();                # Save list of loaded files

    foreach my $file (@$filesref) {
        my $f = $file;
        if ($f =~ /\.gz$/) {
            $f = "/tmp/$$.tmp";
            my $cmd = "gunzip -c $file > $f";
            system($cmd) &&
                die "Unable to run command: $cmd\n";
        }
        #   Read in each XML             
        my $xml = new XML::Simple();
        my $myxml = eval { $xml->XMLin($f) };
        if ($@) {
            if ($f ne $file) { unlink($f); }
            die "XML parse error in '$file': $@\n";
        }
        #   Walk through array of all NWDID entries
        my $k = $#{$myxml->{Study}->{SampleList}->{Sample}};
        if ($opts{verbose}) { print "  Processing $k NWDID entries [$file]\n"; }
        my $loaded = 0;
        my $waiting = 0;
        my $unknown = 0;
        my $errors = 0;
        my $processing = 0;
        for (my $i=0; $i<=$k; $i++) {
            my $phs = $myxml->{Study}->{accession};
            if ($phs =~ /^(phs\d+)/) { $phs = $1; }
            my $r = $myxml->{Study}->{SampleList}->{Sample}[$i];    # For convenience
            my $nwdid = $r->{submitted_sample_id};

            if (ref($r->{SRAData}->{Stats}) ne 'HASH') { $unknown++; next; }
            my $status = $r->{SRAData}->{Stats}->{status} || '';
            if ($status eq 'ready') { $status{$nwdid} = $file; $loaded++; }
            else {
                if (! $status) { $unknown++; next; }
                if ($status =~ /waiting/) { $waiting++; next; }
                if ($status eq 'error') { $processing++; next; }
                if ($status eq 'processing') { $errors++; next; }
                print "NWDID=$nwdid status=$status\n";
            }
        }
        print "    Status counts: loaded=$loaded waiting=$waiting unknown=$unknown error=$errors processing=$processing\n"; 
    }

    #   Now check the status of all delivered bams to see if these were loaded
    #   Apparently we cannot tell anything useful about the primary bam
    #   so we can only look for bams that are delivered and are loaded in the PHS file
    my $sth = DoSQL("SELECT bamid,state_ncbib37,expt_sampleid FROM $opts{bamfiles_table} " .
        "WHERE state_ncbiorig=$DELIVERED");
    my $rowsofdata = $sth->rows();
    my $mixed = 0;
    print "Found $rowsofdata entries marked as delivered\n";
    for (my $i=1; $i<=$rowsofdata; $i++) { 
        my $href = $sth->fetchrow_hashref;
        if ($href->{state_ncbib37} == $COMPLETED || $href->{state_ncbib37} == $DELIVERED) { $mixed++; next; }
        my $nwdid = $href->{expt_sampleid};
        if (! exists($status{$nwdid})) { next; }
        #   This NWDID was loaded, but we don't think software
        print "Delivered '$nwdid' marked as LOADED in '$status{$nwdid}'\n";
    }
    print "Found $mixed entries where secondary is complete, but primary was not\n";
} 

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_phs.pl -  manage the PHS files and database

=head1 SYNOPSIS

  topmed_phs.pl  fetch update
  topmed_phs.pl  fetch
  topmed_phs.pl  -filelist '/tmp/Cleveland.xml.gz /tmp/Jackson.xml.gz' update

=head1 DESCRIPTION

This program fetches the PHS files from http://www.ncbi.nlm.nih.gov/projects/...
and extracts information for each NWDID entry.
The files are saved in $opts{phsdir} for update to use.

The files are read and parsed and certain data is saved in our
database to be used later in creating the 
XML files used when data is delivered to NCBI for archiving.

=head1 OPTIONS

=over 4

=item B<-dryrun>

Do everything except actually update the database. Useful for debugging.

=item B<-help>

Generates this output.

=item B<-filelist file>

If specified this is the file of PHS information from fetch.
This option allows one to separate the normal sequence of fetch and update.

=item B<-phsdir directory>

Specifies the directory where the PHS data is downloaded.
This defaults to B<'/net/topmed/incoming/study.reference/phs>.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<fetch>

Fetch the PHS files from NCBI. This is normally followed by B<update>

=item B<update>

Read the PHS files and update the database.

=item B<verify>

Read the PHS files and try to verify the state of files in the database.
This is generally not used, but was added in an attempt to find lost
entries in the NCBI data.

=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

