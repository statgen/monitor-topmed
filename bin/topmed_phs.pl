#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
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
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/perl5";
use My_DB;
use Getopt::Long;
use XML::Simple;
use POSIX qw(strftime);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    phsconfig => '/net/topmed/incoming/study.reference/study.reference/study.phs.numbers.tab',
    phsdir => '/net/topmed/incoming/study.reference/phs',
    phstable => 'phs',
    phsurl => 'http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetSampleStatus.cgi?study_id=%PHS%&rettype=xml',
    verbose => 0,
    filelist => '',
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose=i filelist=s dryrun
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options]";
    warn "$m fetch|nofetch\n" .
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
    print "Fetched $n PHS files to '$opts{phsdir}'\n";
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
        for (my $i=0; $i<=$k; $i++) {
            my $phs = $myxml->{Study}->{accession};
            if ($phs =~ /^(phs\d+)/) { $phs = $1; }
            my $r = $myxml->{Study}->{SampleList}->{Sample}[$i];    # For convenience
            my $nwdid = $r->{submitted_sample_id};
            my $consent = $r->{consent_short_name} || '';
            my $sra_sample = $r->{sra_sample_id} || '';
            my $sra_details = $r->{sra_data_details} || '';

            #   Update this NWDID in database, else ignore 
            my $sth = DoSQL("SELECT bamid from $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'", 0);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            if ($rowsofdata > 1) { die "Yikes! NWDID '$nwdid' not unique in $opts{bamfiles_table}\n"; }
            my $href = $sth->fetchrow_hashref;

            #   Update with new information
            my $updsql = "UPDATE $opts{bamfiles_table} SET phs='$phs'," .
                "phs_consent_short_name='$consent'," .
                "phs_sra_sample_id='$sra_sample'," .
                "phs_sra_data_details='$sra_details' " .
                "WHERE bamid=$href->{bamid}";
            if ($opts{dryrun}) { print "Did not update with: $updsql\n"; }
            else { DoSQL($updsql); }
            $changes++;
        }
        if ($f ne $file) { unlink($f); }
        print "Updated $changes entries from $file\n";  
    }      
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

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<fetch>

Fetch the PHS files from NCBI. This is normally followed by B<update>

=item B<update>

Read the PHS files and update the database.

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

