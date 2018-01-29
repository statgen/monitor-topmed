#!/usr/bin/perl
###################################################################
#
# Name: topmed_gce_manifest.pl
#
# Description:
#   Use this program to create a manifest for GCE and AWS
#
# ChangeLog:
#   $Log: topmed_gce_manifest.pl,v $
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

use Getopt::Long;
use My_DB;

my $COMPLETED = 20;                 # Task completed successfully

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options] ";
    warn "$m aws | gce\n" .
        "Create manifest for AWS or GCE SHARE repositories\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

DBConnect($opts{realm});

#--------------------------------------------------------------
#   Collect list of samples which have been copied to AWS or GCE
#--------------------------------------------------------------
my $outfile = "/tmp/$fcn.manifest.data-commons-pilot.txt";
if ($fcn eq 'gce') {
    SelectFromDB($fcn, $outfile, 'state_gce38copy', 'gs://topmed-irc-share/genomes');
    exit;
}
if ($fcn eq 'aws') {
    die "$Script - AWS is not ready yet\n";
    SelectFromDB($fcn, $outfile, 'state_gce38copy', 'gs://topmed-irc-share/genomes');
    exit;
}

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   SelectFromDB($cloud, $file, $col, $cramuri)
#
#   Get list of samples we have delivered to the cloud
#==================================================================
sub SelectFromDB {
    my ($cloud, $file, $col, $cramuri) = @_;

    #   Make sure this is a bam we know
    my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE $col=$COMPLETED";
    my $sth = DoSQL($sql);
    if (! $sth) { die "$Script - SQL failure. SQL=$sql\n"; }
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - No samples found for '$cloud'. How can that be?\n"; }

    #   Create manifest
    my ($out, $href, $c, @d);
    open($out, '>' . $file) ||
        die "$Script - Unable to create file '$file': $!\n";
    foreach (1 .. $rowsofdata) {
        $href = $sth->fetchrow_hashref;
        $c = "$cramuri/$href->{expt_sampleid}.b$href->{build}.irc.v1.cram",
        @d = (
            $href->{expt_sampleid},
            $href->{studyname} || 'MISSING',
            $href->{piname} || 'MISSING',
            $href->{phs} || 'notregistered',
            $href->{datayear} || 'MISSING',
            $c,
            "$c.crai"
        );
        print $out join("\t",@d) . "\n";
    }
    close($out);
    print "Created manifest for $rowsofdata '$cloud' samples in file '$file'\n";
}


#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_gce_manifest.pl - Create manifest for files released to the cloud

=head1 SYNOPSIS

  topmed_gce_manifest.pl aws
  topmed_gce_manifest.pl gce

 
=head1 DESCRIPTION

This program creates the manifest file for samples in GCE SHARE or AWS.
It does not copy the file to the correct URI.

See B<perldoc DBIx::Connector> for details defining the database.


=head1 OPTIONS

=over 4

=item B<-help>

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<aws>

Creates the manifest for files that have been copied to AWS.

=item B<gce>

Creates the manifest for files that have been copied to GCE SHARE.

=back


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2018 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

