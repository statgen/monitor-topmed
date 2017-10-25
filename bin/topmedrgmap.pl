#!/usr/bin/perl
###################################################################
#
# Name: topmedrgmap.pl
#
# Description:
#   Use this program to create an rgmap file for cramore so it
#   can rewrite the headers in a cram file.
#   Both this and cramore were written by Hyun
#
# ChangeLog:
#   $Log: topmedrgmap.pl,v $
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

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    samtools => '/usr/cluster/bin/samtools',
    topmedpath => '/usr/cluster/topmed/bin/topmedpath.pl',
    topmedcmd => '/usr/cluster/topmed/bin/topmedcmd.pl',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help checkonly verbose
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options]";
    warn "$m bamid|nwdid [rgmapfile]\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $nwdid = shift @ARGV;
my $rgmapfile = shift @ARGV || '/dev/null';

if ($nwdid !~ /^NWD/) {         # Given bamid, get nwdid
    my $n = `$opts{topmedcmd} show $nwdid expt_sampleid`;
    if (! $n) { die "$Script - unknown NWDID '$nwdid'\n"; }
    chomp($n);
    $nwdid = $n;
}
   
#--------------------------------------------------------------
#   Test if rgmap file needed, create it if needed
#--------------------------------------------------------------
my $cramfile = `$opts{topmedpath} wherefile $nwdid b38`;
if (! $cramfile) { die "$Script - recab.cram was not found for $nwdid\n"; }
chomp($cramfile);

my %h = ();
my $cmd = "$opts{samtools} view -H $cramfile";
my $in;
open($in,"$cmd | grep '^\@RG' |") ||
    die "$Script - Unable to run command: $cmd\n";

my $uniq = 1;
my $atrg = '@RG';
while(<$in>) {
    if (! /^$atrg/ ) { die "$Script - Invalid header line (no $atrg) from $cmd\n"; }
    my ($rg,@F) = split;
    my ($id, $sm) = ('', '');
    foreach my $f (@F) {
	    if ( $f =~ /^ID:/ ) { $id = $f }
	    elsif ( $f =~ /^SM:/ ) { $sm = $f }
    }
    if ($id eq '' || $sm eq '') {
        die "$Script - ERRROR: ID or SM is missing from RG $rg\n$_\n";
    }
    
    my ($smtag,$smid) = split(/:/,$sm);
    if ( $smid ne $nwdid ) {
        die "$Script - ERROR: Sample ID does not match for $nwdid -- $smid observed. This sample must be remapped.\n";
    }
    my $found = 0;
    $id =~ s/^ID://;
    foreach my $key (keys %h) {
	    if ( $id eq $key ) { die "$Script - ERROR: Duplicate RG IDs observed.\n"; }
        #  Jonathan thinks this should be the case:
        if ( $id =~ /^$key\-[A-F0-9]{8}$/ ) { 
	        push(@{$h{$key}}, $id);
	        $found = 1;
	        $uniq = 0;
	        last;
	    }
	    #  Jonathan thinks this should be the case:
	    if ( $key =~ /^$id\-[A-F0-9]{8}$/ ) {
	        die "$Script - ERROR: Strange order of RG observed: $id after $key\n"; }
    }
    if ( $found == 0 ) { $h{$id} = [$id]; }
}
close($in);
if (! %h) { die "$Script - Unable to find any \@RG tag from $nwdid\n"; }

#   In this case the cram header was correct
if ($uniq) {
    print "Header for '$cramfile' [$nwdid] is valid\n";
    exit(0);
}

#   Cram header was wrong, maybe create rgmap file
if ($opts{checkonly}) {
    print "RGMAP for '$cramfile' [$nwdid] is invalid. No RGMAP created.\n";
    exit(1);
}

#   Create rgmap file
my $out;
open($out, '>' . $rgmapfile) ||
    die "$Script - Unable to create file '$rgmapfile' [$nwdid]: $!\n";
foreach my $id (sort keys %h) {
	print $out join("\t",$id,@{$h{$id}}) . "\n";
}
close($out);
print "RGMAP for '$cramfile' [$nwdid] is invalid. Created rgmap file '$rgmapfile'\n";
exit;

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedrgmap.pl - Create an rgmap file to correct headers in a cram

=head1 SYNOPSIS
 
  topmedrgmap.pl 2199 /tmp/2199.rgmap      # Creates rgmap for NWD433184
  topmedrgmap.pl NWD433184 /tmp/2199.rgmap      # Same
  
  topmedrgmap.pl -check 2199 /dev/null     # rc=1 if invalid, rc=0 if valid


=head1 DESCRIPTION

This program is used to create an rgmap file for an invalid cram.

If the input cram (identified by it's NWDID) is valid, no rgmap is created
and the return code is zero.

If the input cram header is NOT valid, an rgmap is created
and the return code is zero.


=head1 OPTIONS

=over 4

=item B<-checkonly>

No rgmap file is created. If the header is valid, the return code is zero, else nonzero.

=item B<-help>

Generates this output.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=item B<nwdid>

Identifies the CRAM header to be checked. This can also be the database bamid.

=item B<rgmap`>

Specifies the path to the rgmap file to be created. This defaults to B</dev/null>.

=back


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by HM Kang and Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut



