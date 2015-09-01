#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name: topmed_nwdid.pl
#
# Description:
#   Use this program to create the NWD id file
#
# ChangeLog:
#   $Log: topmed_nwdid.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my($me, $mepath, $mesuffix) = fileparse($0, '\.pl');
(my $version = '$Revision: 1.0 $ ') =~ tr/[0-9].//cd;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my %opts = (
    samtools => '/net/mario/gotcloud/bin/samtools',
    lookuptable => '/net/topmed/incoming/study.reference/study.reference/nhlbi.463.initial.lookup.table.tab',
    topdir => '/net/topmed/incoming/topmed',
    qcresults => '/net/topmed/incoming/qc.results',
    topmedcmd => '/usr/cluster/monitor/bin/topmedcmd.pl',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose bamid=s nonwdid
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$me$mesuffix [options]";
    warn "$m bamfile\n" .
        "\nVersion $version\n" .
        "Create the NWD id file for this bam\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $bamfile = shift @ARGV;
if ($bamfile !~ /^\//) { die "The path to the bamfile '$bamfile' must begin with '/'\n"; }

#--------------------------------------------------------------
#   samtools view -H filename.bam` returns between 100 and 200
#   lines of header information from filename.bam.  We want all
#   lines which begin with "@RG" in column 1, and only the field
#   which begins "SM:" from each line.
#   The value following "SM:" should be the same in every line
#   from one .bam file, otherwise the sequencing center has made
#   a horrible mistake.  So we want one copy of that value,
#   and the number of "@RG" lines in which it occurs, to be
#   associated with the name of the .bam file.
#--------------------------------------------------------------
my $cmd = $opts{samtools} . ' view -H ' . $bamfile;
my $res = `$cmd`;
my $smvalue = '';
my $smcount = 0;
foreach (split("\n", $res)) {
    if (! /^\@RG/) { next; }
    if (! /SM:(\S+)/) { next; }
    my $s = $1;
    if ($smvalue && $s ne $smvalue) { die "SM: fields in BAM '$bamfile' differ\n\n$res\n"; }
    $smvalue = $s;
    $smcount++;
}
if (! $smvalue) { die "BAM '$bamfile' had no SM: fields\n"; }

#--------------------------------------------------------------
#   In the lookup.table.tab file, get the fields "PI_NAME" and
#   "STUDY" (columns 4 and 5) from the line whose first column
#   matches the value of "SM:"
#--------------------------------------------------------------
my ($pi_name, $study);
open(IN,$opts{lookuptable}) ||
    die "Unable to open file '$opts{lookuptable}': $!\n";
while (<IN>) {
    if (! /^$smvalue/) { next; }
    my @cols = split(' ',$_);
    $pi_name = $cols[3];
    $study = $cols[4];
    last;
}
close(IN);

#   For now, any file whose "SM:" value begins with "LP600"
#   is a special case until they get their act together.
if (! $study && $smvalue =~ /^LP600/) {
    $pi_name = 'Barnes';
    $study = 'Asthma_Afr';
}
if (! $study) { die "Unable to find study for '$smvalue' in '$opts{lookuptable}'\n"; }
if (! $pi_name) { die "Unable to find pi_name for '$smvalue' in '$opts{lookuptable}'\n"; }

#--------------------------------------------------------------
#   Update database with this info if bamid provided
#--------------------------------------------------------------
if ($opts{bamid}) {
    $cmd = "$opts{topmedcmd} set $opts{bamid} piname $pi_name";
    system($cmd);
    $cmd = "$opts{topmedcmd} set $opts{bamid} studyname $study";
    system($cmd);
    $cmd = "$opts{topmedcmd} set $opts{bamid} expt_sampleid $smvalue";
    system($cmd);               # Update NWDID
}

#--------------------------------------------------------------
#   Create the NWDID file
#--------------------------------------------------------------
if (! $opts{nonwdid}) {
    my $qcdir = dirname($bamfile);
    $qcdir =~ s:$opts{topdir}:$opts{qcresults}:;
    $cmd = "mkdir -p $qcdir";
    if (! -d $qcdir) {
        system("mkdir -p $qcdir") && die "Unable to make directory: $qcdir\n";
        if ($opts{verbose}) { print "Created directory: $cmd\n"; }
    }

    my $line = "$bamfile $smvalue $smcount $pi_name $study\n";

    my $nwdidfile = $qcdir . '/' . basename($bamfile, '.bam') . '.nwdid';
    #if ($opts{verbose}) {  print "Create file '$nwdidfile' in that directory with this line:\n  $line\n"; }
    open(OUT,'>' . $nwdidfile) ||
        die "Unable to create file '$nwdidfile': $!\n";
    print OUT $line;
    close(OUT);
    if ($opts{verbose}) {  print "Created file '$nwdidfile'\n"; }
}

exit;



#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_nwdid.pl - Create the NWD id file

=head1 SYNOPSIS

  topmed_nwdid.pl  /incoming/topmed/uw/2015may11.hapmap/89497.bam

=head1 DESCRIPTION

This program creates the NWD id file.
The one-line file will be found in 
/incoming/qc.results/CENTER/RUN/BAMFILE_BASENAME.nwdid

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

