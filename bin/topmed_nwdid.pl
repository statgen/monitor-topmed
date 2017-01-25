#!/usr/bin/perl
###################################################################
#
# Name: topmed_nwdid.pl
#
# Description:
#   Use this program to extract some data from the BAM header
#   and save it in the database. It always tries to create the NWD id file
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
use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    samtools => '/net/mario/gotcloud/bin/samtools',
    lookuptable => '/net/topmed/incoming/study.reference/study.reference/lookup.table.tab',
    topdir => '/net/topmed/incoming/topmed',
    qcresults => '/net/topmed/incoming/qc.results',
    topmedcmd => '/usr/cluster/topmed/bin/topmedcmd.pl',
    nominal_awk => $Bin . '/topmed_get_nominal.awk',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose bamid=s nonwdid
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options]";
    warn "$m bamfile\n" .
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
if (! $smvalue) { die "BAM '$bamfile' had no SM: field (e.g. no NWDID)\n"; }

#--------------------------------------------------------------
#   The field library_name (e.g. LB:Pond-400614) comes from
#   the header too.
#--------------------------------------------------------------
$cmd = $opts{samtools} . ' view -H ' . $bamfile . ' | grep -e \'^@RG\' | head -1 | tr "\t" "\n" | grep LB: | sed -e s/LB://';
my $library_name = `$cmd`;
chomp($library_name);
if (! $library_name) {
    #   Of course not all centers put this in the same place,
    #   so we have a special hack for illumina :-(
    $cmd = $opts{samtools} . ' view -H ' . $bamfile . ' | grep -e \'^@RG\' | head -1 | tr "\t" "\n" | grep DS: | sed -e s/DS://';
    $library_name = `$cmd`;
    chomp($library_name);
}
if (! $library_name) { die "BAM '$bamfile' had no LB: or DS: field for library_name\n"; }

#--------------------------------------------------------------
#   Get the nominal mean and std.dev of mapped insert size from the header.
#   This returns something like:
#       351.82 83.45 400.00 K 1 203752 151.00 2 201823 151.00 ...
#   where column one is nominal_length, then nominal_sdev
#   After the K are triplets of N and two numbers
#   We want the triplet where the N is 1 and from that triplet 
#   we want the last field (e.g. 151.00)
#--------------------------------------------------------------
$cmd = $opts{samtools} . ' view ' . $bamfile . " 2>/dev/null | awk -f $opts{nominal_awk}";
$res = `$cmd`;
if ($res !~ /^\s*(\S+)\s+(\S+)\s+\S+\s+K\s+(.+)/) { die "BAM '$bamfile' nominal mean and std.dev missing\n  CMD=$cmd/n"; }
my $nominal_length = int((($1*10) + 5)/10);
my $nominal_sdev = int((($2*5) + 2.5)/5);
my $base_coord = '';
my @triplets = split(' ', $3);
for (my $i=0; $i<=$#triplets; $i=$i+3) {
    if ($triplets[$i] ne '1') { next; }
    $base_coord = int($triplets[$i+2]) + 1;
    last;
}
if (! $base_coord) { die "BAM '$bamfile' base_coord missing\n  CMD=$cmd/n"; }

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
    $cmd = "$opts{topmedcmd} set $opts{bamid} library_name $library_name";
    system($cmd);
    $cmd = "$opts{topmedcmd} set $opts{bamid} nominal_length $nominal_length";
    system($cmd);
    $cmd = "$opts{topmedcmd} set $opts{bamid} nominal_sdev $nominal_sdev";
    system($cmd);
    $cmd = "$opts{topmedcmd} set $opts{bamid} base_coord $base_coord";
    system($cmd);
    $cmd = "$opts{topmedcmd} set $opts{bamid} piname $pi_name";
    system($cmd);
    $cmd = "$opts{topmedcmd} set $opts{bamid} studyname $study";
    system($cmd);
    $cmd = "$opts{topmedcmd} set $opts{bamid} expt_sampleid $smvalue";
    system($cmd);               # Update NWDID, weird database name
    $cmd = "$opts{topmedcmd} set $opts{bamid} cramname $smvalue.src.cram";
    system($cmd);
}

#--------------------------------------------------------------
#   Create the NWDID file
#--------------------------------------------------------------
if (! $opts{nonwdid}) {
    my $qcdir = dirname($bamfile);
    if ($qcdir !~ /\/([^\/]+)\/([^\/]+)$/) {    # Grab center and run
        die "$Script - Unable to parse center/run/ from '$qcdir'\n";
    }
    $qcdir = $opts{qcresults} . "/$1/$2";
    $cmd = "mkdir -p $qcdir";
    if (! -d $qcdir) {
        system("mkdir -p $qcdir") && die "Unable to make directory: $qcdir\n";
        if ($opts{verbose}) { print "Created directory: $cmd\n"; }
    }

    my $line = "$bamfile $smvalue $smcount $pi_name $study\n";

    my $b = basename($bamfile);
    if ($b =~ /^(\S+)\./) { $b = $1; }
    my $nwdidfile = $qcdir . "/$b.nwdid";
    if ($opts{verbose}) {  print "Create file '$nwdidfile' in that directory with this line:\n  $line\n"; }
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

topmed_nwdid.pl - Create the NWD id file and up probably update database

=head1 SYNOPSIS

  topmed_nwdid.pl  /incoming/topmed/uw/2015may11.hapmap/89497.bam
  topmed_nwdid.pl  -bamid 7622 /incoming/topmed/uw/2015may11.hapmap/89497.bam

=head1 DESCRIPTION

This program creates the NWD id file.
The one-line file will be found in 
/incoming/qc.results/CENTER/RUN/BAMFILE_BASENAME.nwdid

If B-bamid> is provided, these fields from the BAM header
are saved in the database:
  library_name
  nominal_length
  nominal_sdev
  base_coord
  piname
  studyname
  expt_sampleid    (i.e. NWDID)
  cramname

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

