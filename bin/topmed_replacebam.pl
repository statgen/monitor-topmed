#!/usr/bin/perl
###################################################################
#
# Name: topmed_replacebam.pl
#
# Description:
#   Every now and then we are sent a bad BAM/CRAM input file.
#   This script is used to replace the bam file and associated
#   data and cause it all to be reprocessed.
#
# ChangeLog:
#   $Log: topmed_replacebam.pl,v $
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
    netdir => '/net/topmed',
    qcresultsdir => 'incoming/qc.results',
    incomingdir => 'incoming/topmed',
    topmedcmd => '/usr/cluster/topmed/bin/topmedcmd.pl',
    topmedmonitor => '/usr/cluster/topmed/bin/topmed_monitor.pl',
    manifestfile => 'Manifest.txt',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help noprompt backupdir=s findnwdid manifest=s verbose 
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 1 || $opts{help}) {
    warn "$Script [options] id newbamfile [newchecksum]\n" .
        "Replace an existing bamfile with another\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my ($bamid, $newbamfile, $newchecksum) = @ARGV;
my ($nwdid);

#--------------------------------------------------------------
#   Normalize the input, make sure we have the correct bamid
#   If needed get the checksum from the manifest
#--------------------------------------------------------------
if ($opts{findnwdid}) {
    if ($newbamfile !~ /(NWD\d+)/) { die "$Script - Unable to find NWDID in '$newbamfile'\n"; }
    $bamid = $1;
}
my $cmd = "$opts{topmedcmd} whatnwdid $bamid";
my $s = `$cmd` || die "$Script - Unknown NWDID/BAMID '$bamid'\n";
my @c = split(' ', $s);
if ($c[0] !~ /^(NWD\d+)\/(\d+)/) { die "$Script - Unknown NWDID/BAMID '$bamid'\n"; }
($nwdid, $bamid) = ($1, $2);
my $run = $c[6];
$run =~ s/\'//g;

if ($newchecksum && length($newchecksum) != 32) {
    die "$Script - checksum '$newchecksum' looks incorrect to me\n";
}
my $checksumsrc = '   [From command line]';
if ((! $newchecksum) && ($opts{manifest})) {    # Look for file in new Manifest
    open(IN, $opts{manifest}) ||
        die "$Script - Unable to open Manifest file '$opts{manifest}': $!\n";
    my $s = '';
    my $bamf = `basename $newbamfile`;
    chomp($bamf);
    while (<IN>) {
        if (! /^(\w{32})\s+(\S+)/) { next; }
        if ($2 ne $bamf) { next; }
        $newchecksum = $1;
        last;
    }
    $checksumsrc = '   [From Manifest file]';
}

if (! $newchecksum) {
    die "$Script - Unable to find a checksum for '$nwdid' [$bamid]\n";
}

if ($opts{backupdir} && (! -d $opts{backupdir})) {
    die "$Script - Backup directory '$opts{backupdir}' does not exist\n";
}

$cmd = "$opts{topmedcmd} show $nwdid bamname_orig";
my $origname = `$cmd` || 'Not Set';
chomp($origname);

$cmd = "$opts{topmedcmd} show $nwdid bamname";
my $bamname = `$cmd` || 'Not Set';
chomp($bamname);

$cmd = "$opts{topmedcmd} show $nwdid bamsize";
my $bamsize = `$cmd` || 'Not Set';
chomp($bamsize);

$cmd = "$opts{topmedcmd} wherepath $nwdid bam";
my $runpath = `$cmd`;
chomp($runpath);
my $manifestfile = "$runpath/$opts{manifestfile}";
my $bamname_filepath = "$runpath/$bamname";
if (! -f $bamname_filepath) { die "$Script - File '$bamname' not found in '$runpath'\n"; }

$cmd = "$opts{topmedcmd} show $nwdid checksum";
my $checksum = `$cmd` || 'Not Set';
chomp($checksum);
@c = stat($newbamfile);

my $newbamsize = $c[7];
if (! $newbamsize) { die "$Script - File '$newbamfile' error: $!\n"; }

my $baiextension = 'bai';
if ($newbamfile =~ /\.cram$/) { $baiextension = 'crai'; }

my $baifile = '';
if ($opts{bai}) {
    $baifile = "$newbamfile.$baiextension";
    if (! -f $baifile) {
        die "$Script - File '$baifile' does not exist\n";
    }
}

#   Ask if we should proceed
my $backupfile = '';
my $manifestbackupfile = '';
if ($opts{backupdir}) {
    $backupfile = "$opts{backupdir}/$bamname.orig";
    $manifestbackupfile = "$opts{backupdir}/$run.$opts{manifestfile}.orig";
}

print "Replace $nwdid [$bamid] \n" .
    "  Run: $run  is at  $runpath\n" .
    "  Bamname: $bamname\n" .
    "  Original Bamname: $origname\n" .
    "  Size: $bamsize\n" .
    "  Checksum: $checksum\n";
if ($opts{backupdir}) {
    print "  Backups will be done to: $opts{backupdir}\n";
    if (-f $backupfile) { print "  Backup BAM/CRAM already exists\n"; }
    else { print "  Backup BAM/CRAM to: $backupfile\n"; }
    if (-f $manifestbackupfile) { print "  Backup $opts{manifestfile} file already exists\n"; }
    else { print "  Backup $opts{manifestfile} to: $manifestbackupfile\n"; }
}
print "with:\n" .
    "  New file src: $newbamfile\n" .
    "  New Size: $newbamsize\n" .
    "  New Checksum: $newchecksum      $checksumsrc\n";
if ($baifile) { print "  New BAI src: $newbamfile.$baiextension\n"; }
if (! $opts{noprompt}) {
    print "\nEnter Y to continue: ";
    $_ = <STDIN>;
    if (! /^y/i) { die "Nothing done\n"; }
}

#--------------------------------------------------------------
#   We are ready to make the changes
#--------------------------------------------------------------
print "Replacement for $bamname starting.  It'd be a great idea to NOT stop this.\n";

if ($opts{backupdir}) {
    if (! -f $backupfile) {
        print "Backing up $bamname\n";
        $cmd = "mv $bamname_filepath $backupfile";
        DoOrDie($cmd);
        print "Backup completed\n";
    }
    else { print "Backup $bamname already done\n"; }
    if (! -f $manifestbackupfile) {
        print "Backing up $opts{manifestfile}\n";
        $cmd = "cp -p $manifestfile $manifestbackupfile";
        DoOrDie($cmd);
        print "Backup $opts{manifestfile} completed\n";
    }
    else { print "Backup $opts{manifestfile} already done\n"; }
}
else { print "No backup was done\n"; }

print "Copying new file to $bamname  (will be slow)\n";
chmod(0600, $bamname_filepath);
$cmd = "cp -p $newbamfile $bamname_filepath";
DoOrDie($cmd);
chmod(0600, $bamname_filepath);             # We can live if this fails
if ($baifile) {
    chmod(0600, "$bamname_filepath.$baiextension");     # Make writable, JIC. Could fail
    $cmd = "cp -p $baifile $bamname_filepath.$baiextension";
    DoOrDie($cmd);
    chmod(0600, "$bamname_filepath.$baiextension");         # We can live if this fails
}
$_ = "Replaced $bamname";
if ($baifile) { $_ .= " and $bamname.$baiextension"; }
print $_ . "\n";

print "Making database changes\n";
$cmd = "$opts{topmedcmd} set $nwdid checksum $newchecksum";
DoOrDie($cmd);

print "Updating $opts{manifestfile}\n";
open(IN, $manifestfile) ||
    die "$Script - Unable to open Manifest file '$manifestfile': $!\n";
my $ss = '';
my $bamf = `basename $newbamfile`;
while (<IN>) {
    if (/$bamf/) { $_ = "$newchecksum\t$bamname\n"; }
    $ss .= $_;
}
close(IN);
chmod(0600, $manifestfile); 
open(OUT, '>' . $manifestfile) ||
    die "$Script - Unable to write new Manifest file '$manifestfile': $!\n";
print OUT $ss;
close(OUT);

foreach my $fcn (qw(arrived md5verified baid qploted cramed)) {
    $cmd = "$opts{topmedcmd} mark $nwdid $fcn requested";
    DoOrDie($cmd);
}
print "Database updated\n";

#   Kick off the initial steps
$cmd = "$opts{topmedmonitor} -run $run arrive";
DoOrDie($cmd);
$cmd = "$opts{topmedmonitor} -run $run verify";
DoOrDie($cmd);

print "Replacement for $bamname completed.\n\n";
exit;      


#==================================================================
# Subroutine:
#   DoOrDie($cmd)
#
#   Execute a command.  If verbose is set, the command is just shown.
#==================================================================
sub DoOrDie {
    if ($opts{verbose}) { print "  Execute: $_[0]\n"; return; }
    system($cmd) && die "$Script - FAILURE  cmd=$cmd\n";
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_replacebam.pl - Replace a bam/cram incoming file

=head1 SYNOPSIS

  topmed_replacebam.pl NWD508042 NWD508042.macrogen.bam 76c3dc5ee5ebce4d1cac300a7dce24bf
  topmed_replacebam.pl 28829 NWD508042.macrogen.bam 76c3dc5ee5ebce4d1cac300a7dce24bf

  topmed_replacebam.pl -manifest Manifest.txt NWD508042 NWD508042.macrogen.bam
  topmed_replacebam.pl -manifest Manifest.txt -findnwdid xx NWD508042.macrogen.bam

  #    CD to where new replacement files are
  #    Create a place for files to be backed up
  #    Walk all BAMs and invoke topmed_replacebam.pl to replace the file
  #    Do not prompt, use the Manifest.txt in replacement files directory
  #    copy the BAI file, find the NWDID from the file to be copied
  cd /net/topmed4/incoming/topmed/macrogen/2016.0729.mathias.29
  for bam in *.bam; do /tmp/topmed_replacebam.pl / 
     -back /net/topmed4/incoming/topmed/macrogen/.2016.0920.files.replaced \ 
     -bai -fin -man Manifest.txt -nop  xx $bam \ 
     | tee -a /net/topmed4/incoming/topmed/macrogen/.2016.0920.files.replaced/log.txt \ 
  done
 
=head1 DESCRIPTION

This program will replace an existing bam or cram in the incoming tree.
The replaced file has probably been processed at some level, so this
program will force re-processing all steps for the new file.

=head1 OPTIONS
 
=over 4

=item B<-backupdir dirpath>

Specifies a directory where the original file will be copied to.
The backup file will use the name of the BAM/CRAM from the database + '.orig'.
If the bam or cram already exists here, it will not be copied a second time.

=item B<-bai>

Specifies there is a BAI or CRAI file that should be copied also.
This file should be named exactly as the input file with the 
proper extension added. E.g. abc.bam and abc.bam.bai.

=item B<-findnwdid>

Parse the new file path and extract the NWDID. Use this as the target NWDID
for the file to be replaced.


=item B<-help>

Generates this output.

=item B<-manifest filepath>

Specifies the path to the Manifest.txt file with the checksum for the new file. If you do not provide this option, then you must provide the checksum value
on the command line.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<id>

This is the NWDID or BAMID for the existing file to be replaced.
If B<-findnwdid> is provided, the NWDID will be parsed from the 
incoming filepath and this will be ignored.
The NWDID will be used to calculate exactly which existing file will be replaced.

=item B<newbamfilepath>

Path to the new bam or cram file. 

=item B<checksum>

This optional parameter specifies the checksum for the new file.
If the option B<-manifest> is provided, this parameter is ignored.

=back


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut


