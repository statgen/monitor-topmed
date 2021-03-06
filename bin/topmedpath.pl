#!/usr/bin/perl
###################################################################
#
# Name: topmedpath.pl
#
# Description:
#   Use this program to show paths to key directories and files in topmed
#   This was part of topmedcmd.pl but was moved into a new program
#   for stability. This program might return surprising paths
#   to deal with inconsistencies from differing centers.
#   This program can work with topmed and inpsyght
#
# ChangeLog:
#   $Log: topmedpath.pl,v $
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
  qq(/usr/cluster/topmed/lib/perl5),            # So we can run from /tmp
);

use Getopt::Long;
use Topmed::Path;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    verbose => 0,
    datatype => 'genome',
);

Getopt::Long::GetOptions( \%opts,qw(
    help raw realm=s verbose project=s datatype=s
    )) || die "$Script - Failed to parse options\n";
#   Non-typical code for PROJECT so non-project IDs can easily use this
if ($opts{project}) { $ENV{PROJECT} = $opts{project}; }
else {
    if (! exists($ENV{PROJECT})) { $ENV{PROJECT} = 'noproject'; }
}
if (! -d "/usr/cluster/$ENV{PROJECT}") {
   die "$Script - Environment variable PROJECT '$ENV{PROJECT}' incorrect\n";
}

PathInit($ENV{PROJECT}, $opts{datatype});	# Set up project variables, special handling

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [-p project] [-datatype genome|rnaseq|methyl] [options]";
    warn "$m wherepath id KEYWORD\n" .
        "  or\n" .
        "$m whathost id KEYWORD\n" .
        "  or\n" .
        "$m wherefile id KEYWORD\n" .
        "\nWHERE KEYWORD is one of bam|cram|localbackup|remotebackup|remotearchive|" .
        	"qcresults|console|b37|b38|bcf|gceupload|awsupload\n" .
        "  or\n" .
        "WHERE ID is nwdid or sampleid or batchid\n" .
        "WHERE KEYWORD is one of rundir|releasefiles|fileprefix\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'wherepath') { print WherePath(@ARGV) . "\n"; exit; }
if ($fcn eq 'whathost')  { print WhatHost(@ARGV) . "\n"; exit; }
if ($fcn eq 'wherefile') { print WhereFile(@ARGV) . "\n"; exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

1;

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedpath.pl - Show paths for data in the TopMed database

=head1 SYNOPSIS
 
  topmedcmd.pl wherepath 2199 bam          # Returns real path to bam
  topmedcmd.pl wherepath 2199 cram         # Returns path to cram directory
  topmedcmd.pl wherepath 2199 qcresults    # Returns path to directory for qc.results
  topmedcmd.pl wherepath 2199 console      # Returns path to directory for SLURM output
  topmedcmd.pl wherepath 2199 backup       # Returns GCE URI to where backup might be
  topmedcmd.pl wherepath 2199 gceupload    # Returns GCE URI to where remapped CRAMs are copied
  topmedcmd.pl wherepath 2199 awsupload    # Returns AWS URI to where all files are copied

  topmedcmd.pl -proj MYPROJ  whathost 2199 bam   # Returns host for bam
  topmedcmd.pl whathost 2199 cram          # Returns host for cram directory
  topmedcmd.pl whathost 2199 qcresults     # Returns host for directory for qc.results

  topmedcmd.pl wherefile 2199 bam          # Returns path to bam file (may not exist)
  topmedcmd.pl wherefile 2199 cram         # Returns path for cram file (may not exist)
  topmedcmd.pl wherefile 2199 qcresults    # Returns path for qc.results *.vb.SelfSM file (may not exist)
  topmedcmd.pl wherefile 2199 remotelbackup  # Returns GCE URI to where backup file might be
  topmedcmd.pl wherefile 2199 localbackup    # Returns path to local backup (cram)

  topmedcmd.pl -datatype rnaseq wherepath 2199 rundir  # Returns real path to tx_project 

=head1 DESCRIPTION

This program supports simple commands to show the path to key files and directories.
This program provides the path, but does not guarantee the the file/directory exists.

This program requires you to provide the project name (lower case) as either
environment variable or witht he option B<-project>.
Projects must have a directory as /usr/cluster/projectname.

See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-datatype genome|rnaseq|methyl>

Allows one to specify the type of data to work with. 
The default is B<genome>.

=item B<-project projectname>

Allows one to specify the project to work on. You must specify this option
or set the environment variable B<PROJECT}.

=item B<-raw NAME>

Specifies the path returned has not used abs_path(), e.g. has not been verified.
In this case B<-raw> will always return a path, whereas the default behavior
could return a null string if the path is invalid (e.g. NFS mount has failed).

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program are grouped into several groups which are used
to deal with specific sets of information in the monitor databases.
The paths returned may not exist.

B<wherepath bamid|nwdid bam|cram|backup|qcresults|console|b37|b38bcf|gceupload|gcebcfupload|awsupload|awsbucket|awsbucketpath>
If B<bam> was specified, display the path to the real bam file.

If B<cram> or B<localbackup> was specified, display the path to the backup directory.

If B<remotearchive> was specified, display the GCE path (e.g. gs://topmed-archives/...)

If B<remotebackup> was specified, display the GCE path (e.g. gs://topmed-backups/...)

If B<qcresults> was specified, display the path to the directory where
the qc.results for this bamid will be.

If B<console> was specified, display the path to the directory where
the SLURM console output.

If B<b37> was specified, display the path to the directory of remapped data for build 37 results can be found.

If B<b38> was specified, display the path to the directory of remapped data for build 38 results can be found.

If B<bcf> was specified, display the path to the directory of BCF (vt-discover) data.

If B<gceupload> was specified, display the path to the GCE data for remapped CRAMs.

If B<gcebcfupload> was specified, display the path to the GCE data for remapped BCF data.

If B<awsupload> was specified, display the path to the AWS data.
This path may not exist.

If B<awsupload> was specified, display the path to the AWS data.

If B<awsbucket> was specified, display the bucket name for the AWS data.

If B<rundir> was specified for I<datatype=rnaseq>, display the path to the tx_project.

If B<releasefiles> was specified for I<datatype=rnaseq>, display the path to the files for the tx_project.

If B<fileprefix> was specified for I<datatype=rnaseq>, display the path to the files for the tx_project
with the numeric prefix of the files for this sample.


B<whathhost bamid|nwdid key>
returns the host for the key specified.


B<wherefile bamid|nwdid key>
returns the path to the file for the key specified.
This is not valid for I<datatype=rnaseq>.


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

