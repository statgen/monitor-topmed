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
#
# Poor man's regression testing:
#  t=/tmp/topmedpath.pl
#  bams="4890 19099 15495 9702 5870 9670 6948 8615 10466 7740   57141 72522 70732 71329 82982 25804 38013 42819 73094 80438  78764 87128 85385 85873 88522 77264 85869 94934 79050 93186"
#  keys="bam cram localbackup remotebackup remotearchive qcresults console b37 b38 bcf upload"
#  log=/tmp/j
#
#  rm $log; for b in $bams; do  for k in $keys; do echo -n "$b $k  " | tee -a $log; $t whathost $b $k | tee -a $log; done; done
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
  qq(/usr/cluster/topmed/lib/perl5),
  qq(/usr/cluster/topmed/local/lib/perl5),
);

use File::Basename;
use My_DB;
use Getopt::Long;
use Cwd qw(realpath abs_path);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    centers_table => 'centers',
    runs_table => 'runs',
    netdir => '/net/topmed',
    qcresultsdir => 'incoming/qc.results',
    incomingdir => 'incoming/topmed',
    backupsdir => 'working/backups/incoming/topmed',
    bcfsdir => 'working/candidate_variants',
    consoledir => 'working/topmed-output',
    wresults37dir => 'working/schelcj/results',
    iresults37dir => 'incoming/schelcj/results',
    results38dir => 'working/mapping/results',
    gcebackupuri => 'gs://topmed-irc-working/archives',   # Was topmed-backups
    gcearchiveuri => 'gs://topmed-archives',
    #gcebcfuri => 'gs://topmed-bcf',
    gceuploaduri => 'gs://topmed-bcf',
    #gceremapuri => 'gs://topmed-incoming',
    awsbucket => 'nih-nhlbi-datacommons',
    awsbucketpath => 'UofM/crams/b38',
    awsuploaduri => 's3://nih-nhlbi-datacommons/UofM/crams/b38',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options]";
    warn "$m wherepath|where bamid|nwdid KEYWORD\n" .
        "  or\n" .
        "$m whathost bamid|nwdid KEYWORD\n" .
        "  or\n" .
        "$m wherefile bamid|nwdid KEYWORD\n" .
        "\nWHERE KEYWORD is one of bam|cram|localbackup|remotebackup|remotearchive|qcresults|console|b37|b38|bcf|gceupload|awsupload|awsbucket|awsbucketpath\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

DBConnect($opts{realm});

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'where' || $fcn eq 'wherepath') { WherePath(@ARGV); exit; }
if ($fcn eq 'whathost')  { WhatHost(@ARGV); exit; }
if ($fcn eq 'wherefile') { WhereFile(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   $path = WherePath($bamid, $set, $extra1)
#
#   Print paths to various things for bamid based on $set
#==================================================================
sub WherePath {
    my ($bamid, $set, $extra1) = @_;
    if ((! defined($set) || ! $set)) { $set = 'function_missing'; }    # Default

    $bamid = GetBamid($bamid);

    #   Get values of interest from the database
    my $sth = DoSQL("SELECT runid,piname,datayear,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $piname = $href->{piname};
    my $datayear = $href->{datayear};
    my $nwdid = $href->{expt_sampleid};
    my $runid = $href->{runid};
    $sth = DoSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   BAM is in one of those $opts{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamhost = 'none';
        my $bamfdir = abs_path("$opts{netdir}/$opts{incomingdir}/$centername/$rundir") || '';
        print $bamfdir . "\n";
        exit;
    }
 
    #   Find where the cram file lives
    if ($set eq 'cram' || $set eq 'localbackup') {
        my $bakbamfdir = abs_path("$opts{netdir}/$opts{backupsdir}/$centername/$rundir") || '';
        print $bakbamfdir . "\n";
        exit;
    }

    #   Find where the GCE backup cram file lives
    if ($set eq 'remotebackup') {
        my $gcebackup = "$opts{gcebackupuri}/$centername/$rundir";
        print $gcebackup . "\n";
        exit;
    }

    #   Find where the GCE archive bam file lives
    if ($set eq 'remotearchive') {
        my $gcebackup = "$opts{gcearchiveuri}/$centername/$rundir";
        print $gcebackup . "\n";
        exit;
    }

    #   Print where qc.results are 
    if ($set eq 'qcresults') {
        my $qcdir = abs_path("$opts{netdir}/$opts{qcresultsdir}/$centername/$rundir");
        print $qcdir . "\n";
        exit;
    }
 
    #   Print where SLURM console output can be found
    if ($set eq 'console') {
        my $dir = "$opts{netdir}/$opts{consoledir}";
        print $dir . "\n";
        exit;
    }

    #   Try to guess where the b37 remapped CRAM lives
    if ($set eq 'b37') {
        my $file = FindB37($centername, $piname, $nwdid, $bamid);
        print dirname($file) . "\n";
        exit;
    }

    if ($set eq 'b38') {
        my $file = FindB38($centername, $piname, $nwdid, $bamid, $datayear);
        print dirname($file)  . "\n";
        exit;
    }
 
    if ($set eq 'bcf') {
        my $host = 'topmed8';
        my $dir = "/net/$host/$opts{bcfsdir}/$piname";
        mkdir $dir,0755 ||
            die "$Script - Unable to create bcf path for '$bamid' to '$dir': $!\n";
        print $dir . "\n";
        exit;
    }
 
    if ($set eq 'gceupload' || $set eq 'upload') {
        print "$opts{gceuploaduri}/$nwdid\n";
        exit;
    }

    if ($set eq 'awsupload') {
        print "$opts{awsuploaduri}\n";
        exit;
    }

    if ($set eq 'awsbucket') {
        print "$opts{awsbucket}\n";
        exit;
    }

    if ($set eq 'awsbucketpath') {
        print "$opts{awsbucketpath}\n";
        exit;
    }

    die "$Script - Unknown WherePath option '$set'\n";
}

#==================================================================
# Subroutine:
#   $host = WhatHost($bamid, $set)
#
#   Print host name for the path to various things for bamid based on $set
#==================================================================
sub WhatHost {
    my ($bamid, $set, $extra1) = @_;
    if ((! defined($set) || ! $set)) { $set = 'unset'; }    # Default

    $bamid = GetBamid($bamid);

    #   Get values of interest from the database
    my $sth = DoSQL("SELECT runid,piname,datayear,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $nwdid = $href->{expt_sampleid};
    my $runid = $href->{runid};
    my $piname = $href->{piname};
    my $datayear = $href->{datayear};
    $sth = DoSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   BAM is in one of those $opts{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamhost = 'none';
        my $bamfdir = abs_path("$opts{netdir}/$opts{incomingdir}/$centername/$rundir") || '';
        if ($bamfdir =~ /\/net\/([^\/]+)\//) { $bamhost = $1; }
        print $bamhost . "\n";
        exit;
    }
 
    #   Find where the backup CRAM lives and show host
    if ($set eq 'cram' || $set eq 'localbackup') {
        my $backuphost = 'none';
        my $bakbamfdir = abs_path("$opts{netdir}/$opts{backupsdir}/$centername/$rundir");
        if ($bakbamfdir =~ /\/net\/([^\/]+)\//) { $backuphost = $1; }
        print $backuphost . "\n";
        exit;
    }

    #   Print where qc.results are and show host
    if ($set eq 'qcresults') {
        my $qcdir = abs_path("$opts{netdir}/$opts{qcresultsdir}/$centername/$rundir/$nwdid");
        if ($qcdir && $qcdir =~ /\/net\/([^\/]+)\//) { print $1 . "\n"; }
        exit;
    }

    if ($set eq 'b37') {
        my $file = FindB37($centername, $piname, $nwdid, $bamid);
        if ($file =~ /\/net\/([^\/]+)\//) { print $1 . "\n"; }
        exit;
    }

    if ($set eq 'b38') {
        my $file = FindB38($centername, $piname, $nwdid, $bamid, $datayear);
        if ($file =~ /\/net\/([^\/]+)\//){ print $1 . "\n"; }
        exit;
    }

    if ($set eq 'bcf') {
        print 'topmed8' . "\n";
        exit;
    }

    if ($set eq 'console') {
        print 'topmed' . "\n";
        exit;
    }

    die "$Script - Unknown WhatHost option '$set'\n";
}

#==================================================================
# Subroutine:
#   $filepath = WhereFile($bamid, $set)
#
#   Print paths to various files for bamid based on $set
#==================================================================
sub WhereFile {
    my ($bamid, $set) = @_;
    if ((! defined($set) || ! $set)) { $set = 'unset'; }    # Default
    $bamid = GetBamid($bamid);

    #   Get values of interest from the database
    my $sth = DoSQL("SELECT runid,bamname,cramname,piname,datayear,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $bamname = $href->{bamname};
    my $cramname = $href->{cramname} || 'CRAMNAME_not_SET';
    my $piname = $href->{piname};
    my $datayear = $href->{datayear};
    my $nwdid = $href->{expt_sampleid};
    my $runid = $href->{runid};
    $sth = DoSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   BAM is in one of those $opts{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamfile = abs_path("$opts{netdir}/$opts{incomingdir}/$centername/$rundir") || '';
        if ($bamfile) { print $bamfile . "/$bamname" . "\n"; }
        exit;
    }
 
    #   Find where the local backup CRAM lives
    if ($set eq 'cram' || $set eq 'localbackup') {
        my $bakfile = abs_path("$opts{netdir}/$opts{backupsdir}/$centername/$rundir") || '';
        if ($bakfile) { print $bakfile . "/$cramname" . "\n"; }
        exit;
    }

    #   Find where the backup GCE cram file lives
    if ($set eq 'remotebackup') {
        my $gcebackup = "$opts{gcebackupuri}/$centername/$rundir/$cramname";
        print $gcebackup . "\n";
        exit;
    }

    #   Find where the archive GCE bam file lives
    if ($set eq 'remotearchive') {
        my $gcearchive = "$opts{gcearchiveuri}/$centername/$rundir/$bamname";
        print $gcearchive . "\n";
        exit;
    }

    #   Print where qc.results we are interested are
    if ($set eq 'qcresults') {
        my $qcdir = abs_path("$opts{netdir}/$opts{qcresultsdir}/$centername/$rundir");
        $bamname =~ s/.bam//;               # Instead of nwdid, maybe it is original bamname
        if (-f "$qcdir/$bamname.vb.selfSM") {
            print "$qcdir/$bamname.vb.selfSM" . "\n";
            exit;
        }
        print "$qcdir/$nwdid.vb.selfSM" . "\n"; # Here is what we really want it to be
        exit;
    }

    #   Try to guess where the b37 remapped CRAM lives
    if ($set eq 'b37') {
        print FindB37($centername, $piname, $nwdid, $bamid) . "\n";
        exit;
    }

    #   Try to guess where the b38 remapped CRAM lives
    if ($set eq 'b38') {
        print FindB38($centername, $piname, $nwdid, $bamid, $datayear) . "\n";
        exit;
    }
 
    if ($set eq 'bcf') {
        my $host = 'topmed8';
        my $dir = "/net/$host/$opts{bcfsdir}/$piname";
        mkdir $dir,0755 ||
            die "$Script - Unable to create bcf path for '$bamid' to '$dir': $!\n";
        print $dir . "/$nwdid.bcf" . "\n";
        exit;
    }

    die "$Script - Unknown WhereFile option '$set'\n";
}

#==================================================================
# Subroutine:
#   GetBamid($bamid)
#
#   Return bamid for bamid or expt_sampleid
#==================================================================
sub GetBamid {
    my ($bamid) = @_;
    my ($sth, $rowsofdata, $href);
    if ($bamid =~ /^\d+$/) { return $bamid; }

    if ($bamid =~ /^NWD/){
        $sth = DoSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE expt_sampleid='$bamid'");
        $rowsofdata = $sth->rows();
        if ($rowsofdata) {
            $href = $sth->fetchrow_hashref;
            return $href->{bamid};
        }
    }
    if ($bamid !~ /^\d+$/){
        die "$Script - Invalid bamid or NWDID ($bamid). Try '$Script -help'\n";
    }
    return $bamid;
}

#==================================================================
# Subroutine:
#   FindB37($centername, $piname, $nwdid, $bamid)
#
#   Return full path to B37 cram or null
#==================================================================
sub FindB37 {
    my ($centername, $piname, $nwdid, $bamid) = @_;
    die "$Script - there are no local B37 files. See gs://topmed-irc-working/remapping/b37\n";
    my %files = ();
    foreach my $n ('', '2', '3', '4', '5', '6', '7', '9', '10') {     # All possible topmed hosts
        my $p = abs_path("$opts{netdir}$n/$opts{wresults37dir}/$centername/$piname/$nwdid/bams/$nwdid.recal.cram") || '';
        if ($p) { $files{$p} = 1; next; }
        $p = abs_path("$opts{netdir}$n/$opts{iresults37dir}/$centername/$piname/$nwdid/bams/$nwdid.recal.cram") || '';
        if ($p) { $files{$p} = 1; next; }
    }
    if (! %files) { return ''; }                # Nothing found
    my @filekeys = keys %files;
    if ($#filekeys == 0) { return $filekeys[0]; }   # One file found, return it
    #   Error condition - too many files found
    foreach (@filekeys) { print "   "; system("ls -l $_"); }
    die "$Script - Found " . scalar(@filekeys) . " B37 files fir $bamid\n";
}

#==================================================================
# Subroutine:
#   FindB38($centername, $piname, $nwdid, $bamid, $datayear)
#
#   The path returned here is more convoluted than we want
#   because we started allocating b38 files on 9 and 10
#   and then after allocating a large number of samples
#   realized we needed more space.
#   Hence we must support two schemes to determine the path.
#
#   Then came a decision to remap year 3 samples, requiring
#   even more space. Ugh!
#
#   Return full path to B38 cram or null
#   Caller must make the directory tree
#==================================================================
sub FindB38 {
    my ($centername, $piname, $nwdid, $bamid, $datayear) = @_;
    if (! defined($datayear)) { $datayear=''; }

    if ("$datayear" ne '3') {
        my @host_partialpath = (
            [ qw/topmed10 topmed9 topmed6  topmed7  topmed9  topmed10/ ],
            [ qw/working  working incoming incoming incoming incoming/ ]
        );
        # First determine the old path and if that directory exists, use it
        my $mod = $bamid % 2;
        my $file = "/net/$host_partialpath[0][$mod]/$host_partialpath[1][$mod]/" .
            "mapping/results/$centername/$piname/b38/$nwdid/$nwdid.recab.cram";
        if ( -f $file) { return $file; }        # If file exists, return that

        # File does not exist, allocate using the new scheme
        $mod = $bamid % 6;
        $file = "/net/$host_partialpath[0][$mod]/$host_partialpath[1][$mod]/" .
            "mapping/results/$centername/$piname/b38/$nwdid/$nwdid.recab.cram";
        return $file;
    }
    #   Yet another brand new scheme invented for year 3
    if ("$datayear" eq '3') {
        my @host_partialpath = (
            [ qw/topmed  topmed2 topmed3 topmed7 topmed5  topmed6  topmed7 topmed9/ ],
            [ qw/working working working working incoming incoming incoming incoming/ ]
        );
        # First determine the old path and if that directory exists, use it
        my $mod = $bamid % 8;
        my $file = "/net/$host_partialpath[0][$mod]/$host_partialpath[1][$mod]/" .
            "mapping/results/$centername/$piname/b38/$nwdid/$nwdid.recab.cram";
        return $file;
    }
    die "$Script - FindB38 did not know about datayear=$datayear\n";
}

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
  topmedcmd.pl wherepath 2199 gceupload    # Returns GCE URI to where all files are copied
  topmedcmd.pl wherepath 2199 awsupload    # Returns AWS URI to where all files are copied

  topmedcmd.pl whathost 2199 bam           # Returns host for bam
  topmedcmd.pl whathost 2199 cram          # Returns host for cram directory
  topmedcmd.pl whathost 2199 qcresults     # Returns host for directory for qc.results

  topmedcmd.pl wherefile 2199 bam          # Returns path to bam file (may not exist)
  topmedcmd.pl wherefile 2199 cram         # Returns path for cram file (may not exist)
  topmedcmd.pl wherefile 2199 qcresults    # Returns path for qc.results *.vb.SelfSM file (may not exist)
  topmedcmd.pl wherefile 2199 remotelbackup  # Returns GCE URI to where backup file might be
  topmedcmd.pl wherefile 2199 localbackup    # Returns path to local backup (cram)

=head1 DESCRIPTION

This program supports simple commands to show the path to key files and directories.
This program provides the path, but does not guarantee the the file/directory exists.

See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

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

B<wherepath bamid|nwdid bam|cram|backup|qcresults|console|b37|b38bcf|gceupload|awsupload|awsbucket|awsbucketpath>
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

If B<gceupload> was specified, display the path to the GCE data.

If B<awsupload> was specified, display the path to the AWS data.
This path may not exist.

If B<awsupload> was specified, display the path to the AWS data.

If B<awsbucket> was specified, display the bucket name for the AWS data.

If B<awsbucketpath> was specified, display the path of data in the AWS bucket.

B<whathhost bamid|nwdid key>
returns the host for the key specified.

B<wherefile bamid|nwdid key>
returns the path to the file for the key specified.


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

