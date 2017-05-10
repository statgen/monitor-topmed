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
use feature qw(say);

use TopMed_Get;
use My_DB;
use Getopt::Long;
use Cwd qw(realpath abs_path);

use Topmed::DB;

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
    gcebackupuri => 'gs://topmed-backups',
    gcebcfuri => 'gs://topmed-bcf',
    gceuploaduri => 'gs://topmed-bcf',
    gceremapuri => 'gs://topmed-incoming',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options] [-persist]";
    warn "$m wherepath|where bamid|nwdid KEYWORD\n" .
        "  or\n" .
        "$m whathost bamid|nwdid KEYWORD\n" .
        "  or\n" .
        "$m wherefile bamid|nwdid KEYWORD\n" .
        "  or\n" .
        "$m whatbamid bamname\n" .
        "WHERE KEYWORD is one of bam|cram|localbackup|remotebackup|qcresults|console|b37|b38|bcf|upload\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

DBConnect($opts{realm});
my $schema = Topmed::DB->new();

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
    my $sample = $schema->resultset('Bamfile')->find($bamid);

    #   Get values of interest from the database
    my $sth = ExecSQL("SELECT runid,piname,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $piname = $href->{piname};
    my $nwdid = $href->{expt_sampleid};
    my $runid = $href->{runid};
    $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   BAM is in one of those $opts{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamhost = 'none';
        my $bamfdir = abs_path("$opts{netdir}/$opts{incomingdir}/$centername/$rundir") || '';
        print "$bamfdir\n";
        exit;
    }
 
    #   Find where the cram file lives
    if ($set eq 'cram' || $set eq 'localbackup') {
        my $bakbamfdir = abs_path("$opts{netdir}/$opts{backupsdir}/$centername/$rundir") || '';
        print "$bakbamfdir\n";
        exit;
    }

    #   Find where the GCE backup cram file lives
    if ($set eq 'remotebackup') {
        my $gcebackup = "$opts{gcebackupuri}/$centername/$rundir";
        print "$gcebackup\n";
        exit;
    }

    #   Print where qc.results are 
    if ($set eq 'qcresults') {
        my $qcdir = abs_path("$opts{netdir}/$opts{qcresultsdir}/$centername/$rundir") || '';
        print "$qcdir\n";                       # File might not really exist
        exit;
    }
 
    #   Print where SLURM console output can be found
    if ($set eq 'console') {
        my $dir = "$opts{netdir}/$opts{consoledir}";
        print "$dir\n";
        exit;
    }

    if ($set eq 'b37') {
        die "$Script - Unable to provide path for b37 right now\n";
        exit;
    }

    if ($set eq 'b38') {
        my $meth = qq{${set}_mapped_path};
        my $f = $sample->$meth . ".$nwdid.recab.cram";
        say $sample->$meth;
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
 
    if ($set eq 'upload') {
        my $gceupload = "$opts{gceuploaduri}/$nwdid";
        print "$gceupload\n";
        exit;
    }

    die "$Script - Unknown Where option '$set'\n";
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
    my $sample = $schema->resultset('Bamfile')->find($bamid);

    #   Get values of interest from the database
    my $sth = ExecSQL("SELECT runid,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $nwdid = $href->{expt_sampleid};
    my $runid = $href->{runid};
    $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   BAM is in one of those $opts{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamhost = 'none';
        my $bamfdir = abs_path("$opts{netdir}/$opts{incomingdir}/$centername/$rundir") || '';
        if ($bamfdir =~ /\/net\/([^\/]+)\//) { $bamhost = $1; }
        print "$bamhost\n";
        exit;
    }
 
    #   Find where the backup CRAM lives and show host
    if ($set eq 'cram' || $set eq 'localbackup') {
        my $backuphost = 'none';
        my $bakbamfdir = abs_path("$opts{netdir}/$opts{backupsdir}/$centername/$rundir");
        if ($bakbamfdir =~ /\/net\/([^\/]+)\//) { $backuphost = $1; }
        print "$backuphost\n";
        exit;
    }

    #   Print where qc.results are and show host
    if ($set eq 'qcresults') {
        my $qchost = 'none';
        my $qcdir = abs_path("$opts{netdir}/$opts{qcresultsdir}/$centername/$rundir") || '';
        if ($qcdir =~ /\/net\/([^\/]+)\//) { $qcdir = $1; }
        print "$qcdir\n";
        exit;
    }

    if ($set eq 'b37') {
      say $sample->host;
      exit;
    }

    if ($set eq 'b38') {
        my $meth = qq{${set}_mapped_path};
        my $f = $sample->$meth . ".$nwdid.recab.cram";
        my $h = 'none';
        if ($f =~ /\/net\/([^\/]+)\//) { $h = $1; }
        say $h;
        exit;
    }

    if ($set eq 'bcf') {
        my $host = 'topmed8';
        print $host . "\n";
        exit;
    }

    die "$Script - Unknown Where option '$set'\n";
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
    my $sample = $schema->resultset('Bamfile')->find($bamid);

    #   Get values of interest from the database
    my $sth = ExecSQL("SELECT runid,bamname,cramname,piname,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $bamname = $href->{bamname};
    my $cramname = $href->{cramname} || 'CRAMNAME_not_SET';
    my $piname = $href->{piname};
    my $nwdid = $href->{expt_sampleid};
    my $runid = $href->{runid};
    $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   BAM is in one of those $opts{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamfile = abs_path("$opts{netdir}/$opts{incomingdir}/$centername/$rundir") || '';
        if ($bamfile) {
            print $bamfile . "/$bamname\n";
        }
        else {
            warn "WhereFile: abs_path($opts{netdir}/$opts{incomingdir}/$centername/$rundir) failed\n";
        }
        exit;
    }
 
    #   Find where the local backup CRAM lives
    if ($set eq 'cram' || $set eq 'localbackup') {
        my $bakfile = abs_path("$opts{netdir}/$opts{backupsdir}/$centername/$rundir") || '';
        if ($bakfile) {
            print $bakfile . "/$cramname\n";
        }
        else {
            warn "WhereFile: abs_path($opts{netdir}/$opts{backupsdir}/$centername/$rundir) failed\n";
        }
        exit;
    }

    #   Find where the backup GCE cram file lives
    if ($set eq 'remotebackup') {
        my $gcebackup = "$opts{gcebackupuri}/$centername/$rundir/$nwdid.src.cram";
        print "$gcebackup\n";
        exit;
    }

    #   Print where qc.results we are interested are
    if ($set eq 'qcresults') {
        my $qcdir = abs_path("$opts{netdir}/$opts{qcresultsdir}/$centername/$rundir") || '';
        if ($qcdir) { 
            $bamname =~ s/\.bam//;                   # Remove the extension
            $bamname =~ s/\.cram//;
            print "$qcdir/$bamname.vb.selfSM\n";
        }
        else {
            warn "WhereFile: abs_path($opts{netdir}/$opts{qcresultsdir}/$centername/$rundir) failed\n";
        }
        exit;
    }

    #   Try to guess where the b37 remapped CRAM lives
    if ($set eq 'b37') {
        my %files = ();
        foreach ('', '2', '3', '4', '5', '6', '7', '8', '9', '10') {
            my $p = abs_path("$opts{netdir}$_/$opts{wresults37dir}/$centername/$piname/$nwdid/bams/$nwdid.recal.cram") || '';
            if ($p) { $files{$p} = 1; next; }
            $p = abs_path("$opts{netdir}$_/$opts{iresults37dir}/$centername/$piname/$nwdid/bams/$nwdid.recal.cram") || '';
            if ($p) { $files{$p} = 1; next; }
        }
        my @found = keys %files;
        if (! @found) { exit; }           # Nothing found
        if ($#found == 0) { print $found[0]; exit; }
        print "$Script - Found " . scalar(@found) . " files: \n";
        foreach (@found) { print "   "; system("ls -l $_"); }
        exit(1);
    }

    #   Try to guess where the b38 remapped CRAM lives
    if ($set eq 'b38') {
        my $meth = qq{${set}_mapped_path};
        my $f = $sample->$meth . "/$nwdid.recab.cram";
        say $f;
        exit;
    }
 
    if ($set eq 'bcf') {
        my $host = 'topmed8';
        my $dir = "/net/$host/$opts{bcfsdir}/$piname";
        mkdir $dir,0755 ||
            die "$Script - Unable to create bcf path for '$bamid' to '$dir': $!\n";
        print $dir . "/$nwdid.bcf\n";
        exit;
    }

    die "$Script - Unknown Where option '$set'\n";
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
        $sth = ExecSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE expt_sampleid='$bamid'");
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
#   ExecSQL($sql, $die)
#
#   Execute SQL.  Keep trying if connection lost.
#
#   Returns handle for SQL
#==================================================================
sub ExecSQL {
    my ($sql, $die) = @_;
    if ($opts{persist}) { return PersistDoSQL($opts{realm}, $sql); }
    return DoSQL($sql, $die);
}

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
  topmedcmd.pl wherepath 2199 upload       # Returns GCE URI to where all files are copied

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

B<wherepath bamid|nwdid bam|cram|backup|qcresults|console|b37|b38>
If B<bam> was specified, display the path to the real bam file.

If B<cram> or B<localbackup> was specified, display the path to the backup directory.

If B<remotebackup> was specified, display the GCE path (e.g. gs://topmed-backups/...)

If B<qcresults> was specified, display the path to the directory where
the qc.results for this bamid will be.

If B<console> was specified, display the path to the directory where
the SLURM console output.

If B<b37> was specified, display the path to the directory of remapped data for build 37 results can be found.

If B<b38> was specified, display the path to the directory of remapped data for build 37 results can be found.

B<whathhost bamid|nwdid bam|cram|qcresults>
returns the host for the bam, backup, cram or qc.results for the bam file.

B<wherefile bamid|nwdid bam|cram|backup|qcresults>
returns the path to the file for the bam, backup. cram or qc.results. This file may not exist


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

