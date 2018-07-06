package Topmed::Path;
#==================================================================
# Path.pm
#   This is nearly the whole of topmedpath
#   By making this an include, it can be combined with topmedcheck
#   making it run far faster.
#==================================================================
use strict;
use warnings;

use vars qw(@EXPORT_OK @EXPORT @ISA);
use Exporter;

@ISA = qw(Exporter);
@EXPORT_OK = qw(PATHOPTS PathInit WherePath WhatHost WhereFile);
@EXPORT    = qw(PATHOPTS PathInit WherePath WhatHost WhereFile);

use File::Basename;
use My_DB;
use Cwd qw(realpath abs_path);

our %PATHOPTS = (
    bamfiles_table => 'bamfiles',
    centers_table => 'centers',
    runs_table => 'runs',
    netdir => '/net/',
    qcresultsdir => 'incoming/qc.results',
    incomingdir => 'incoming/',
    backupsdir => 'working/backups/incoming/',
    bcfsdir => 'working/candidate_variants',
    consoledir => 'working/',
    wresults37dir => 'working/schelcj/results',
    iresults37dir => 'incoming/schelcj/results',
    results38dir => 'working/mapping/results',
    gcebackupuri => 'gs://',                        # Was topmed-backups
    gcearchiveuri => 'gs://',
    gcebcfuploaduri => 'gs://',                     # For remapped BCF data
    gceuploaduri => 'gs://',                        # For remapped CRAM data
    awsbucket => 'nih-nhlbi-datacommons',
    awsbucketpath => 'UofM/crams/b38',
    awsuploaduri => 's3://nih-nhlbi-datacommons/UofM/crams/b38',
);

#==================================================================
# Subroutine:
#   PathInit($project)
#
#   Updates the %pathopts variables for this project
#==================================================================
sub PathInit {
    my ($project) = @_;
    $PATHOPTS{realm} = "/usr/cluster/$project/etc/.db_connections/$project";
    $PATHOPTS{netdir} .= $project;
    $PATHOPTS{backupsdir} .= $project;
    $PATHOPTS{incomingdir} .= $project;
    $PATHOPTS{consoledir} .= $project . '-output';
    $PATHOPTS{gcebackupuri} .= $project . '-backups';
    $PATHOPTS{gcearchiveuri} .= $project . '-archives';
    $PATHOPTS{gcebcfuploaduri} .= $project . '-bcf';
    $PATHOPTS{gceuploaduri} .= $project . '-irc-share/genomes';
    My_DB::DBConnect($PATHOPTS{realm});        # Set up database
}

#==================================================================
# Subroutine:
#   $path = WherePath($bamid, $set)
#
#   Print paths to various things for bamid based on $set
#==================================================================
my ($bamid, $bamname, $cramname, $piname, $datayear, $nwdid, $rundir, $centername); # For other functions
sub WherePath {
    my ($b, $set) = @_;
    if ((! defined($set) || ! $set)) { $set = 'unset'; }
    $bamid = GetBamid($b);

    #   Get values of interest from the database
    my $sth = My_DB::DoSQL("SELECT * FROM $PATHOPTS{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$main::Script - BAM '$b' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    $bamname = $href->{bamname};
    $cramname = $href->{cramname} || 'CRAMNAME_not_SET';
    $piname = $href->{piname};
    $datayear = $href->{datayear};
    $nwdid = $href->{expt_sampleid};
    my $runid = $href->{runid};
    $sth = My_DB::DoSQL("SELECT centerid,dirname FROM $PATHOPTS{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$main::Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = My_DB::DoSQL("SELECT centername FROM $PATHOPTS{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$main::Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    $centername = $href->{centername};

    #   BAM is in one of those $PATHOPTS{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamfdir = AbsPath("$PATHOPTS{netdir}/$PATHOPTS{incomingdir}/$centername/$rundir");
        return $bamfdir;
    }
 
    #   Find where the cram file lives
    if ($set eq 'cram' || $set eq 'localbackup') {
        my $bakbamfdir = AbsPath("$PATHOPTS{netdir}/$PATHOPTS{backupsdir}/$centername/$rundir");
        return $bakbamfdir;
    }

    #   Find where the GCE backup cram file lives
    if ($set eq 'remotebackup') {
        my $gcebackup = "$PATHOPTS{gcebackupuri}/$centername/$rundir";
        return $gcebackup;
    }

    #   Find where the GCE archive bam file lives
    if ($set eq 'remotearchive') {
        my $gcebackup = "$PATHOPTS{gcearchiveuri}/$centername/$rundir";
        return $gcebackup;
    }

    #   Print where qc.results are 
    if ($set eq 'qcresults') {
        my $qcdir = AbsPath("$PATHOPTS{netdir}/$PATHOPTS{qcresultsdir}/$centername/$rundir");
        return $qcdir;
    }
 
    #   Print where SLURM console output can be found
    if ($set eq 'console') {
        my $dir = "$PATHOPTS{netdir}/$PATHOPTS{consoledir}";
        return $dir;
    }

    #   Try to guess where the b37 remapped CRAM lives
    if ($set eq 'b37') {
        my $file = FindB37($centername, $piname, $nwdid, $bamid);
        return dirname($file);
    }

    if ($set eq 'b38') {
        my $file = FindB38($centername, $piname, $nwdid, $bamid, $datayear);
        return dirname($file);
    }
 
    if ($set eq 'bcf') {
        my $host = 'topmed8';
        my $dir = "/net/$host/$PATHOPTS{bcfsdir}/$piname";
        mkdir $dir,0755 ||
            die "$main::Script - Unable to create bcf path for '$bamid' to '$dir': $!\n";
        return $dir;
    }
 
    if ($set eq 'gcebcfupload') {
        return "$PATHOPTS{gcebcfuploaduri}/$nwdid";
    }

    if ($set eq 'gceupload') {
        return $PATHOPTS{gceuploaduri};
    }

    if ($set eq 'awsupload') {
        return $PATHOPTS{awsuploaduri};
    }

    if ($set eq 'awsbucket') {
        return $PATHOPTS{awsbucket};
    }

    if ($set eq 'awsbucketpath') {
        return $PATHOPTS{awsbucketpath};
    }

    die "$main::Script - Unknown WherePath option '$set'\n";
}

#==================================================================
# Subroutine:
#   $host = WhatHost($bamid, $set)
#
#   Print host name for the path to various things for bamid based on $set
#==================================================================
sub WhatHost {
    my ($bamid, $set) = @_;
    $bamid = GetBamid($bamid);
    my $path = WherePath(@_);
    if (! $path) { return ''; }

    if ($set eq 'bcf') { return 'topmed8'; }

    if ($set eq 'console') { return 'topmed'; }

    #   All these are done the same way
    my %sets = (
        bam => 1,
        cram => 1,
        localbackup => 1,
        qcresults => 1,
        b37 => 1,
        b38 => 1,
        bcf => 1,
        remotebackup => 1,
        remotearchive => 1,
    );
    if (exists($sets{$set})) {
        my $host = $path;
        if ($path =~ /\/net\/([^\/]+)\//) { $host = $1; }
        return $host;
    }

    die "$main::Script - Unknown WhatHost option '$set'\n";
}

#==================================================================
# Subroutine:
#   $filepath = WhereFile($bamid, $set)
#
#   Print paths to various files for bamid based on $set
#==================================================================
sub WhereFile {
    my ($bamid, $set) = @_;
    $bamid = GetBamid($bamid);
    my $path = WherePath(@_);
    if (! $path) { return ''; }

    if ($set eq 'bam') { return $path . "/$bamname"; }
 
    if ($set eq 'cram' || $set eq 'localbackup') { return $path . "/$cramname"; }

    if ($set eq 'remotebackup') { return $path . "/$cramname"; }

    if ($set eq 'remotearchive') { return $path . "/$cramname"; }

    if ($set eq 'qcresults') {
        $bamname =~ s/.bam//;               # Instead of nwdid, maybe it is original bamname
        if (-f "$path/$bamname.vb.selfSM") { return "$path/$bamname.vb.selfSM"; }
        return "$path/$nwdid.vb.selfSM";   # Here is what we really want it to be
    }

    if ($set eq 'b37') { return $path . "/$nwdid.recal.cram"; }

    if ($set eq 'b38') { return $path . "/$nwdid.recab.cram"; }
 
    if ($set eq 'bcf') { return $path . "/$nwdid.bcf"; }

    die "$main::Script - Unknown WhereFile option '$set'\n";
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

    $sth = My_DB::DoSQL("SELECT bamid FROM $PATHOPTS{bamfiles_table} WHERE expt_sampleid='$bamid'");
    $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        $href = $sth->fetchrow_hashref;
        return $href->{bamid};
    }
    die "$main::Script - Invalid bamid or NWDID ($bamid). Try '$main::Script -help'\n";
}

#==================================================================
# Subroutine:
#   FindB37($centername, $piname, $nwdid, $bamid)
#
#   Return full path to B37 cram or null
#==================================================================
sub FindB37 {
    my ($centername, $piname, $nwdid, $bamid) = @_;
    die "$main::Script - there are no local B37 files. See gs://topmed-irc-working/remapping/b37\n";
    my %files = ();
    foreach my $n ('', '2', '3', '4', '5', '6', '7', '9', '10') {     # All possible topmed hosts
        my $p = AbsPath("$PATHOPTS{netdir}$n/$PATHOPTS{wresults37dir}/$centername/$piname/$nwdid/bams/$nwdid.recal.cram");
        if ($p) { $files{$p} = 1; next; }
        $p = AbsPath("$PATHOPTS{netdir}$n/$PATHOPTS{iresults37dir}/$centername/$piname/$nwdid/bams/$nwdid.recal.cram");
        if ($p) { $files{$p} = 1; next; }
    }
    if (! %files) { return ''; }                # Nothing found
    my @filekeys = keys %files;
    if ($#filekeys == 0) { return $filekeys[0]; }   # One file found, return it
    #   Error condition - too many files found
    foreach (@filekeys) { print "   "; system("ls -l $_"); }
    die "$main::Script - Found " . scalar(@filekeys) . " B37 files fir $bamid\n";
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

    if ($ENV{PROJECT} eq 'inpsyght') {
        return "/net/$ENV{PROJECT}/working/mapping/results/$centername/$piname/b38/$nwdid/$nwdid.recab.cram";
    }

    if ($ENV{PROJECT} eq 'topmed' && $datayear ne '3') {
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
    if ($ENV{PROJECT} eq 'topmed' && $datayear ge '3') {
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
    die "$main::Script - FindB38 did not know about datayear=$datayear\n";
}

#==================================================================
# Subroutine:
#   AbsPath($path)
#
#   Returns the absolute path for a path. This is done
#   manually by get the value of a symlink and then
#   replacing any leading ../../ prefix with /net/
#
#   Return full path or null string
#==================================================================
sub AbsPath {
    my ($p) = @_;
    if (-l $p) {
        #   This is symlink, replace prefix with absolute path
        if ($p =~ /^\.\.[\.\/]+([a-z].+)/) {
            $p = '/net/' . $1;
        }
    }
    if (! $PATHOPTS{raw}) {
        return abs_path($p) || '';
    }
    return $p;
}

1;

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

Get.pm

=head1 SYNOPSIS

  use Topmed::Path;
  
=head1 DESCRIPTION

This is nearly the whole of topmedpath.pl.
By making this an include, it can be combined with topmedcheck.pl
making it run far faster.

=head1 Functions

=over 4

=item B<PathInit>

Initializes the project options.

=item B<WherePath()>

Returns a string of the path to a sample.

=item B<WhatHost()>

Returns a string of the host where a sample source exists.

=item B<WhereFile()>

Returns a string of the path to a sample file.

=back

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut
