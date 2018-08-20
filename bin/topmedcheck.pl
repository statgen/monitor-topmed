#!/usr/bin/perl
###################################################################
#
# Name: topmedcheck.pl
#
# Description:
#   Do a serious check of data for a set of samples
#
# ChangeLog:
#   $Log: topmedcheck.pl,v $
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
use Topmed::Constants qw(:states);
use Topmed::Get;
use Topmed::Path;

#my $NOTSET = 0;                     # Not set
#my $REQUESTED = 1;                  # Task requested
#my $SUBMITTED = 2;                  # Task submitted to be run
#my $STARTED   = 3;                  # Task started
#my $DELIVERED = 19;                 # Data delivered, but not confirmed
#my $COMPLETED = 20;                 # Task completed successfully
#my $CANCELLED = 89;                 # Task cancelled
#my $FAILED    = 99;                 # Task failed

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    bamfiles_table => 'bamfiles',
    centers_table => 'centers',
    runs_table => 'runs',
    gceuri => 'gs://topmed-bcf',
    gsutil => '/usr/bin/gsutil',
    b37gceuri => 'gs://topmed-irc-working/remapping/b37',
    b38gcebackup => 'gs://topmed-backups',      # Current place for backups
    prevgcebackup => 'gs://topmed-archives',    # Previous place for backups (why 2? this is dumb)
    gcels => '/usr/bin/gsutil ls',   # Add $opts{gceuri}/NWDxxxxxxx
    redotool => '/usr/cluster/topmed/bin/topmed_redo.sh',
    samtools => '/usr/cluster/bin/samtools',
    cacheage => 60*60*24,            # If cache file older than this, refetch
    cachefile => "/tmp/$Script.gcefiles", 
    remapstatus => 'curl --silent --insecure https://104.198.71.226/api/sample-status\\?ids=',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help register fixer=s execfixer max=i redo verbose nfsmounts replace cachefile=s
    all center=s run=s piname=s studyname=s datayear=i sample=s build=i 
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($opts{help} || ($#ARGV<0) && (! $opts{nfsmounts})) {
    my $m = "$Script [options] -subset value  db|localfiles|awsfiles|gcefiles|b37|recabfiles\n" .
        "Where \n" .
        "  [-subset] is -center, -run or -sample\n" .
        "\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $action = shift @ARGV;
if ($opts{redo}) { $opts{fixer} = $opts{redotool}; }
if ($0 =~ /\/(\w+)check/) { PathInit($1); } # Set up project variables;


#--------------------------------------------------------------
#   Get list of bamids to check
#--------------------------------------------------------------
my @bamidrange = ();
my @nwdids = ();
my $morewhere = '';
if ($opts{piname})    { $morewhere .= " AND piname='$opts{piname}'"; }
if ($opts{studyname}) { $morewhere .= " AND studyname='$opts{studyname}'"; }
if ($opts{datayear})  { $morewhere .= " AND datayear=$opts{datayear}"; }
if ($opts{build})     { $morewhere .= " AND build=$opts{build}"; }
if ($opts{b37files})  { $morewhere .= " AND state_b37=20"; }
#if ($opts{recabfiles}){ $morewhere .= " AND state_gce38push=20 AND state_gce38pull!=20"; }
if ($opts{max})       { $morewhere .= " LIMIT $opts{max}"; }

#==================================================================
#   Resume normal processing
#==================================================================
if ($opts{nfsmounts}) { Check_NFSMounts(); exit; }  # Special case processing

if ($opts{sample}) {                    # Maybe one or a range of bamids specified
    if ($opts{samplel} eq 'all') {
        my $sql = "SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} " .
            "WHERE state_arrive=20" . $morewhere;
        my $sth = My_DB::DoSQL($sql);
        GetArrayOfSamples($sth);
    }
    else {
        my $bamid = $opts{sample};
        my @input = ();
        if ($bamid =~/^(\d+)\-(\d+)/) { @input = $1 .. $2; }
        else { push @input, $bamid; }
        for my $b (@input) {
            my $sql = "SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} WHERE ";
            if ($b =~ /^NWD/) { $sql .= "expt_sampleid='$b'"; }
            else { $sql .= "bamid=$b"; }
            $sql .= $morewhere;
            my $sth = My_DB::DoSQL($sql);
            GetArrayOfSamples($sth);
        }
    }
}
if ($opts{center}) {                        # Range specified by center or run
    my $centersref = GetCenters();          # Only returns one center
    foreach my $cid (keys %{$centersref}) {
        $opts{centername} = $centersref->{$cid};    # Save center name
        my $runsref = GetRuns($cid) || next;
        #   Get bamid for every run
        foreach my $runid (keys %{$runsref}) {
            my $sql = "SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} " .
                "WHERE runid=$runid AND state_cram=20" . $morewhere;
            my $sth = My_DB::DoSQL($sql);
            GetArrayOfSamples($sth);
        }
    }
}
#   Get all bamids for this run
if ($opts{run}) {
    if ($opts{run}=~/\D+/) {               # Get runid from dirname
        my $sql = "SELECT runid FROM $opts{runs_table} WHERE dirname='$opts{run}'";
        my $sth = My_DB::DoSQL($sql);
        if (! $sth) { die "$Script - Unknown run '$opts{run}'\n"; }
        my $href = $sth->fetchrow_hashref;
        $opts{run} = $href->{runid};
    }
    #   Get bamid for this run
    my $sql = "SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} " .
        "WHERE state_cram=20 AND runid=$opts{run}" . $morewhere;
    my $sth = My_DB::DoSQL($sql);
    GetArrayOfSamples($sth);
}

if (! @bamidrange) { die "$Script - No samples were found\n"; }
if ($opts{max}) { splice(@bamidrange,$opts{max}); } # We might have selected too many samples          
print "Checking " . scalar(@bamidrange) . " samples\n";

#--------------------------------------------------------------
#   Do varying checks over a set of bamids
#--------------------------------------------------------------
my %emptyok = ();                   # Cache directories with missing data
my $didsomething = 0;
if ($action eq 'db') {
    $didsomething++;
    print "Verify database values\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_DB($bamidrange[$i], $nwdids[$i]);
    }
    print "Completed check of db\n";
}
if ($action eq 'recabfiles') {
    $didsomething++;
    print "Verify remapped files GCE recabs\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_RecabFiles($bamidrange[$i], $nwdids[$i]);
    }
    print "Completed check of recabfiles\n";
}
if ($action eq 'localfiles') {
    $didsomething++;
    print "Verify local files\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_LocalFiles($bamidrange[$i], $nwdids[$i]);
    }
    if (%emptyok) {
        my @ds = sort keys %emptyok;
        print "Warning: " . scalar(@ds) . " directories were unexectedly empty:\n";
        if ($opts{verbose}) { print "  " . join("\n  ", @ds) . "\n"; }
    }
    print "Completed check of localfiles\n";
}
if ($action eq 'gcefiles') {
    $didsomething++;
    print "Verify GCE files\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_GCEFiles($bamidrange[$i], $nwdids[$i]);
    }
    print "Completed check of gcefiles\n";
}
if ($action eq 'awsfiles') {
    $didsomething++;
    print "Verify Amazon files\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_AWSFiles($bamidrange[$i], $nwdids[$i]);
    }
    print "Completed check of awsfiles\n";
}
if ($action eq 'b37') {
    $didsomething++;
    print "Verify B37 GCE files\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_B37Files($bamidrange[$i], $nwdids[$i]);
    }
    print "Completed check of b37files\n";
}
if ($action eq 'backups') {
    $didsomething++;
    print "Verify backups are in GCE\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_Backups($bamidrange[$i], $nwdids[$i]);
    }
    print "Completed check of backups\n";
}
#if ($action eq 'checkfiles') {
#    $didsomething++;
#    print "Verify the local files (md5 and flagstat) - patience!\n";
#    for (my $i=0; $i<=$#bamidrange; $i++) {
#        Check_Files($bamidrange[$i], $nwdids[$i]);
#    }
#    print "Completed checkfiles\n";
#}

#if ($opts{checkbcf}) {
#    print "Verify the local BCF files and index - patience!\n";
#    for (my $i=0; $i<=$#bamidrange; $i++) {
#        Check_BCFFiles($bamidrange[$i], $nwdids[$i]);
#    }
#}

#$ bcftools index -f /net/topmed8/working/candidate_variants/Ellinor/NWD694535.bcf
#[E::bgzf_uncompress] inflate failed: progress temporarily not possible, or in() / out() returned an error
#[E::bgzf_read_block] inflate_block error -1
#[E::bgzf_read] bgzf_read_block error -1 after 9 of 32 bytes
#index: failed to create index for "/net/topmed8/working/candidate_variants/Ellinor/NWD694535.bcf"
#$ echo $?
#255

if (! $didsomething) { die "$Script - Unable to figure out what to check. Try $Script -help\n"; }
exit;

#==================================================================
# Subroutine:
#   Check_DB()
#   Verify database columns (e.g. not bookkeeping)
#==================================================================
sub Check_DB {
    my ($bamid, $nwdid) = @_;
    my @mustexist = (
        'bamname',
        'cramname',
        'base_coord',
        'library_name',
        'cramchecksum',
        'studyname',
        'piname',
        'bamflagstat',
        'cramflagstat',
        'checksum',
        'expt_sampleid',
        'nwdid_known',
        'poorquality',
        'bamsize',
        'datayear',
        'build',
        'offsite'
    );
    my @forcenull = (
        'bamflagstat',
        'cramflagstat',
        'b27flagstat',
        'b38flagstat',
        'cramchecksum',
        'b37cramchecksum',
        'b38cramchecksum',
    );
    my @error = ();

    #   Get all data for this sample
    my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE bamid=$bamid";
    my $sth = My_DB::DoSQL($sql);
    if (! $sth) { die "$Script - Unable to get data for $sql\n"; }
    my $href = $sth->fetchrow_hashref;
    my $year = $href->{datayear};

    #   Certain columns may not be blank
    foreach my $c (@mustexist) {
        if (! defined($href->{$c})) { push @error,"Column $c may not be NULL"; }
    }
    #   Force undefined columns to blank to avoid Perl issue
    foreach my $c (@forcenull) {
        if (! defined($href->{$c})) { $href->{$c} = ''; }
    }
    
    #   Certain columns must have specific values
    my $c = 'expt_sampleid';
    if (! defined($href->{$c})) { $href->{$c} = ''; }
    if ($href->{$c} !~ /^NWD\d{6}$/) {
        push @error,"Column $c $href->{$c} must be NWDnnnnnn";
    }

    $c = 'bamflagstat';
    if ($href->{$c} !~ /^\d+$/) { push @error,"Column $c $href->{$c} must be numeric"; }

    if ($href->{state_cram} == $COMPLETED) {
        $c = 'cramflagstat';
        if ($href->{$c} !~ /^\d+$/) { push @error,"Column $c $href->{$c} must be numeric"; }
        if ($href->{bamflagstat} ne $href->{$c}) {
            push @error,"Column $c $href->{$c} must match bamflagstat $href->{bamflagstat}";
        }
    }

    if ($href->{state_b37} == $COMPLETED) {
        $c = 'b37flagstat';
        if (! defined($href->{$c})) { $href->{$c} = ''; }
        if ($href->{$c} ne $href->{bamflagstat}) {
            if ($href->{$c} eq '') { $href->{$c} = 'x'; }
            push @error,"Column $c $href->{$c} must match bamflagstat $href->{bamflagstat}";
        }
    }

    if ($href->{state_b38} == $COMPLETED) {
        $c = 'b38flagstat';
        if (! defined($href->{$c})) { $href->{$c} = ''; }
        if ($href->{$c} ne $href->{bamflagstat}) {
            if ($href->{$c} eq '') { $href->{$c} = 'x'; }
            push @error,"Column $c $href->{$c} must match bamflagstat $href->{bamflagstat}";
        }
    }

    foreach my $c (qw /nwdid_known poorquality offsite/) {
        if ($href->{$c} !~ /^Y|N|D$/) {
            push @error,"Column $c $href->{$c} must be N or Y or D";
        }
    }

    $c = 'bamsize';
    if (! defined($href->{$c})) { $href->{$c} = ''; }
    my $minbamsize = 10000000000;
    if ($href->{$c} !~ /^\d+/ && $href->{$c} < $minbamsize) {
        push @error,"Column $c $href->{$c} must be > $minbamsize";
    }

    $c = 'datayear';
    if (! defined$href->{$c}) { $href->{$c} = ''; }
    if ($href->{$c} !~ /^1|2|3|4$/) {
        push @error,"Column $c $href->{$c} must be 1, 2 or 3";
    }

    $c = 'build';
    if (! defined($href->{$c})) { $href->{$c} = ''; }
    if ($href->{$c} !~ /^37|38$/) {
        push @error,"Column $c $href->{$c} must be 37 or 38";
    }

    #   All four "phs" cols depend on NCBI data, not exist until study is registered
    if (! $href->{phs} && $opts{register}) { print "BAMID=$bamid NWD=$nwdid is not registered\n"; }

    #   If there were problems, show messages
    if (@error) { ShowErrors($bamid, $nwdid, \@error); }
    else {
        if ($opts{verbose}) { print "  BAMID=$bamid NWD=$nwdid - Database clean\n"; }
    }
}

#==================================================================
# Subroutine:
#   Check_LocalFiles()
#   Verify local files exist
#==================================================================
sub Check_LocalFiles {
    my ($bamid, $nwdid) = @_;
    my @error = ();

    #   Get some data for this sample
    my $sql = "SELECT datayear,state_b37,state_b38,state_qplot FROM $opts{bamfiles_table} WHERE bamid=$bamid";
    my $sth = My_DB::DoSQL($sql);
    if (! $sth) { die "$Script - Unable to get data for $sql\n"; }
    my $href = $sth->fetchrow_hashref;
    my $year = $href->{datayear};

    #   Check original file (could be bam or cram)
    my $f = WhereFile($bamid, 'bam');
    my $e = VerifyFile($f);
    if ($e) { push @error,"BAM cram/index is invalid: $e: $f"; }
    else {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid BAM OK: $f\n"; }
    }
    
    #   Check cram, either original file or converted from bam
    $f = WhereFile($bamid, 'cram');
    $e = VerifyFile($f);
    if ($e) { push @error,"CRAM index is invalid: $e: $f"; }
    else {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid CRAM OK: $f\n"; }
    }

    #   QPlot files should exist
    if ($href->{state_qplot} == $COMPLETED) {
        $f = WhereFile($bamid, 'qcresults');
        if (! $f) { die "Unable to determine qcresults file for $bamid\n"; }
        if ( -f $f ) {
            if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid QCRESULTS OK: $f\n"; }
        }
        else {
            #   We normally expect qcresults to exist, but sometimes they get deleted
            my $d = dirname($f);
            if (! exists($emptyok{$d})) {
                if (-f "$d/.emptyok") { $emptyok{$d} = 1; }
                else { push @error,"QCRESULTS missing: $f"; }
            }
        }
    }

    #   B37 crams should exist for ALL year 1 samples
    #   Some other year samples were also mapped
    #   Note, b37 files have been moved to GCE gs://topmed-irc-working/remapping/b37
    if (0 && $year eq '1' && $href->{state_b37} != $COMPLETED) {
        push @error,"Warning: Year $year sample was not mapped for build B37";
    }
    if (0 && $href->{state_b37} == $COMPLETED) {
        $f = WhereFile($bamid, 'b37');
        my $e = VerifyFile($f);
        if ($e) { push @error,"B37 cram/index is invalid: $e: $f"; }
        else {
            if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid B38 OK: $f\n"; }
        }
    }

    #   B38 crams must exist for all samples
    if ($href->{state_b38} == $COMPLETED) {
        $f = WhereFile($bamid, 'b38');
        my $e = VerifyFile($f);
        if ($e) { push @error,"B38 cram/index is invalid: $e: $f"; }
        else {
            if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid B38 OK: $f\n"; }
        }
    }

    #   BCF must exist for all samples
    if (! defined($href->{state_gce38bcf})) { $href->{state_gce38bcf} = 0; }
    if ($href->{state_gce38bcf} == $COMPLETED) {
        $f = WhereFile($bamid, 'bcf');
        my $e = VerifyFile($f);
        if ($e) { push @error,"BCF/index is invalid: $e: $f"; }
        else {
            if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid BCF OK: $f\n"; }
        }
    }

    #   If there were problems, show messages
    if (@error) { ShowErrors($bamid, $nwdid, \@error); }
    else {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid  - Local files all exist\n"; }
    }
}

#==================================================================
# Subroutine:
#   Get_Cache_Data()
#   Data of interest looks like
#     1454900  2018-02-04T09:30:18Z  gs://topmed-bcf/NWD103302/NWD103302.bcf.csi
#   Get GCE listing of files in local cache file
#
#   Return:
#       Reference to hash of array of files for each sample NWD
#==================================================================
our %cachelines = ();                # Saved cache data here
our $current_cache_file = '';
sub Get_Cache_Data {
    my ($file, $gsuri) = @_;    

    #   Nothing changed and we already did this
    if (! defined($current_cache_file)) { $current_cache_file = ''; }
    if ($file eq $current_cache_file && %cachelines) { return \%cachelines; }
    %cachelines = ();

    #   If file is old, redo
    if ($opts{replace}) { unlink($file); }  # Force rebuild of cache file
    if (-f $file) {
        my @s = stat($file);
        if ($s[9] > (time() + $opts{cacheage})) {
            print "#### Shall I replace cache file '$file' (y/n): ";
            my $a = <STDIN>;
            if ($a !~ /^y/) { die "$Script - cache file not removed. Stopiing\n"; }
            unlink($file);
        }
    }
    if (! -f $file) {
        my $k = 0;
        my $IN;
        my $cmd = $opts{gcels} . ' -r ' . $gsuri . '/';
        my $t = time();
        print "Fetching cache file $file using '$cmd'\n";
        open($IN, "$cmd |") ||
            die "$Script - Unable to fetch GCE data: CMD=$cmd\n";
        while (my $l = <$IN>) {
            chomp($l);
            if ((!defined($l)) || (! $l)) { next; }
            if ($l =~ /INFO/) { next; }     # gsutil had to retry
            if ($l =~ /[:\/]$/) { next; }
            if ($l =~ /\/([^\/]+)\/([^\/]+)$/) {
                my ($d, $f) = ($1, $2);
                if ($d !~ /^NWD/) {             # Studyname, get NWDID from file
                    if ($f =~ /(\w+)\./) { $d = $1; }
                }
                push @{$cachelines{$d}},$f;
                $k++;
                #if ($k > 1000) { last; }       # Handy for debugging
                next;
            }
            warn "Unrecognized line: $l\n";
        }
        close($IN);
        my $OUT;
        open($OUT, '>' . $file) ||
            die "$Script - Unable to create file '$file': $!\n";
        foreach my $ky (keys %cachelines) {
            print $OUT $ky . ' ' . join(' ', @{$cachelines{$ky}}) . "\n";
        }
        close($OUT);
        $t = time() - $t;
        print "Created cache file '$file' with $k files in $t secs\n";
    }
    else {
        my $IN;
        open($IN, $file) ||
            die "$Script - Unable to read GCE cache from '$file': $!\n";    
        my $k = 0;
        while (my $l = <$IN>) {
            $k++;
            my @c = split(' ', $l);
            my $k = shift(@c);
            $cachelines{$k} = \@c;
        }
        close($IN);
        print "Read cache file '$file' with $k samples\n";

    }
    $current_cache_file = $file;        # Remember what we just did
    return \%cachelines;
}


#==================================================================
# Subroutine:
#   Check_GCEFiles()
#   Verify files in Google Cloud Environment
#==================================================================
sub Check_GCEFiles {
    my ($bamid, $nwdid) = @_;
    my @error = ();

    #   Do not bother if we have not copied data to GCE
    my $sql = "SELECT state_gce38copy,state_gce38cpbcf FROM $opts{bamfiles_table} WHERE bamid=$bamid";
    my $sth = My_DB::DoSQL($sql);
    my $href = $sth->fetchrow_hashref;
    if ($href->{state_gce38copy} != $COMPLETED) {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid not copied to GCE yet\n"; }
        return;
    }
    my %requiredextensions = ( 'recab.cram' => 1, 'recab.cram.crai' => 1);
    if ($href->{state_gce38cpbcf} != $COMPLETED) {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid BCF data not copied to GCE yet, continuing\n"; }
    }
    else {
        $requiredextensions{bcf} = 1;
        $requiredextensions{'bcf.csi'} = 1;
    }

    #   All data is supposed to be in GCE
    my $cacheref = Get_Cache_Data($opts{cachefile}, $opts{gceuri});
    if (! exists($cacheref->{$nwdid})) {
        push @error,"No $nwdid GCE files found";
        @{$cacheref->{$nwdid}} = ();        # Create empty array so things work
    }    
    foreach my $f (@{$cacheref->{$nwdid}}) {
        chomp($f);
        if ($f !~ /^$nwdid\.(\S+)/) {
        push @error,"GCE file is invalid: Unknown file: $f"; }
        else {
            my $x = $1;                 # Extension of file
            if (exists($requiredextensions{$x})) { $requiredextensions{$x} = 0; }
            else { push @error,"GCE: Invalid extension: $f"; }
        }
    }
    foreach my $x (keys %requiredextensions) {
        if ($requiredextensions{$x}) { push @error,"GCE file missing: $nwdid.$x"; }
    }

    #   If there were problems, show messages
    if (@error) { ShowErrors($bamid, $nwdid, \@error); }
    else {
        if ($opts{verbose}) { print "  BAMID=$bamid NWD=$nwdid - GCE files all exist\n"; }
    }
}


#==================================================================
# Subroutine:
#   Check_B37Files()
#   Verify all the B37 files exist
#==================================================================
sub Check_B37Files {
    my ($bamid, $nwdid) = @_;
    my @error = ();

    my %requiredextensions = ( 'recal.cram' => 1, 'recal.cram.crai' => 1, 'recal.cram.md5' => 1, 'recal.cram.flagstat' => 1);

    #   All data is supposed to be in GCE
    my $cacheref = Get_Cache_Data($opts{cachefile}, $opts{b37gceuri});
    if (! exists($cacheref->{$nwdid})) {
        push @error,"No $nwdid GCE files found";
        @{$cacheref->{$nwdid}} = ();        # Create empty array so things work
    }    
    foreach my $f (@{$cacheref->{$nwdid}}) {
        chomp($f);
        if ($f !~ /^$nwdid\.(\S+)/) {
        push @error,"GCE file is invalid: Unknown file: $f"; }
        else {
            my $x = $1;                 # Extension of file
            if (exists($requiredextensions{$x})) { $requiredextensions{$x} = 0; }
            else { push @error,"GCE: Invalid extension: $f"; }
        }
    }
    foreach my $x (keys %requiredextensions) {
        if ($requiredextensions{$x}) { push @error,"GCE file missing: $nwdid.$x"; }
    }

    #   If there were problems, show messages
    if (@error) { ShowErrors($bamid, $nwdid, \@error); }
    else {
        if ($opts{verbose}) { print "  BAMID=$bamid NWD=$nwdid - GCE files all exist\n"; }
    }
}

#==================================================================
# Subroutine:
#   Check_Backups()
#   Make sure backup files are in GCE as we expect
#==================================================================
sub Check_Backups {
    my ($bamid, $nwdid) = @_;
    my @error = ();

    #   Original files can be backed up to either of two places
    my $sql = "SELECT state_backup,bamname,build,datayear,runid FROM $opts{bamfiles_table} WHERE bamid=$bamid";
    my $sth = My_DB::DoSQL($sql);
    my $href = $sth->fetchrow_hashref;
    if ($href->{state_backup} != $COMPLETED) {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid was not backed up yet\n"; }
        return;
    }
    my $bamname = $href->{bamname};
    my $build = $href->{build};
    my $datayear = $href->{datayear};
    my $runid = $href->{runid};
    my $ext = '???';
    if ($bamname =~ /\.(\w+)$/) { $ext = $1; }
    my $gsuri = $opts{prevgcebackup};                # Assume previous backup
    my $cachefile = $opts{cachefile} . '.prev';
    # Same logic as found in topmed_gcebackup.sh
    if ($ext eq 'cram' && $datayear >= 3 && $build == 38) {
        $gsuri = $opts{b38gcebackup};               # Use recent backup
        $cachefile = $opts{cachefile};
    }

    #   All data is supposed to be in GCE
    my $backedup = 0;
    my $cacheref = Get_Cache_Data($cachefile, $gsuri);
        if (exists($cacheref->{$nwdid})) {
        foreach my $f (@{$cacheref->{$nwdid}}) {    # All files for NWD
            if ($f =~ /\.$ext$/) { $backedup++; last; }
        }
    }
    if (! $backedup) { push @error,"BACKUP original file not found [$bamname $build $datayear $runid $gsuri]"; }

    #   If there were problems, show messages
    if (@error) { ShowErrors($bamid, $nwdid, \@error); }
    else {
        if ($opts{verbose}) { print "  BAMID=$bamid NWD=$nwdid - GCE backup file exists\n"; }
    }
}

#==================================================================
# Subroutine:
#   Check_RecabFiles()
#   Verify files to be remapped
#==================================================================
sub Check_RecabFiles {
    my ($bamid, $nwdid) = @_;
    my @error = ();
    die "$Script - this is not ready yet\n";

    my $tmpfile = '/run/shm/topmedcheck.tmp';

    my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE bamid=$bamid";
    my $sth = My_DB::DoSQL($sql);
    my $href = $sth->fetchrow_hashref;

    my ($center, $dirname) = GetCenterRun($href->{runid});

    my %summary = ( bamid => $bamid, nwdid => $nwdid );

    if ($href->{state_gce38push} = $COMPLETED) { $summary{gcepush} = 'pushed'; }
    else { $summary{gcepush} = 'not_pushed'; }
    
    if ($href->{state_gce38pull} = $COMPLETED) {$summary{gcepull} = 'pulled'; }
    else { $summary{gcepull} = 'not_pulled'; }

    my $s = $opts{remapstatus} . $nwdid . " | grep NWD";
    $s = `$s`;
    if (! $s) { $summary{remapstatus} = $nwdid . "_unknown"; }
    else {
        if ($s =~ /post-aligned/) { $summary{remapstatus} = 'remapped_done'; }
        else {
            if ($s =~ /failed-pre-align/) { $summary{remapstatus} = "failed-pre-align"; }
            else { $summary{remapstatus} = "unknown_status"; }
        }
    }

    my $f = WhereFile($bamid, 'b38');
    if (-f $f) { $summary{local_file} = 'local_exists'; }
    else { $summary{local_file} = 'local_missing'; }

    my $rc = system("$opts{gcels} gs://topmed-incoming/$center/$dirname/$nwdid.src.cram > $tmpfile 2>/dev/null");
    if ($rc == 0) { $summary{to_be_remapped} = 'incoming_exists';}
    else { $summary{to_be_remapped} = 'incoming_missing'; }

    $rc = system("$opts{gcels} gs://topmed-recabs/$nwdid > $tmpfile 2>/dev/null");
    if ($rc == 0) { $summary{to_be_pulled} = 'recab_exists'; }
    else { $summary{to_be_pulled} = 'recab_missing'; }

    $rc = system("$opts{gcels} gs://topmed-bcf/$nwdid > $tmpfile 2>/dev/null");
    if ($rc == 0) { $summary{gcecopy} = 'gcecopy_exists'; }
    else { $summary{gcecopy} = 'gcecopy_missing'; }

    my @cols = sort keys %summary;
    print join(',',@cols) . "\n";
    $s = '';
    foreach my $k (@cols) { $s .= $summary{$k} . ','; }
    chop($s);
    print "$s\n";
    return;

    print "    gcepush state:  $href->{state_gce38push}\n";
    print "    gcepull state:  $href->{state_gce38pull}\n";
    print "     remap status: "; system($opts{remapstatus} . $nwdid);
    print "       local file: "; system("ls -l $f");
    print "   to be remapped: ";
      system("$opts{gcels} gs://topmed-incoming/$center/$dirname/$nwdid.src.cram");
    print "     to be pulled: "; system("$opts{gcels} gs://topmed-recabs/$nwdid");
    print "          gcecopy: "; system("$opts{gcels} gs://topmed-bcf/$nwdid\\*");
    print "\n";
}

#==================================================================
# Subroutine:
#   Check_Files()
#   Calculate the md5sum and flagstat value for files in this sample
#==================================================================
sub Check_Files {
    my ($bamid, $nwdid) = @_;
    my @error = ();
    die "$Script - this is not ready yet\n";
}

#==================================================================
# Subroutine:
#   Check_AWSFiles()
#   Verify files in Amazon Environment
#==================================================================
sub Check_AWSFiles {
    my ($bamid, $nwdid) = @_;
    my @error = ();
    die "$Script - this is not ready yet\n";

}

#==================================================================
# Subroutine:
#   GetArrayOfSamples($sth)      Follows $sth = My_DB::DoSQL($sql)
#
#   Sets values in @bamidrange and @nwdids
#==================================================================
sub GetArrayOfSamples {
    my ($sth) = @_;

    #my $rowsofdata = $sth->rows() || die "$Script - Unable to get sample data\n";
    my $rowsofdata = $sth->rows() || return;
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        push @bamidrange, $href->{bamid};
        push @nwdids, $href->{expt_sampleid};
    }
}

#==================================================================
# Subroutine:
#   VerifyFile($file)
#
#   To verify integrity of a .cram file, see awk script
#     /net/topmed/incoming/study.reference/current.6.2016/nhlbi.1530.cram.end.region.awk
#   File /net/topmed/incoming/study.reference/study.reference/NWD850567.trunc.test.cram
#   provides a bad example on which the test should fail
#
#   There must exist an index file for this file... and it should be
#   dated after the cram
#
#   Return error msg if failure
#==================================================================
sub VerifyFile {
    my ($file) = @_;
    if (! $file) { return 0; }
    if ($file !~ /\.(\w+)$/) { die "Unable to parse extension for '$file'\n"; }
    my $ext = $1;

    my $index = 'crai';
    if ($ext eq 'bam') { $index = 'bai'; }
    if ($ext eq 'bcf') { $index = 'csi'; }

    #   Check for the index file and that it is younger than the file
    my @s = stat($file);
    if (! @s) {
        if ($file =~ /\.bam$/) { return ''; }    # Missing BAM is OK
        return 'File not found';                 # But not ok otherwise
    }
    my $filetime = $s[9];
    my $filesize = $s[7];
    if (($ext eq 'bam' || $ext eq 'cram') && $filesize < 7288450000) {
        return "File seems too small ($filesize)";
    }

    @s = stat("$file.$index");
    if (! @s) { return "Index missing"; }
    my $indextime = $s[9];
    if ($indextime < $filetime) { return "Index is older"; }

    #   Shall we check if the cram is truncated ?
    if ($opts{checkfiles} && $ext eq 'cram') {
        #my $cmd = "$opts{samtools} view -H $file | awk -f /net/topmed/incoming/study.reference/current.6.2016/nhlbi.1530.cram.end.region.awk";
        #   This command does not behave as Tom claims

        #do ValidateIndex bamid bamorcram
    }
    return '';
}

#==================================================================
# Subroutine:
#   ($center, $dirname) = GetCenterRun($runid)
#
#   Return center and run name for a runid
#==================================================================
sub GetCenterRun {
    my ($runid) = @_;

    my $sth = My_DB::DoSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Run '$runid' unknown, how'd that happen?\n"; }
    my $href = $sth->fetchrow_hashref;
    my $dirname = $href->{dirname};
    $sth = My_DB::DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Center for '$dirname', how'd that happen?\n"; }
    $href = $sth->fetchrow_hashref;
    return ($href->{centername}, $dirname);
}

#==================================================================
# Subroutine:
#   ShowErrors($b, $nwd, $aref)
#
#   Print the array of error messages
#==================================================================
sub ShowErrors {
    my ($b, $nwd, $aref) = @_;
    #print "  Errors for BAMID=$b  NWD=$nwd -\n";
    foreach my $e (@{$aref}) {
        if ($opts{fixer}) {
            my $cmd = $opts{fixer} . ' ' . $b . ' ' . $nwd . ' ' . $e;
            if ($opts{execfixer}) { 
                print "  EXECUTING: $cmd\n";
                system($cmd);
            }
            else { print "  $cmd\n"; }
        }
        else {
            print $b . ' ' . $e . "\n";
        }
    }
}

#==================================================================
# Subroutine:
#   Check_NFSMounts()
#   Special case - check one sample in each run for each center
#   We have to do this because NFS is so flakey. Run this on each host.
#==================================================================
sub Check_NFSMounts {
    my $centersref = GetCenters();          # Only returns one center
    foreach my $cid (keys %{$centersref}) {
        $opts{centername} = $centersref->{$cid};    # Save center name
        my $runsref = GetRuns($cid) || next;
        #   Get one bamid for every run
        foreach my $runid (keys %{$runsref}) {
            my $sql = "SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} " .
                "WHERE runid=$runid AND state_cram=20 ORDER BY RAND() LIMIT 1";
            my $sth = My_DB::DoSQL($sql);
            GetArrayOfSamples($sth);
        }
    }
    #@bamidrange=(44872, 127984);
    #   We now have a list of one sample from each run for all centers
    if ($opts{verbose}) { print "Verify NFS mounts by checking " .
        scalar(@bamidrange) . " samples\n"; }
    my $errs = 0;
    my $success = 0;
    my $h = `hostname`;
    chomp($h);
    foreach my $bamid (@bamidrange) {
        my $p = WhereFile($bamid, 'cram');
        if (! $p) {
            print "Unable to find path for cram file $bamid\n";
            next;
        }
        if (-f $p) { $success++; next; }
        print "Host '$h', failed find file $p\n";
        $errs++;
    }
    if (! $errs) {
	if ($opts{verbose}) { print "NFS mounts on '$h' seem to be okay  (found $success files)\n"; }
    }
    else { print "NFS mounts on '$h' found $success files, failed to find $errs files\n"; }
    exit($errs);
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedcheck.pl - Check all the monitor information for a set of samples

=head1 SYNOPSIS

  topmedcheck.pl -sample 34567 gcefiles       # Check files in GCE for a single sample
  topmedcheck.pl -sample 34567-24577 awsfiles # Check files in AWS for a small set of samples
  topmedcheck.pl -run 2015Sep05 localfiles    # Check localfiles for samples in a run
  topmedcheck.pl -center nygc db              # Check database for samples in a center
  topmedcheck.pl -center nygc backups         # Check backup files still exist
  topmedcheck.pl -center nygc -datayear 1 db  # Check database  for samples in center/year

  topmedcheck.pl -nfsmounts                   # Special case - check for broken NFS mounts
 
=head1 DESCRIPTION

This program provides a means to sanity check all the details for a sample
or set of samples.
It can check that the database columns appear to be correct.
It can check that the files in the local filesystem or remote filesystems are correct.

You must specify either B<-center>, B<-run> or <-sample>.
Of that set of samples, you may further subset the data checked with the other options
like B<-build> etc.

The messages it generates should be easily parsable to you can provide the
output to a shell script to correct the failures.
See the B<-redo> and B<-fixer> options.

=head1 OPTIONS

=over 4

=item B<-build N>

Selects samples from a particular build.

=item B<-center NAME>

Selects samples from a particular center.

=item B<-datayear N>

Selects samples from a particular year.

=item B<-execfixer>

If set this will force the B<-fixer> command to be run as the errors are detected.

=item B<-fixer COMMANDSTRING>

If specified a string will be prepended to error messages.
This provides a convenient way to build shell commands that could be used
to correct the error. See B<-redo> below for more detail.

=item B<-help>

Generates this output.

=item B<-max N>

Reduce the number of samples that would normally be selected to N samples.

=item B<-nfsmounts>

Check all the NFS mounts. This is a special case to ensure the NFS mounts for this host still work.

=item B<-run dirname>

Selects samples from a particular run.

=item B<-piname NAME>

Selects samples from a particular PI.

=item B<-redo>

Short hand for B<-fixer /usr/cluster/topmed/bin/topmed_redo.sh>.
This is a convenient way to build a command which can correct database or data problems.
You might use it like this:

  topmedcheck.pl -redo 87654 | grep redo.sh > /tmp/errors
  .../topmed_redo.sh 87654 Column b38flagstat x must match bamflagstat 0
  bash /tmp/errors

Of course you will carefully check /tmp/errors before running it, but this trick
allows one to capture all the error messages and build a command
(e.g. topmed_redo.sh or your own command) to correct the problem.

=item B<-replace>

Force the cache file to be rebuilt.

=item B<-runs runname>

Selects samples from a run. You maybe specify the dirname or the runid.

=item B<-studyname NAME>

Selects samples from a studyname.

=item B<-verbose>

Provided for developers to see additional information.

=back


=head1 PARAMETERS

The only parameter is used to specify the action to take -- e.g. what to check.
Possibilities are to check the specified samples for:

  db           Check the database entries
  backups      Check the backup file is in GCE
  localfiles   Check local file system files
  awsfiles     Check the finished files in AWS
  gcefiles     Check the finished files in GCE
  b37          Check the remapped B37 files in GCE
  recabfiles   Check the GCE bucket for remapped, but not finished files


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

