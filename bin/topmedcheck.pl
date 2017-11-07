#!/usr/bin/perl
###################################################################
#
# Name: topmedcheck.pl
#
# Description:
#    Do a serious check of data for a set of bamids
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
use Cwd qw(realpath abs_path);
use My_DB;
use Topmed::Constants qw(:states);
use Topmed::Get;

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
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    centers_table => 'centers',
    runs_table => 'runs',
    gceuri => 'gs://topmed-bcf',
    topmedpath => '/usr/cluster/topmed/bin/topmedpath.pl',
    gcels => 'gsutil ls',   # Add $opts{gceuri}/NWDxxxxxxx
    redotool => '/usr/cluster/topmed/bin/topmed_redo.sh',
    samtools => '/usr/cluster/bin/samtools',
    #db => 1,                # Always check the database values
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help db localfiles awsfiles gcefiles checkfiles b37files
    register checkfiles fixer=s execfixer random max=i redo verbose
    center=s run=s piname=s studyname=s datayear=i
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($opts{help} || ($#ARGV<0 && (! $opts{run}) && (! $opts{center}))) {
    my $m = "$Script [options]";
    warn "$m bamid|bamid-bamid|nwdid\n" .
        " or\n" .
        "$m -subsetkey value  (e.g. -center, -runs, -piname, -studyname, -datayear)\n" .
        "\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
if (! ($opts{db} || $opts{localfiles} || $opts{awsfiles} || $opts{gcefiles} || $opts{b37files})) {
    die "$Script - Specify the type of check:  -db -localfiles -awsfiles or -gcefiles\n";
}
if ($opts{redo}) { $opts{fixer} = $opts{redotool}; }

DBConnect($opts{realm});

#--------------------------------------------------------------
#   Get list of bamids to check
#--------------------------------------------------------------
my @bamidrange = ();
my @nwdids = ();
my $morewhere = '';
if ($opts{piname})    { $morewhere .= " AND piname='$opts{piname}'"; }
if ($opts{studyname}) { $morewhere .= " AND studyname='$opts{studyname}'"; }
if ($opts{datayear})  { $morewhere .= " AND datayear=$opts{datayear}"; }
if ($opts{b37files})  { $morewhere .= " AND state_b37=20"; }
if ($opts{random})    { $morewhere .= ' ORDER BY RAND()'; }
if ($opts{max})       { $morewhere .= " LIMIT $opts{max}"; }

if (@ARGV) {                # Maybe one or a range of bamids specified
    my $bamid = shift @ARGV;
    my @input = ();
    if ($bamid =~/^(\d+)\-(\d+)/) { @input = $1 .. $2; }
    else { push @input, $bamid; }
    for my $b (@input) {
        my $sql = "SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} WHERE ";
        if ($b =~ /^NWD/) { $sql .= "expt_sampleid='$b'"; }
        else { $sql .= "bamid=$b"; }
        $sql .= $morewhere;
        my $sth = DoSQL($sql);
        GetArrayOfSamples($sth);
    }
}
else {                      # Range specified by center or run
    if ($opts{center}) {
        my $centersref = GetCenters();          # Only returns one center
        foreach my $cid (keys %{$centersref}) {
            $opts{centername} = $centersref->{$cid};    # Save center name
            my $runsref = GetRuns($cid) || next;
            #   Get bamid for every run
            foreach my $runid (keys %{$runsref}) {
                my $sql = "SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} " .
                    "WHERE runid=$runid AND state_cram=20" . $morewhere;
                my $sth = DoSQL($sql);
                GetArrayOfSamples($sth);
            }
        }
    }
    #   Get all bamids for this run
    if ($opts{run}) {
        if ($opts{run}=~/\D+/) {               # Get runid from dirname
            my $sql = "SELECT runid FROM $opts{runs_table} WHERE dirname='$opts{run}'";
            my $sth = DoSQL($sql);
            if (! $sth) { die "$Script - Unknown run '$opts{run}'\n"; }
            my $href = $sth->fetchrow_hashref;
            $opts{run} = $href->{runid};
        }
        #   Get bamid for this run
        my $sql = "SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} " .
            "WHERE state_cram=20 AND runid=$opts{run}" . $morewhere;
        my $sth = DoSQL($sql);
        GetArrayOfSamples($sth);
    }
}
if (! @bamidrange) { die "$Script - no Samples found\n"; }
if ($opts{max}) { splice(@bamidrange,$opts{max}); } # We might have selected too many samples          
print "Checking " . scalar(@bamidrange) . " samples\n";

#--------------------------------------------------------------
#   Do varying checks over a set of bamids
#--------------------------------------------------------------
if ($opts{db}) { 
    print "Verify database values\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_DB($bamidrange[$i], $nwdids[$i]);
    }
}
if ($opts{localfiles}) {
    print "Verify local files\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_LocalFiles($bamidrange[$i], $nwdids[$i]);
    }
}
if ($opts{gcefiles}) {
    print "Verify GCE files\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_GCEFiles($bamidrange[$i], $nwdids[$i]);
    }
}
if ($opts{awsfiles}) {
    print "Verify Amazon files\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_AWSFiles($bamidrange[$i], $nwdids[$i]);
    }
}
if ($opts{b37files}) {
    print "Verify the local files for remapped b37 - patience!\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_B37Files($bamidrange[$i], $nwdids[$i]);
    }
}
if ($opts{checkfiles}) {
    print "Verify the local files (md5 and flagstat) - patience!\n";
    for (my $i=0; $i<=$#bamidrange; $i++) {
        Check_Files($bamidrange[$i], $nwdids[$i]);
    }
}
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
    my $sth = DoSQL($sql);
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
    if ($href->{$c} !~ /^1|2|3$/) {
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
    my $filecmd = $opts{topmedpath} . " wherefile $bamid ";

    #   Get some data for this sample
    my $sql = "SELECT datayear,state_b37,state_b38 FROM $opts{bamfiles_table} WHERE bamid=$bamid";
    my $sth = DoSQL($sql);
    if (! $sth) { die "$Script - Unable to get data for $sql\n"; }
    my $href = $sth->fetchrow_hashref;
    my $year = $href->{datayear};

    #   Check original file (could be bam or cram)
    my $f = `$filecmd bam`;
    chomp($f);
    my $e = VerifyFile($f);
    if ($e) { push @error,"BAM cram/index is invalid: $e: $f"; }
    else {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid BAM OK: $f\n"; }
    }
    
    #   Check cram, either original file or converted from bam
    $f = `$filecmd cram`;
    chomp($f);
    $e = VerifyFile($f);
    if ($e) { push @error,"CRAM index is invalid: $e: $f"; }
    else {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid CRAM OK: $f\n"; }
    }

    #   QPlot files should exist
    $f = `$filecmd qcresults`;
    if (! $f) { die "Unable to determine qcresults file for $bamid\n"; }
    chomp($f);
    if ( -f $f ) {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid QCRESULTS OK: $f\n"; }
    }
    else { push @error,"QCRESULTS missing: $f"; }

    #   B37 crams should exist for ALL year 1 samples
    #   Some other year samples were also mapped
    if ($year eq '1' && $href->{state_b37} != $COMPLETED) {
        push @error,"Warning: Year $year sample was not mapped for build B37";
    }
    if ($href->{state_b37} == $COMPLETED) {
        $f = `$filecmd b37`;
        chomp($f);
        my $e = VerifyFile($f);
        if ($e) { push @error,"B37 cram/index is invalid: $e: $f"; }
        else {
            if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid B38 OK: $f\n"; }
        }
    }

    #   B38 crams must exist for all samples
    if ($href->{state_b38} == $COMPLETED) {
        $f = `$filecmd b38`;
        chomp($f);
        my $e = VerifyFile($f);
        if ($e) { push @error,"B38 cram/index is invalid: $e: $f"; }
        else {
            if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid B38 OK: $f\n"; }
        }
    }

    #   BCF must exist for all samples
    if (! defined($href->{state_gce38bcf})) { $href->{state_gce38bcf} = 0; }
    if ($href->{state_gce38bcf} == $COMPLETED) {
        $f = `$filecmd bcf`;
        chomp($f);
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
#   Check_GCEFiles()
#   Verify files in Google Cloud Environment
#==================================================================
sub Check_GCEFiles {
    my ($bamid, $nwdid) = @_;
    my @error = ();

    #   Do not bother if we have not copied data to GCE
    my $sql = "SELECT state_gce38copy FROM $opts{bamfiles_table} WHERE bamid=$bamid";
    my $sth = DoSQL($sql);
    my $href = $sth->fetchrow_hashref;
    if ($href->{state_gce38copy} != $COMPLETED) {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid not copied to GCE yet\n"; }
        return;
    }

    #   All data is supposed to be in GCE
    my $cmd = $opts{gcels} . ' ' . $opts{gceuri} . "/$nwdid/";
    my %gcefiles = ();
    my %requiredextensions = ( 'bcf' => 1, 'bcf.csi' => 1, 'recab.cram' => 1,
        'recab.cram.crai' => 1, 'recab.cram.md5' => 1, 'recab.cram.flagstat' => 1 );
    my $rmleng = length($opts{gceuri}) + length($nwdid) + 2;
    my @lines = split("\n", `$cmd`);
    foreach my $l (@lines) {
        chomp($l);
        if ((!defined($l)) || (! $l)) { next; }
        if ($l =~ /INFO/) { next; }     # gsutil had to retry
        my $f = substr($l,$rmleng);     # Isolate just the filename
        if (! defined($f)) {
            warn "Undefined here $bamid $nwdid\n" . @lines . "\n\n";
            next;
        }
        if ($f !~ /^$nwdid\.(\S+)/) { push @error,"GCE file is invalid: Unknown file: '$f'"; }
        else {
            my $x = $1;                 # Extension of file
            if (exists($requiredextensions{$x})) { $requiredextensions{$x} = 0; }
            else { push @error,"GCE NWDID file is invalid: Invalid extention $x: $f"; }
        }
    }
    foreach my $x (keys %requiredextensions) {
        if ($requiredextensions{$x}) { push @error,"GCE $x file missing: $nwdid.$x"; }
    }

    #   If there were problems, show messages
    if (@error) { ShowErrors($bamid, $nwdid, \@error); }
    else {
        if ($opts{verbose}) { print "  BAMID=$bamid NWD=$nwdid - GCE files all exist\n"; }
    }
}

#==================================================================
# Subroutine:
#   Check_AWSFiles()
#   Verify files in Amazon Environment
#==================================================================
sub Check_AWSFiles {
    my ($bamid, $nwdid) = @_;
    my @error = ();
    my $filecmd = $opts{topmedpath} . " wherefile $bamid ";

}

#==================================================================
# Subroutine:
#   Check_B37Files()
#   Verify all the B37 files exist
#==================================================================
sub Check_B37Files {
    my ($bamid, $nwdid) = @_;
    my @error = ();
    my $filecmd = $opts{topmedpath} . " wherefile $bamid ";

    #   Check local file 
    my $f = `$filecmd b37`;
    chomp($f);
    my $e = VerifyFile($f);
    if ($e) { push @error,"B37 cram/index is invalid: $e: $f"; }

    #   If there were problems, show messages
    if (@error) { ShowErrors($bamid, $nwdid, \@error); }
    else {
        if ($opts{verbose}) { print "BAMID=$bamid NWD=$nwdid B37 CRAM OK: $f\n"; }
    }

}

#==================================================================
# Subroutine:
#   Check_Files()
#   Calculate the md5sum and flagstat value for files in this sample
#==================================================================
sub Check_Files {
    my ($bamid, $nwdid) = @_;
    my @error = ();
    my $filecmd = $opts{topmedpath} . " wherefile $bamid ";

}

#==================================================================
# Subroutine:
#   GetArrayOfSamples($sth)      Follows $sth = DoSQL($sql)
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
#   ShowErrors($b, $nwd, $aref)
#
#   Print the array of error messages
#==================================================================
sub ShowErrors {
    my ($b, $nwd, $aref) = @_;
    print "  Errors for BAMID=$b  NWD=$nwd -\n";
    foreach my $e (@{$aref}) {
        if ($opts{fixer}) {
            my $cmd = $opts{fixer} . ' ' . $b . ' ' . $e;
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
my $a=$b;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedcheck.pl - Check all the monitor information for a set of samples

=head1 SYNOPSIS

  topmedcheck.pl 34567                      # Check single sample
  topmedcheck.pl 34567-24577                # Check a small set of samples
  topmedcheck.pl -run 2015Sep05             # Check samples in a run
  topmedcheck.pl -center nygc               # Check samples in a center
  topmedcheck.pl -center nygc -pi MESA      # Check samples in a center for a PI
  topmedcheck.pl -center nygc -datayear 1   # Check samples in a center for a year

  topmedcheck.pl -local 34567               # Check local files for a sample
  topmedcheck.pl -aws 34567                 # Check AWS files for a sample

 
=head1 DESCRIPTION

This program provides a means to sanity check all the details for a sample
or set of samples.
It can check that the database columns appear to be correct.
It can check that the files in the local filesystem or remote filesystems are correct.

The messages it generates should be easily parsable to you can provide the
output to a shell script to correct the failures.
See the B<-redo> and B<-fixer> options.


=head1 OPTIONS

=over 4

=item B<-awsfiles>

Specifies that files in the AWS filesystem should be checked.

=item B<-center NAME>

Selects samples from a particular center.

=item B<-checkfiles>

Specifies that more thorough checking should be done on local CRAM and index files.
This will slow down the check quite a bit.

=item B<-datayear N>

Selects samples from a particular year.

=item B<-fixer COMMANDSTRING>

If specified a string will be prepended to error messages.
This provides a convenient way to build shell commands that could be used
to correct the error. See B<-redo> below for more detail.

=item B<-gcefiles>

Specifies that files in the GCE filesystem should be checked.

=item B<-help>

Generates this output.

=item B<-localfiles>

Specifies that files in the local filesystem should be checked.

=item B<-max N>

Reduce the number of samples that would normally be selected to N samples.

=item B<-register>

Will generate a message when a sample is not registered in the NCBI database.

=item B<-piname NAME>

Selects samples from a particular PI.

=item B<-random>

Randomly select data to be processed.

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

=item B<-runs runname>

Selects samples from a run. You maybe specify the dirname or the runid.

=item B<-studyname NAME>

Selects samples from a studyname.

=item B<-verbose>

Provided for developers to see additional information.

=item B<-verifyfiles>

Specifies that the local files should be verified as completely as possible.
This will mean a md5sum and flagstat value will be calculated for each
source file and compared to the database value.
B<This is extremely disk intensive and will take a long time.> 

=back


=head1 PARAMETERS

Parameters to this program can be specified in several ways:

  * Single sampleid (NWDID or bamid).
  * As a range of bamids (e.g.  50-60)

If no parameter is provided, you must provide some other way to select the
samples, B<-center>, B<-run>, B<-piname>, B<-studyname> or B<-datayear>.


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

