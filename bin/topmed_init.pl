#!/usr/bin/perl
###################################################################
#
# Name:	topmed_init.pl
#
# Description:
#   Use this program to check for files that are arriving
#   and initialize the NHLBI TOPMED database
#
# ChangeLog:
#   $Log: topmed_init.pl,v $
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
use My_DB;
use TopMed_Get;
use Getopt::Long;

use POSIX qw(strftime);
use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my $NOTSET    = 0;            # Not set
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    topdir => '/net/topmed/incoming/topmed',
    runcount => 0,
    bamcount => 0,
    bamcountruns => '',
    arrivedsecs => 86400*7,         # If no new bam in a week, stop looking at run
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s run=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] updatedb\n" .
        "Monitor NHLBI data arriving in a directory (default=$opts{topdir}').\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

my $dbh = DBConnect($opts{realm});
my $nowdate = time();
if ($fcn ne 'updatedb') { die "$Script  - Invalid function '$fcn'\n"; }
chdir($opts{topdir}) ||
    die "$Script Unable to CD to '$opts{topdir}': $!\n";

#--------------------------------------------------------------
#   For each center watch for a new run to arrive
#--------------------------------------------------------------
my $centersref = GetCenters();
foreach my $centerid (keys %{$centersref}) {
    my $c = $centersref->{$centerid};
    my $d = $opts{topdir} . '/' . $c;
    if (! chdir($d)) {
        warn "$Script Unable to CD to '$d': $!\n";
        next;
    }
    #   Get all the known batch runs for this center that are not fully arrived
    my $sql = "SELECT runid,dirname,arrived FROM $opts{runs_table} WHERE centerid=$centerid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %knownruns = ();
    my $dir;
    for (my $i=1; $i<=$rowsofdata; $i++) {
        foreach my $href ($sth->fetchrow_hashref) {
            $dir = $href->{dirname};            
            $knownruns{$dir}{runid} = $href->{runid};
            $knownruns{$dir}{arrived} = $href->{arrived};
        }
    }
    #   Get list of all runs for this center
    #   Find new ones and add details to database
    my $dirsref = GetDirs('.');
    my $runid;
    foreach my $d (@{$dirsref}) {
        $runid = $knownruns{$d}{runid} || CreateRun($centerid, $d);
        if (! defined($runid)) {
            if ($opts{verbose}) { warn "$Script RUNID not defined for dir=$d\n"; }
            next;
        }
        if ($opts{run} && $opts{run} ne $d) { next; }
        #   Check if this run has arrived, no need to look at it further
        if (defined($knownruns{$d}{arrived}) && $knownruns{$d}{arrived} eq 'Y') { next; }
        if ($opts{verbose}) { print "Try to add BAMs in '$d' [$runid]\n"; }
        AddBams($runid, $d);
    }
}

$nowdate = strftime('%Y/%m/%d %H:%M', localtime);

if ($opts{runcount}) { print "$nowdate  Added $opts{runcount} runs\n"; }
if ($opts{bamcount}) { print "$nowdate  Added $opts{bamcount} bams from: $opts{bamcountruns}\n"; }
exit;

#==================================================================
# Subroutine:
#   CreateRun - Add details on this run to the database
#   Make sure this directory is runnable. If not, no new runid
#
# Arguments:
#   d - directory (e.g. run name)
#   cid - center id
#
# Returns:
#   runid or undef
#==================================================================
sub CreateRun {
    my ($cid, $d) = @_;
    #   Runs with a magic name are ignored
    if ($d eq 'upload' || $d eq 'slots') { return undef(); }

    #   Try to write in this directory. Can't trust the users
    my $touchfile = '.test';
    if ( ! open(TOUCH, '>' . $touchfile)) {
        if ($opts{verbose}) { print "$Script - Ignoring non-writable directory '$d'\n"; }
        return undef();
    }
    if ( ! print TOUCH 'test if writable') {
        if ($opts{verbose}) { print "$Script - Ignoring non-writable directory '$d'\n"; }
        return undef();
    }
    close(TOUCH);
    unlink($touchfile);

    #   Directory is writable, create SQL record
    my $sql = "INSERT INTO $opts{runs_table} " .
        "(centerid,dirname,comments,bamcount,dateinit) " .
        "VALUES($cid,'$d','',0,'$nowdate')";
    my $sth = DoSQL($sql);
    my $runid = $sth->{mysql_insertid};
    print "$Script - Added run '$d' [ $runid ]\n";
    $opts{runcount}++;
    return $runid;
}

#==================================================================
# Subroutine:
#   AddBams - Add details for bams in this directory to the database
#       In order to know if we've found all the bam files, we look
#       for any BAMs older than the MD5 files. If so, we parse the
#       md5 files again.
#
# Arguments:
#   runid - run id
#   d - directory
#
# Returns:
#   Boolean if any were added or not
#==================================================================
sub AddBams {
    my ($runid, $d) = @_;
    if (! -d $d) {
        print "$Script - Unable to read directory '$d': $!\n";
        return 0;
    }

    #   Get all the known bams for this directory/run   Get NWD name and original name
    my $sql = "SELECT bamname,bamname_orig FROM $opts{bamfiles_table} WHERE runid=$runid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %knownbams = ();
    my %knownorigbam = ();
    for (my $i=1; $i<=$rowsofdata; $i++) {
        foreach my $href ($sth->fetchrow_hashref) {
            $knownbams{$href->{bamname}} = 1;
            if (exists($href->{bamname_orig}) && defined($href->{bamname_orig}) &&
                $href->{bamname_orig} ne '' && $href->{bamname_orig} ne ' ') {
                $knownorigbam{$href->{bamname_orig}} = 1;
            }
        }
    }

    #   Get list of all MD5 files and reprocess them. Don't get excited when
    #   the BAM has already been added to bamfiles. This will not go on forever.
    my @md5files = ();
    if ( -r "$d/Manifest.txt") { push @md5files,'Manifest.txt'; }
    else {
        opendir(my $dh, $d) ||
            die "$Script - Unable to read directory '$d'\n";
        while (readdir $dh) {
            if (! /\.md5$/) { next; }
            push @md5files,$_;
        }
        closedir $dh;
    }
    if (! @md5files) {
        if ($opts{verbose}) { print "$Script - No new MD5 files found in $d\n"; }
        return 0;
    }

    #   There is no consistency what people do here.
    #   Foreach md5 file, get the BAM name and create the bamfiles record
    #   and rename the md5 file so we do not process it again
    my $newbams = 0;
    my ($fn, $checksum);
    my $md5lines = '';
    foreach my $origf (@md5files) {
        my $f = "$d/$origf";
        if (! open(IN, $f)) {
            warn "$Script - Unable to read MD5 file '$f': $!\n";
            next;
        }
        my $badmd5 = 0;
        while (my $l = <IN>) {          # Read md5 checksum file
            ($fn, $checksum) = NormalizeMD5Line($l, $f);
            if (! $checksum) { $badmd5++; next; }
            #   Don't make duplicate records
            if (exists($knownbams{$fn})) { next; }       # Skip known bams
            if (exists($knownorigbam{$fn})) { next; }    # Skip known original bams

            #   New BAM, create database record for it. Actual BAM may not exist yet
            $sql = "INSERT INTO $opts{bamfiles_table} " .
                "(runid,bamname,checksum,dateinit) " .
                "VALUES($runid,'$fn','$checksum', $nowdate)";
            if ($opts{verbose}) { print "SQL=$sql\n"; }
            $sth = DoSQL($sql);
            $newbams++;
        }
        close(IN);
        if ($badmd5) { next; }              # This MD5 was in error
    }

    #   If we added bams, change the bamcount
    if ($newbams) {
        print "$Script - $newbams new bams found in '$d'\n";
        $opts{bamcount} += $newbams;            # Stats for ending msg
        $opts{bamcountruns} .= $d . ' ';
    }

    #   Get number of database records
    $sql = "SELECT bamid FROM $opts{bamfiles_table} WHERE runid=$runid";
    $sth = DoSQL($sql);
    my $numbamrecords = $sth->rows();
    $sql = "UPDATE $opts{runs_table}  SET bamcount=$numbamrecords WHERE runid=$runid";
    $sth = DoSQL($sql);

    #   Last sanity check, see if number of BAM files matches records
    #   This might not always be accurate, but ...
    if ($newbams) {
        my $n = `ls $d/N*.bam $d/N*.cram 2>/dev/null | wc -l`;
        chomp($n);
        if ($n eq $numbamrecords) { print "$Script - Congratulations, # bams = # database records\n"; }
        else {
            print "$Script - Warning, # bams [$n] != # database records [$numbamrecords].  " .
            "If data is incoming, this might be OK\n";
        }
    }

    #   If SOME ALL the sample has been processed as arrived, maybe we do not
    #   need to look at this run any more.
    $sql = "SELECT bamid from $opts{bamfiles_table} WHERE state_arrive!=$NOTSET";
    $sth = DoSQL($sql);
    my $numberarrived = $sth->rows();
    if (! $numberarrived) { return 1; }         # No sample processed, keep looking

    #   At least one sample was marked as arrived, see if there has been little changed
    #   in this directory in a long time, then make this run as 'arrived' so we'll
    #   look at this run again
    my $oldestbamdate = OldestBAM($d);
    if ((time() - $oldestbamdate) > $opts{arrivedsecs}) {
        $sql = "UPDATE $opts{runs_table}  SET arrived='Y' WHERE runid=$runid";
        $sth = DoSQL($sql);
        print "$Script - Run '$d' has finally arrived. Look at it no more\n";
    }
    return 1;
}

#==================================================================
# Subroutine:
#   Get date of oldest BAM in this directory
#   filedate = OldestBAM($d)
#
# Arguments:
#   d - directory to search
#
# Returns:
#   Date of oldest BAM or zero
#==================================================================
sub OldestBAM {
    my ($d) = @_;
    my $oldestbamdate = 0;
    opendir(my $dh, $d) ||
        die "$Script - Unable to read directory '$d'\n";
    while (readdir $dh) {
        if ((! /\.bam$/) && (! /\.cram$/)) { next; }
        my @stats = stat("$d/$_");
        if ($oldestbamdate < $stats[9]) { $oldestbamdate = $stats[9]; }
    }
    closedir $dh;
    return $oldestbamdate;
}

#==================================================================
# Subroutine:
#   NormalizeMD5Line - Return filename and checksum
#
# Arguments:
#   l - line from an md5 file (well, what they pretend is an md5 file)
#   f - name of file (for error msgs only)
#
# Returns:
#   array of filename and checksum or NULL
#==================================================================
sub NormalizeMD5Line {
    my ($l, $f) = @_;
    my @retvals = ();
    if ($l =~ /^#/) { return @retvals; }

    my ($checksum, $fn) = split(' ',$l);
    #   Do sanity check trying to guess the format of their md5 file
    if ($checksum =~ /\./) { ($fn, $checksum) = ($checksum, $fn); }
    if (length($checksum) != 32) {
        print "$Script - Invalid checksum '$checksum' in '$f'. Line: $l";
        return @retvals;
    }
    if ($fn =~ /\//) { $fn = basename($fn); }   # Path should not be qualified, but ...
    if ($fn !~ /bam$/ && $fn !~ /cram$/) {      # File better be a bam or cram
        #   Only generate error message for true errors, not dumb ones :-(
        if ($fn !~ /bai$/) { print "$Script - Invalid BAM name '$fn' in '$f'. Line: $l"; }
        return @retvals;
    }
    @retvals = ($fn, $checksum);
    return @retvals;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_init.pl - check for files that are arriving and initialize the database

=head1 SYNOPSIS

  topmed_init.pl updatedb

=head1 DESCRIPTION

This program monitors directories for incoming data and then
creates records for new run directories and/or new BAMs
(based on new MD5 files being discovered).

=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<topmed>.

=item B<-runs NAME>

Specifies a specific run on which to run the action.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<updatedb>

This directs the program to monitor the B<-dir> directory for changes
and update the database table specified by B<-realm>.

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

