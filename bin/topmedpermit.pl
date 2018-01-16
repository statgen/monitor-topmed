#!/usr/bin/perl -I /usr/cluster/boehnke/lib/perl5
###################################################################
#
# Name: topmedpermit.pl
#
# Description:
#   Use this program to determine if an action is permitted
#   This program is tightly coupled with topmed/html/index.php
#
#   This was originally part of topmedcmd.pl, but had to be broken
#   into a separate program because of all the dependencies
#   that topmedcmd.pl had.
#
#   This has the minimum dependencies possible, even though
#   that means there is duplicated code here.
#
# ChangeLog:
#   $Log: topmedpermit.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use FindBin qw($Bin $Script);
use Getopt::Long;
use DBIx::Connector;

use POSIX qw(strftime);

our ($DBC, $DBH);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    centers_table => 'centers',
    runs_table => 'runs',
    permissions_table => 'permissions',
    squeuedata => '/run/shm/squeue.results',
    verbose => 0,
);
my %VALIDOPS = (                    # Used for permit
    all => 1,
    verify => 1,
    gcebackup => 1,
    qplot => 1,
    cram => 1,
    sb37 => 1,
    gcepush => 1,
    gcepull => 1,
    sb38 => 1,
    bcf => 1,
    gcecopy => 1,
    awscopy => 1,
    sexpt => 1,
    sorig => 1,
);

Getopt::Long::GetOptions( \%opts,qw(
    help conf=s verbose
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] test action bamid [host]\n" .
        "$Script [options] remove id\n" .
        "$Script [options] add action datayear center run\n" .
        "\n" .
        "Return code indicates if an action is permitted or not.\n" .
        "You may also add or remove a permit-rule.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}

my $fcn = shift @ARGV;
DBConnect($opts{realm});

if ($fcn eq 'test')    { Test(@ARGV); exit; }
if ($fcn eq 'remove')  { Remove(@ARGV); exit; }
if ($fcn eq 'add')     { Add(@ARGV); exit; }

die "$Script - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   Test ($op, $bamid, $run)
#
#   Is operation allowed
#       e.g. topmedpermit.pl permit test qplot 4955 [HOST]
#   Can fail because too many jobs for host or database rule#
#
#   Returns zero if an operation for a center/run is allowed, else failure
#==================================================================
sub Test {
    my ($op, $bamid, $h) = @_;

    #   Sometimes we just do not want to check for limits
    if ($ENV{IGNORE_PERMIT}) { exit; }

    #   To start, get the the system and host maximums permitted
    #       e.g. action= bcf sysmax= 350 hostmax= 35
    my $systemmax = 0;
    my $hostmax = 0;
    my ($in);
    if ($h && open($in, $opts{squeuedata})) {
        while (<$in>) {
            if (/JOBID/) { last; }
            my @c = split(' ', $_);
            if ($c[1] eq $op && $c[2] eq 'sysmax=' && $c[4] eq 'hostmax=') {
                $systemmax = $c[3];
                $hostmax = $c[5];
            }
        }
        #print "op=$op host=$h systemmax=$systemmax hostmax=$hostmax\n";
        #   Data we read is like this:
        #   3896779  topmed-working topmed 91398-bcf  topmed PD 0:00 1 (None)   topmed4
        #   3798393  topmed-working topmed 104073-bcf topmed  R 3:24:57 1 topmed7   topmed7
        my $queued = 0;                 # How many are queued on this host
        my $running = 0;                # How many are running on this host
        my $opcount = 0;                # How many of these operations are queued/running anywhere
        while (<$in>) {
            my @c = split(' ', $_);
            if ($c[1] ne 'topmed-working') { next; }
            if ($c[3] !~ /\d+-$op/) { next; }
            $opcount++;
            if ($c[9] ne $h) { next; }
            if ($c[5] eq 'R') { $running++; }
            if ($c[5] eq 'PD') { $queued++; }
        }
        close($in);
        #print "op=$op host=$h queued=$queued running=$running opcount=$opcount\n";
        
        #   See if this would overload the host or the system
        my $now = strftime('%Y/%m/%d %H:%M', localtime);
        if ($opcount >= $systemmax) {
            print "$now Too many $op jobs [$opcount] running/queued system wide [$systemmax]\n";
            exit(5);                        # Don't ask again, special return code
        }
        my $n = $queued + $running;
        if ($n >= $hostmax) {
            print "$now Too many $op jobs [$n] running/queued on $h [$hostmax]\n";
            exit(4);                        # Operation not permitted, special return code
        }
    }

    #   This is a small table, read it all
    my $sth = DoSQL("SELECT runid,centerid,datayear,operation FROM $opts{permissions_table}");
    my $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        #   Get runid and centerid for this bam
        my ($centerid, $runid, $datayear) = GetBamidInfo($bamid);
        for (my $i=1; $i<=$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref;
            if ($centerid eq $href->{centerid} || $href->{centerid} eq '0') {
                if ($runid eq $href->{runid} || $href->{runid} eq '0') {
                    if ($datayear eq $href->{datayear} || $href->{datayear} eq '0') {
                        if ($op eq $href->{operation} || $href->{operation} eq 'all') {
                            print "Operation $op blocked by $opts{permissions_table} table\n";
                            exit(4);        # Operation not permitted, special return code
                        }
                    }
                }
            }
        }
    }
    #sleep(2);                               # Try to avoid possible race in SLURM ??
    exit;                                   # Permitted
}

#==================================================================
# Subroutine:
#   Remove($id)
#
#   Enable a permission by deleting a database entry
#   e.g. topmedcmd.pl permit remove id
#==================================================================
sub Remove {
    my ($id) = @_;
    my $sql = "SELECT * FROM $opts{permissions_table} WHERE id=$id";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "Permission '$id' does not exist\n"; }
    my $href = $sth->fetchrow_hashref;
    my $s = "$href->{datayear} /$href->{centername} / $href->{dirname} / $href->{operation}";
    $sql = "DELETE FROM $opts{permissions_table} WHERE id='$id'";
    $sth = DoSQL($sql);
    if (! $sth) { die "Failed to remove permission '$fcn' for 'id=$id $s'\nSQL=$sql"; }
    print "Deleted permission control for '$s'\n";
    exit;
}

#==================================================================
# Subroutine:
#   Add($op, $datayear, $center, $run)
#
#   Disable a permission by adding to the database entry
#   e.g. ... add [qplot [broad [2015sep18]]]
#==================================================================
sub Add {
    my ($op, $datayear, $center, $run) = @_;
    my ($sql, $sth, $href, $rowsofdata);

    if (! $op)       { $op = 'all'; }           # Set defaults
    if (! $datayear) { $datayear = 'all'; }
    if (! $center)   { $center = 'all'; }
    if (! $run)      { $run = 'all'; }
    my ($centerid, $runid) = (0, 0);

    if (! exists($VALIDOPS{$op})) { die "Operation '$op' is not known\n"; }

    #   Verify op and run and center setting runid and centerid
    if ($center ne 'all') {
        $sql = "SELECT centerid FROM $opts{centers_table} WHERE centername='$center'";
        $sth = DoSQL($sql);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "Center '$center' is not known\n"; }
        $href = $sth->fetchrow_hashref;
        $centerid = $href->{centerid};
    }

    if ($run ne 'all') {
        $sql = "SELECT runid,dirname FROM $opts{runs_table} WHERE dirname='$run' OR runid=$run";
        $sth = DoSQL($sql);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "Run '$run' is not known\n"; }
        $href = $sth->fetchrow_hashref;
        $runid = $href->{runid};
        $run = $href->{dirname};
    }

    if ($datayear ne 'all') {
        if ($datayear !~ /^\d+$/) { die "datayear '$datayear' must be a number\n"; }
    }
    else { $datayear = 0; }

    $sql = "INSERT INTO $opts{permissions_table} " .
        "(datayear,centername,dirname,centerid,runid,operation) " .
        "VALUES ($datayear, '$center', '$run', $centerid, $runid, '$op')";
    $sth = DoSQL($sql);
    if (! $sth) { die "Failed to add permission '$fcn' for '$op'\nSQL=$sql"; }
    print "Added permission '$datayear / $center / $runid / $op'\n";
    exit;
}

#==================================================================
# Subroutine:
#   ($centerid, $runid, $datayear) = GetBamidInfo($bamid)
#
#   Return id for the center and run for a bamid
#==================================================================
sub GetBamidInfo {
    my ($bamid) = @_;

    if ($bamid !~ /^\d+$/) { return (0,0); }     # No bamid, no ids
    my $sth = DoSQL("SELECT runid,datayear FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};
    my $datayear = $href->{datayear};
    $sth = DoSQL("SELECT centerid FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "BAM '$bamid' has no center?  How'd that happen?\n"; }
    $href = $sth->fetchrow_hashref;
    return ($href->{centerid}, $runid, $datayear);
}

#==================================================================
# Subroutine:
#   DBConnect($realm)
#
#   Connect to our database using realm '$realm'. Return a DB handle.
#   Get the connection information from DBIx::Connector
#   Fully qualified realm file may be provided
#   Subsequent SQL errors will generate an error message.
#==================================================================
sub DBConnect {
    my ($realm) = @_;
    if (! $realm) { warn "DBConnect: No REALM provided\n";  return 0; }
    #   Get the connection information FROM DBIx::Connector
    #   Fully qualified realm file may be provided
    if ($realm =~ /^(\/.+)\/([^\/]+)$/) {
        my $d = $1;
        $realm = $2;
        $DBC = new DBIx::Connector(-realm => $realm, -connection_dir => $d,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    else {
        $DBC = new DBIx::Connector(-realm => $realm,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    $DBH = $DBC->connect();
    return $DBH;
}

#==================================================================
# Subroutine:
#   DoSQL - Execute an SQL command
#
# Arguments:
#   sql - String of SQL to run
#   die - boolean if we should die on error
#
# Returns:
#   SQL handle for subsequent MySQL actions
#   Does not return if error detected
#==================================================================
sub DoSQL {
    my ($sql, $die) = @_;
    if (! defined($die)) { $die = 1; }
    if ($main::opts{verbose} > 1) { warn "DEBUG: SQL=$sql\n"; }
    my $sth = $DBH->prepare($sql);
    if (! $die) { $DBH->{RaiseError} = 0; $DBH->{PrintError} = 0; }
    $sth->execute();
    if ($DBI::err) {
        if (! $die) { return 0; }
        die "SQL failure: $DBI::errstr\n  SQL=$sql\n";
    }
    return $sth;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedpermit.pl - Control if one may submit a job for an action

=head1 SYNOPSIS

  #  Stop qplot job submissions for a particular run
  topmedpermit.pl add all qplot broad 2015oct18
  #  Stop qplot job submissions for a particular center
  topmedpermit.pl add all qplot broad
  #  Stop qplot job submissions for a particular datayear
  topmedpermit.pl add 3 qplot 

  topmedpermit.pl remove 12          # Remove a permit control
  topmedpermit.pl test cram 4567     # Test if we should submit a cram job for this sampleid
  topmedpermit.pl test cram 4567 topmed5   # Test if submit okay for topmed5

=head1 OPTIONS

=over 4

=item B<-conf PATH>

Specifies the path to a topmedthrottle configuration file

=item B<-help>

Generates this output.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

B<add/remove enable/disable action center run>
Use this to control the database which allows one enable or disable topmed actions
(e.g. backup, verify etc) for a center or run.
Use B<all> for all datayears, centers or all runs or all actions.

B<test action bamid host>
Use this to test if an action (e.g. backup, verify etc) may be submitted 
for a particular bam. Returns 4 for not permitted.

=head1 EXIT

When testing, a return code of 4 means the task is not permitted.
A return code of 5 means no more of these are permitted.
When not testing an error returns a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut
