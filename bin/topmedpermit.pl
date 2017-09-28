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
#   that topmedcmd.pl now has.
#
#   This has the minimum dependencies possible.
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
    help verbose
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] permit test action bamid\n" .
        "$Script [options] permit remove id\n" .
        "$Script [options] permit add action datayear center run\n" .
        "Returns boolean if an action is permitted or not\n" .
        "or can be used to add or remove a permit-rule\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}

my $fcn = shift @ARGV;
DBConnect($opts{realm});
if ($fcn eq 'permit')    { Permit(@ARGV); exit; }
die "$Script - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   Permit ($fcn, $op, $center, $run)
#
#   Controls a database of enabled or disabled operations for a center/run
#
#   Permit ('test', $op, $bamid)
#
#   Returns a boolean if an operation for a center/run is allowed
#==================================================================
sub Permit {
    my ($fcn) = @_;
    my ($sql, $sth, $href, $rowsofdata);
    shift(@_);

    #   Test is some operation is allowed
    #   topmedcmd.pl permit test qlpot 4955
    if ($fcn eq 'test') {
        my ($op, $bamid) = @_;

        #   Get runid and centerid for this bam
        my ($centerid, $runid, $datayear) = GetBamidInfo($bamid);

        #   This is a small table, read it all
        $sth = DoSQL("SELECT runid,centerid,datayear,operation FROM $opts{permissions_table}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { exit(1); }             # Table empty, permitted
        for (my $i=1; $i<=$rowsofdata; $i++) {
            $href = $sth->fetchrow_hashref;
            if ($centerid eq $href->{centerid} || $href->{centerid} eq '0') {
                if ($runid eq $href->{runid} || $href->{runid} eq '0') {
                    if ($datayear eq $href->{datayear} || $href->{datayear} eq '0') {
                        if ($op eq $href->{operation} || $href->{operation} eq 'all') {
                            exit;                   # Operation not permitted
                        }
                    }
                }
            }
        }
        exit(1);                                    # TRUE, permitted  
    }

    #   Enable a permission by deleting a database entry
    #   e.g. topmedcmd.pl permit remove id
    if ($fcn eq 'remove') {
        my ($id) = @_;
        $sql = "SELECT * FROM $opts{permissions_table} WHERE id=$id";
        $sth = DoSQL($sql);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "Permission '$id' does not exist\n"; }
        $href = $sth->fetchrow_hashref;
        my $s = "$href->{datayear} /$href->{centername} / $href->{dirname} / $href->{operation}";
        $sql = "DELETE FROM $opts{permissions_table} WHERE id='$id'";
        $sth = DoSQL($sql);
        if (! $sth) { die "Failed to remove permission '$fcn' for 'id=$id $s'\nSQL=$sql"; }
        print "Deleted permission control for '$s'\n";
        exit;
    }

    if ($fcn ne 'add') { die "Unknown function '$fcn'\n"; }

    #   Disable a permission by adding to the database entry
    #   e.g. # topmedcmd.pl permit add [qplot [broad [2015sep18]]]
    my ($op, $datayear, $center, $run) = @_;
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

topmedpermit.pl - Determine if an action is allowed

=head1 SYNOPSIS

  #  Stop qplot job submissions for a particular run
  topmedpermit.pl permit add all qplot broad 2015oct18
  #  Stop qplot job submissions for a particular center
  topmedpermit.pl permit add all qplot broad
  #  Stop qplot job submissions for a particular datayear
  topmedpermit.pl permit add 3 qplot 

  topmedpermit.pl permit remove 12          # Remove a permit control
  topmedpermit.pl permit test cram 4567     # Test if we should submit a cram job for this sampleid

=head1 OPTIONS

No options are supported.

=head1 PARAMETERS

B<permit enable/disable action center run>
Use this to control the database which allows one enable or disable topmed actions
(e.g. backup, verify etc) for a center or run.
Use B<all> for all datayears, centers or all runs or all actions.

B<permit test action bamid>
Use this to test if an action (e.g. backup, verify etc) may be submitted 
for a particular bam.

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut
