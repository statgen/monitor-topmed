#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name:	monitor.pl
#
# Description:
#   Use this program to monitor data arriving in a directory.
#   Database tables are updated with changes
#   If the directory contains the file 'donotmonitor', we ignore it.
#
#   Requires MIME:Lite    (apt-get install libmime-lite-perl)
#
# ChangeLog:
#   $Log: monitor.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use File::Basename;
use Cwd;
use Getopt::Long;
use DBIx::Connector;

my($me, $mepath, $mesuffix) = fileparse($0, '\.pl');
(my $version = '$Revision: 1.6 $ ') =~ tr/[0-9].//cd;
if (substr($mepath,0,1) ne '/') {
    my $d = getcwd();
    chdir($mepath) ||
        die "$me$mesuffix Unable to CD to '$d': $!\n";
    $mepath = getcwd();
    chdir($d);
}
my $dbc;
my %MONITORDIRS = ();               # Directories to SQL id

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my %opts = (
    moncmd => '/usr/cluster/monitor/bin/moncmd.pl',
    realm => '/usr/cluster/monitor/etc/.db_connections/monitor',
    table => 'status',
    dir => '/exports/CoreDump',
    timeup => '6h',
    waitsecs => 3600,
    verbose => 0,
    donotmonitor => 'do_not_monitor_me',    # If this exists, we ignore the directory
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose daemon waitsecs=i dir=s timeup=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$me$mesuffix [options] monitor\n" .
        "\nVersion $version\n" .
        "Monitor data arriving in a directory (default=$opts{dir}').\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;
if ($opts{timeup} =~ /^(\d+)m/i) { $opts{timeup} = 60*$1; }
if ($opts{timeup} =~ /^(\d+)h/i) { $opts{timeup} = 60*60*$1; }

my $dbh = DBConnect($opts{realm});
if ($opts{verbose}) { print "$me$mesuffix Version $version, realm=$opts{realm}\n"; }

#--------------------------------------------------------------
#   Monitor begins here
#--------------------------------------------------------------
if ($fcn ne 'monitor') { die "$me$mesuffix  - Invalid function '$fcn'\n"; }
chdir($opts{dir}) ||
    die "$me$mesuffix Unable to CD to '$opts{dir}': $!\n";

#--------------------------------------------------------------
#   Get list of directories and notice any new ones
#   If a new directory is found, update database with details
#   Rinse and repeat
#--------------------------------------------------------------
my $timetoend = time();
if ($opts{daemon}) { $timetoend += $opts{timeup}; }
else { $opts{waitsecs} = 1; }
my $newdir = 0;
do {
    my $dirsref = GetDirs('.');
    foreach my $d (@{$dirsref}) {
        if (Monitor($d)) { $newdir++; }
    }
    sleep($opts{waitsecs});
} while ( time() < $timetoend);

#--------------------------------------------------------------
#   All done,  finish up
#--------------------------------------------------------------
if ($opts{verbose}) { print "$me$mesuffix terminating, $newdir new directories found\n"; }
exit;

#==================================================================
# Subroutine:
#   DBConnect($realm)
#
#   Connect to our database using realm '$realm'. Return a DB handle.
#   Get the connection information from DBIx::Connector
#   Fully qualified realm file may be provided
#==================================================================
sub DBConnect {
    my ($realm) = @_;
    if (! $realm) { return 0; }
    #   Get the connection information from DBIx::Connector
    #   Fully qualified realm file may be provided
    if ($realm =~ /^(\/.+)\/([^\/]+)$/) {
        my $d = $1;
        $realm = $2;
        $dbc = new DBIx::Connector(-realm => $realm, -connection_dir => $d,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    else {
        $dbc = new DBIx::Connector(-realm => $realm,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    return $dbc->connect();
}

#==================================================================
# Subroutine:
#   GetDirs - Get list of directories to check
#
# Arguments:
#   dirname
#
# Returns:
#   Reference to array of dir names
#==================================================================
sub GetDirs {
    my ($d) = @_;

    opendir(DIR, $d) ||
        die "Unable to read directory '$d': $!";
    my @dirlist = grep { (/^\w/ && -d "$d/$_") } readdir(DIR);
    closedir DIR;
    return \@dirlist;
}

#==================================================================
# Subroutine:
#   Monitor - Check this directory for any changes
#
# Arguments:
#   d - directory to monitor
#
# Returns:
#   Boolean if this was new or not
#==================================================================
sub Monitor {
    my ($d) = @_;

    #   Perhaps we have seen this and it is done?
    if ($MONITORDIRS{$d}) { return 0; }     # Already seen this
    if (-f "$d/$opts{donotmonitor}") { return 0; }  # Dir marked as do not monitor
    my $sth = DoSQL("SELECT * FROM $opts{table} WHERE dirname='$d'");
    my $rowsofdata = $sth->rows();
    if ($rowsofdata) { return 0; }          # Already seen this
    #   If this does not exist, create it so we can do UPDATES only
    my @finfo = stat($d);
    my $t = time();
    DoSQL("INSERT INTO $opts{table} (dirname,status,size,dateinit,datedirname,_lastchange) " .
        "VALUES('$d','init',0,'$t','" . $finfo[9] . "','$t')");
    my $cmd = "$opts{moncmd} incoming $d started";
    system($cmd) &&
        warn "Unable to update the monitor events database. Cmd=$cmd\n";           
    $MONITORDIRS{$d} =  1;                  # Remember we have seen this before
    return 1;
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
    if ($opts{verbose}) { warn "DEBUG: SQL=$sql\n"; }
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    if ($DBI::err) {
        if (! $die) { return 0; }
        die "$me$mesuffix SQL failure: $DBI::errstr\n  SQL=$sql\n";
    }
    return $sth;
}

#==================================================================
# Subroutine:
#   FmtTime - Return current date and time
#
# Arguments:
#   t - optionaltime since epoch
#
# Returns:
#   string
#==================================================================
sub FmtTime {
    my ($t) = @_;
    if (! defined($t)) { $t = time(); return $t; }
    my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my ($sec, $min, $hr, $mday, $month, $year, $wday, $yday, $isdst) = localtime($t);
    $year = 1900 + $year;
    $month++;
    return sprintf("$year/%02d/%02d %02d:%02d:%02d", $month, $mday, $hr, $min, $sec);
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

monitor.pl - Monitor

=head1 SYNOPSIS

  monitor.pl monitor            # Check directories, runs once
  monitor.pl -daemon monitor    # Check directories, runs repeatedly

=head1 DESCRIPTION

This program monitors directories for incoming data and then
updates a database with various status values.
The program B<moncmd.pl> supports simple commands to set key elements of the database.

The typical mode is to run this in 'daemon mode' by launching
it with a crontab entry, perhaps like this:

  01 06 * * * /usr/local/bin/monitor.pl -daemon monitor 2>&1 /dev/null

This will run for some time, checking the status of the directory specified
in by the B<-dir> option and updating the status in a database every
once in a while.

The queue of tasks is kept in a MySQL database.
See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-daemon>

Specifies the script should not exit after checking the directory.
It will wait for B<waitsecs> seconds and then relook at the directory again.
The program exits only after running for B<timeup>.

=item B<-dir>

Specifies the directory to be monitored.
This defaults to /data/CoreDump

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<webcd>.

=item B<-timeup nnn[MH]>

Specifies how long this program will remain running before exitting.
The number may have 'M' (minutes) or 'H' (hours) appended.
This defaults to B<6H>.
This is only used if B<daemon> is specified.

=item B<-verbose>

Provided for developers to see additional information.

=item B<-waitsecs nnn>

Specifies how many seconds this program will wait before checking the database queue.
This defaults to B<300>.
This is only used if B<daemon> is specified.

=back

=head1 PARAMETERS

=over 4

=item B<monitor>

This directs the program to monitor the B<-dir> directory for changes
and update the database table specified by B<-realm>.

=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2010 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

