package My_DB;
#==================================================================
# My_DB.pm
#   Common DB functions for the Perl programs
#   See DBIx::Connector for details defining the database
#==================================================================
use strict;
use warnings;
use DBIx::Connector;

use vars qw(@EXPORT_OK @EXPORT @ISA $DBC $DBH);
use Exporter;

@ISA = qw(Exporter);
@EXPORT_OK = qw(DBConnect DoSQL SQL_Last_Insert PersistDBConnect PersistDoSQL);
@EXPORT    = qw(DBConnect DoSQL SQL_Last_Insert PersistDBConnect PersistDoSQL);

our ($DBC, $DBH);

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
#   PersistDBConnect($realm)
#
#   Connect to our database using realm '$realm'. Return a DB handle.
#   Get the connection information from DBIx::Connector
#   Fully qualified realm file may be provided
#   If the connection fails, keep trying for quite some time.
#   Subsequent SQL errors will NOT generate an error message.
#==================================================================
sub PersistDBConnect {
    my ($realm) = @_;
    if (! $realm) { die "PersistDBConnect: No REALM provided\n"; }

    #   Get the connection information FROM DBIx::Connector
    #   Fully qualified realm file may be provided
    if ($realm =~ /^(\/.+)\/([^\/]+)$/) {
        my $d = $1;
        $realm = $2;
        $DBC = new DBIx::Connector(-realm => $realm, -connection_dir => $d,
            -max_wait => '20min', -max_interval => '20sec', -callback => \&PersistCallBack,
            -dbi_options => {RaiseError => 0, PrintError => 0});
    }
    else {
        $DBC = new DBIx::Connector(-realm => $realm,
            -max_wait => '20min', -max_interval => '20sec', -callback => \&PersistCallBack,
            -dbi_options => {RaiseError => 0, PrintError => 0});
    }
    if ($DBC) { $DBH = $DBC->connect(); }
    return $DBH;
}

#==================================================================
# Subroutine:
#   PersistCallBack - Invoked by DBI when a connect fails
#   See DBIx::Connector for more details
#
# Arguments:
#   Hash as described in 
#   realm - realm name DBIx::Connector
#
#==================================================================
sub PersistCallBack {
    warn "PersistDBConnect: Datbase connection failed for realm '$_[1]', " ,
        "waiting $_[5] secs to retry. Emsg=$_[9]\n"; 
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
    if (defined($main::opts{verbose}) && $main::opts{verbose} > 1) { warn "DEBUG: SQL=$sql\n"; }
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
# Subroutine:
#   PersistDoSQL($realm, $sql)
#   Connect to our database using realm '$realm' and Execute SQL.
#   If this fails for some environment error (e.g. lost connection)
#   persist for a long time.
#   Our intent is that this SQL command simply must complete
#   unless the SQL itself is invalid.
#
# Arguments:
#   realm - realm giving DB connection information
#   sql - String of SQL to run
#
# Returns:
#   SQL handle for subsequent MySQL actions
#   Does not return if error detected
#==================================================================
sub PersistDoSQL {
    my ($realm, $sql) = @_;
    my $emsg;
    if (defined($main::opts{verbose}) && $main::opts{verbose} > 1) { warn "DEBUG: SQL=$sql\n"; }

    for (my $maxsleep=120; $maxsleep>0; $maxsleep-=5) {   
        if (! $DBH) { PersistDBConnect($realm); }
        my $sth = $DBH->prepare($sql);
        if ($sth) {
            $sth->execute();
            if (! $DBI::err) { return $sth; }   # Success !
        }
        #   Failed. Keep trying for certain kinds of errors, else fail
        $emsg = $DBI::errstr;
        if ($emsg !~ /Lost connection/) {
            die "SQL failure: $emsg\n  SQL=$sql\n";
        }
        undef($DBH);                    # Force reconnect to server
        if (defined($main::opts{verbose}) && $main::opts{verbose} > 1) { warn "PersistDoSQL: $emsg, wait and retry\n"; }
        sleep(5);
    }
     die "SQL failure: $emsg\n  SQL=$sql\n";
}

#==================================================================
# Subroutine:
#   SQL_Last_Insert - returns the autoincrement of the last INSERT
#
# Arguments:
#   sth - STH handle returned by DoSQL
#
# Returns:
#   autoincrement id
#==================================================================
sub SQL_Last_Insert {
    my ($sth) = @_;
    return $sth->{mysql_insertid};
}

1;

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

My_DB.pm

=head1 SYNOPSIS

  use My_DB;
  
=head1 DESCRIPTION

This provides common DB convenience functions.

See B<perldoc DBIx::Connector> for details defining the database.

=head1 Functions

=over 4

=item B<DBConnect($realm)>

Connect to the database. Uses a REALM file as defined in DBIx::Connector.
This can be a fully qualified path to a realm file or just the name
of the realm file (e.g. abc.realm).

=item B<DoSQL($sql, $die)

Executes an SQL command and will normally die unless $die is FALSE.

=back

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2010 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut
