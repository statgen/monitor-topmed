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
@EXPORT_OK = qw(DBConnect DoSQL);
@EXPORT    = qw(DBConnect DoSQL);

our ($DBC, $DBH);

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
    {
        if (! $die) { $DBH->{RaiseError} = 0; $DBH->{PrintError} = 0; }
        $sth->execute();
    }
    if ($DBI::err) {
        if (! $die) { return 0; }
        die "SQL failure: $DBI::errstr\n  SQL=$sql\n";
    }
    return $sth;
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
