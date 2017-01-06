#!/usr/bin/env perl

use FindBin qw($Bin);
use lib (qq($Bin/../lib/perl5), qq($Bin/../local/lib/perl5));

use DBIx::Connector;
use DBIx::Class::Schema::Loader qw(make_schema_at);

$ENV{CAG_ILLEGAL_ACCESSOR_NAME_OK} = 1;

my $ctx = DBIx::Connector->new(realm => 'topmed');
my $realm = $ctx->read_realm_file();

make_schema_at(
  'Topmed::DB::Schema', {
    debug          => 1,
    dump_directory => qq($Bin/../lib/perl5),
    components     => [qw(InflateColumn::DateTime)],
  },
  [qq{dbi:mysql:database=$realm->{DATABASE};$realm->{SERVER}}, $realm->{USER}, $realm->{PASS}]
);
