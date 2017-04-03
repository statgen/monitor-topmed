#!/usr/bin/env perl

use FindBin qw($Bin);
use lib (qq($Bin/../lib/perl5), qq($Bin/../local/lib/perl5));

use DBIx::Connector;
use DBIx::Class::Schema::Loader qw(make_schema_at);

my $ctx = DBIx::Connector->new(realm => 'topmed');
my $realm = $ctx->read_realm_file();

make_schema_at(
  'Topmed::DB::Schema', {
    debug          => 0,
    dump_directory => qq($Bin/../lib/perl5),
    constraint     => qr/^(?:actions|bamfiles|bamfiles_actions|builds|centers|freezes|qc_results|runs|states)$/,
    components     => [qw(InflateColumn::DateTime)],
  },
  [qq{dbi:mysql:database=$realm->{DATABASE};$realm->{SERVER}}, $realm->{USER}, $realm->{PASS}]
);
