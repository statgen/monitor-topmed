package Topmed::DB;

use base qw(Topmed::DB::Schema);

use DBIx::Connector;

sub new {
  #  Not only hardcode the realm, provide the path too
  my $ctx = DBIx::Connector->new(realm => 'topmed', connection_dir => '/usr/cluster/topmed/etc/.db_connections');
  my $realm = $ctx->read_realm_file();

  return __PACKAGE__->connect(
    qq{dbi:mysql:database=$realm->{DATABASE};$realm->{SERVER}},
    $realm->{USER},
    $realm->{PASS}
  );
}

1;
