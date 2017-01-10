package Topmed::DB;

use base qw(Topmed::DB::Schema);

use DBIx::Connector;

sub new {
  my $ctx   = DBIx::Connector->new(realm => 'topmed');
  my $realm = $ctx->read_realm_file();

  return __PACKAGE__->connect(
    qq{dbi:mysql:database=$realm->{DATABASE};$realm->{SERVER}},
    $realm->{USER},
    $realm->{PASS}
  );
}

1;
