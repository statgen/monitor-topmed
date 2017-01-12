package Topmed::DB::Schema::ResultSet::Bamfile;

use base qw(DBIx::Class::ResultSet);

sub find_by_nwdid {
  my ($self, $nwdid) = @_;
  return $self->find({expt_sampleid => $nwdid});
}

1;
