package Topmed::DB::Schema::ResultSet::Bamfile;

use base qw(DBIx::Class::ResultSet);

sub find_by_nwdid {
  my ($self, $nwdid) = @_;
  return $self->find({expt_sampleid => $nwdid});
}

sub find_gce_uploads {
  return shift->search(
    {
      state_gce38push => 20,
      state_gce38pull => 0,
    }
  );
}

1;
