package Topmed::DB::Schema::ResultSet::Bamfile;

use base qw(DBIx::Class::ResultSet);

use Topmed::Constants qw(:states);

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

sub completed_for_build {
  my ($self, %params) = @_;

  my $state_col_map = {
    b37 => 'state_b37',
    b38 => 'state_b38',
  };

  return unless exists $state_col_map->{$params{build}};

  return $self->search(
    {
      $state_col_map->{$params{build}} => $COMPLETED,
    }, {
      order_by => 'expt_sampleid',
    }
  );
}

1;
