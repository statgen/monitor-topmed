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
    b37 => 'me.state_b37',
    b38 => 'me.state_b38',
  };

  return unless exists $state_col_map->{$params{build}};

  return $self->search(
    {
      $state_col_map->{$params{build}} => $COMPLETED,
      -and => [
        {'qc_results.pct_freemix'     => {'<' => 3}},
        {'qc_results.pct_genome_dp10' => {'>', 95}},
      ],
    }, {
      join     => 'qc_results',
      order_by => 'expt_sampleid',
    }
  );
}

sub find_gce_bcf_samples {
  return shift->search(
    {
      state_verify    => $COMPLETED,
      state_gce38push => $COMPLETED,
      state_gce38bcf  => $NOTSET,
    }, {
      order_by => 'RAND()',
    }
  );
}

sub check_gce_bcf_status {
  return shift->search(
    {
      state_gce38bcf_push => $COMPLETED,
      state_gce38bcf_pull => $NOTSET,
      state_gce38bcf      => $SUBMITTED,
      gce38bcf_opid       => {-not => undef},
      expt_sampleid       => {-not => undef},
    }, {
      order_by => 'RAND()',
    }
  );
}

1;
