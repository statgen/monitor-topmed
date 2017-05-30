package Topmed::GCE::Pipelines;

use Topmed::Base;
use Topmed::Constants qw(:all);

use Moose;
use namespace::autoclean;

has 'sample_id' => (
  is       => 'ro',
  isa      => 'Str',
  required => 1,
);

has 'pipeline' => (
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_pipeline',
);

has 'logging' => (
  is      => 'rw',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_logging',
);

has 'inputs' => (
  is      => 'rw',
  isa     => 'ArrayRef[HashRef]',
  handles => {
    add_input   => 'push',
    join_inputs => 'join',
  },
);

has 'outputs' => (
  is      => 'rw',
  isa     => 'ArrayRef[HashRef]',
  handles => {
    add_output   => 'push',
    join_outputs => 'join',
  },
);

has 'labels' => (
  is      => 'rw',
  isa     => 'ArrayRef[HashRef]',
  handles => {
    add_label   => 'push',
    join_labels => 'join',
  },
);

around 'join_labels' => sub {
  my $orig = shift;
  my $self = shift;

  return lc($self->$orig(@_));
};

sub _build_pipeline {
  my ($self, $arg) = @_;

  my $file = qq{$FindBin::Bin/../etc/pipelines/mapping/${arg}.yaml};
  die 'pipeline yaml file does not exist'; # TODO - exception

  return $file;
}

sub _build_logging {
  # TODO - not sure about where this should live yet
  return 'gs://topmed-mapping/logs/' . shift->sample_id;
}

sub run {
}

__PACKAGE__->meta->make_immutable;

1;
