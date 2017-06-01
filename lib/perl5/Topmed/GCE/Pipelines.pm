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

has 'stage' => (
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
  is      => 'ro',
  isa     => 'Str',
  lazy    => 1,
  builder => '_build_logging',
);

has 'inputs' => (
  traits  => ['Hash'],
  is      => 'rw',
  isa     => 'HashRef',
  handles => {
    add_input   => 'set',
    input_pairs => 'kv',
  },
);

has 'outputs' => (
  traits  => ['Hash'],
  is      => 'rw',
  isa     => 'HashRef',
  handles => {
    add_output   => 'set',
    output_pairs => 'kv',
  },
);

has 'labels' => (
  traits  => ['Hash'],
  is      => 'rw',
  isa     => 'HashRef',
  handles => {
    add_label   => 'set',
    label_pairs => 'kv',
  },
);

around 'add_label' => sub {
  my ($orig, $self, %params) = @_;

  for (keys %params) {
    unless ($params{$_} =~ /[a-z]([-a-z0-9]*[a-z0-9])?/i) {
      die "invalid label value: $params{$_}";
    }
  }

  return $self->$orig(%params);
};

sub _build_pipeline {
  my ($self, $arg) = @_;

  my $file = qq{$FindBin::Bin/../etc/pipelines/mapping/} . $self->stage . '.yaml';
  die 'pipeline yaml file does not exist' unless -e $file;    # TODO - throw exception, maybe

  return $file;
}

sub _cmd {
  my $self = shift;

  my @args = ($GOOGLE_PIPELINE_CMD, 'run');
  push @args, '--logging ' . $self->logging;
  push @args, '--pipeline-file ' . $self->pipeline;
  push @args, '--inputs ' . $self->inputs_stringify;
  push @args, '--outputs ' . $self->outputs_stringify;
  push @args, '--labels ' . $self->labels_stringify;

  return join($SPACE, @args);
}

sub _build_logging {
  return qq($GOOGLE_BUCKETS{mapping}/logs/) . shift->sample_id;
}

sub _stringify {
  my ($self, @pairs) = @_;
  return join($COMMA, map {$_->[0] . $EQUAL . $_->[1]} sort {$a->[0] cmp $b->[0]} @pairs);
}

sub inputs_stringify {
  my $self = shift;
  return $self->_stringify($self->input_pairs);
}

sub outputs_stringify {
  my $self = shift;
  return $self->_stringify($self->output_pairs);
}

sub labels_stringify {
  my $self = shift;
  return lc($self->_stringify($self->label_pairs));
}

sub run {
  my $self = shift;

  my $cmd = $self->_cmd;
  return;

  chomp(my $output = capture($cmd));
  (my $operation_id = $output) =~ s|Running \[operations/([^\]]+)\]\.|$1|g;

  return $operation_id;
}

__PACKAGE__->meta->make_immutable;

1;
