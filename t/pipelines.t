#!/usr/bin/env perl

use FindBin;
use lib (
  qq($FindBin::Bin/../lib/perl5),
  qq($FindBin::Bin/../local/lib/perl5),
);

use Test::Most qw(no_plan);
use Data::Dumper;

use Topmed::GCE::Pipelines;

my $pipeline = Topmed::GCE::Pipelines->new(sample_id => 'NWD12346', stage => 'pre-align');

isa_ok($pipeline, 'Topmed::GCE::Pipelines');
ok($pipeline->add_input(foo => 'bar'), 'added an input parameter');
ok($pipeline->add_output(foo => 'bar'), 'added an output parameter');
ok($pipeline->add_label(foo => 'bar'), 'added an label parameter');
ok($pipeline->add_label(baz => 'bar'), 'added an label parameter');

print STDERR Dumper($pipeline->_cmd);
