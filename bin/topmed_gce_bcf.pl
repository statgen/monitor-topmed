#!/usr/bin/env perl

use FindBin;
use lib (
  qq($FindBin::Bin/../local/lib/perl5),
  qq($FindBin::Bin/../lib/perl5),
);

use Topmed::Base;
use Topmed::Constants qw(:all);
use Topmed::DB;

my $schema = Topmed::DB->new();

for my $sample ($schema->resultset('Bamfiles')->find_gce_bcf_samples) {

}
