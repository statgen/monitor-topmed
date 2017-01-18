#!/usr/bin/env perl

use FindBin;
use lib (qq($FindBin::Bin/../local/lib/perl5), qq($FindBin::Bin/../lib/perl5),);

use Topmed::Base qw(cmds);
use Topmed::Constants qw(:states);
use Topmed::DB;

my $schema       = Topmed::DB->new();
my $sample_rs    = $schema->resultset('Bamfile')->find_gce_uploads();
my $flagstat_ptn = YASF->new('gsutil ls {uri}/{nwdid}.recab.cram.flagstat 2> /dev/null');

for my $sample ($sample_rs->all) {
  my $cmd = $flagstat_ptn % {uri => $sample->gce_recab_uri, nwdid => $sample->expt_sampleid};
  capture([0..1], $cmd);

  unless ($EXITVAL) {
    # TODO - set state_gce38pull to $REQUESTED
    #
    # $sample->update({state_gce38pull => $REQUESTED});
    #
    say 'SAMPLE : ' . $sample->expt_sampleid;
  }
}
