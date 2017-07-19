#!/usr/bin/env perl

use FindBin;
use lib (qq($FindBin::Bin/../lib/perl5), qq($FindBin::Bin/../local/lib/perl5),);

use Test::Most qw(no_plan);
use Data::Dumper;

use Topmed::GCE::Pipelines;

# TODO - build pre-align test
{
  my $cmd = 'gcloud alpha genomics pipelines run --logging gs://topmed-mapping/logs/NWD12346 --pipeline-file /home/schelcj/src/monitor-topmed/t/../etc/pipelines/mapping/pre-align.yaml --inputs INPUT_FILE=gs://topmed-incoming/NWD123456/NWD123456.src.cram,NWDID=NWD12346 --outputs OUTPUT_PATH=gs://topmed-mapping/pre-align/NWD123456/ --labels build=b38,nwdid=nwd123456,stage=pre-align,year=y3';

  my $pipeline = Topmed::GCE::Pipelines->new(sample_id => 'NWD12346', stage => 'pre-align');

  ok($pipeline->add_input(NWDID => $pipeline->sample_id), 'added nwdid input param');
  ok($pipeline->add_input(INPUT_FILE => 'gs://topmed-incoming/NWD123456/NWD123456.src.cram'), 'added input_file param');

  ok($pipeline->add_output(OUTPUT_PATH => 'gs://topmed-mapping/pre-align/NWD123456/'), 'added output_path param');

  ok($pipeline->add_label(NWDID => 'NWD123456'),      'added nwdid label');
  ok($pipeline->add_label(STAGE => $pipeline->stage), 'added stage label');
  ok($pipeline->add_label(BUILD => 'b38'),            'added build label');
  ok($pipeline->add_label(YEAR  => 'y3'),             'added year label');

  is($pipeline->_cmd, $cmd, 'command matches');
}

# TODO - build align test
{
  my $cmd = 'gcloud alpha genomics pipelines run --logging gs://topmed-mapping/logs/NWD12346 --pipeline-file /home/schelcj/src/monitor-topmed/t/../etc/pipelines/mapping/align.yaml --inputs INPUT_FILE=gs://topmed-mapping/pre-align/NWD123456/NWD123456.0.foo.fastq.gz,READ_GROUP=Foo --outputs FASTQ_READS_FILE=gs://topmed-mapping/align/NWD123456.0.foo.fastq.gz.reads,FLAGSTAT_FILE=gs://topmed-mapping/align/NWD123456.0.foo.cram.flagstat,OK_FILE=gs://topmed-mapping/align/NWD123456.0.foo.cram.ok,OUTPUT_FILE=gs://topmed-mapping/align/NWD123456.0.foo.cram --labels build=b38,stage=align,year=y3';

  my $pipeline = Topmed::GCE::Pipelines->new(sample_id => 'NWD12346', stage => 'align');

  ok($pipeline->add_input(INPUT_FILE => 'gs://topmed-mapping/pre-align/NWD123456/NWD123456.0.foo.fastq.gz'), 'added input_file param');
  ok($pipeline->add_input(READ_GROUP => 'Foo'), 'added read_group param'); # TODO - this will likely be a problem with real data

  ok($pipeline->add_output(OUTPUT_FILE => 'gs://topmed-mapping/align/NWD123456.0.foo.cram'),    'added output_file param');
  ok($pipeline->add_output(OK_FILE     => 'gs://topmed-mapping/align/NWD123456.0.foo.cram.ok'), 'added ok_file param');
  ok($pipeline->add_output(FASTQ_READS_FILE => 'gs://topmed-mapping/align/NWD123456.0.foo.fastq.gz.reads'), 'added fastq_reads_file param');
  ok($pipeline->add_output(FLAGSTAT_FILE => 'gs://topmed-mapping/align/NWD123456.0.foo.cram.flagstat'), 'added flagstat_file param');

  ok($pipeline->add_label(STAGE => $pipeline->stage), 'added stage label');
  ok($pipeline->add_label(BUILD => 'b38'),            'added build label');
  ok($pipeline->add_label(YEAR  => 'y3'),             'added year label');

  is($pipeline->_cmd, $cmd, 'command matches');
}

# TODO - build post-align test
{
  my $cmd = 'gcloud alpha genomics pipelines run --logging gs://topmed-mapping/logs/NWD12346 --pipeline-file /home/schelcj/src/monitor-topmed/t/../etc/pipelines/mapping/post-align.yaml --inputs INPUT_PATH=gs://topmed-mapping/align/NWD123456/*.cram,NWDID=NWD12346 --outputs BCF=gs://topmed-mapping/results/NWD123456/NWD123456.bcf,CRAI=gs://topmed-mapping/results/NWD123456/NWD123456.recab.cram.crai,CRAM=gs://topmed-mapping/results/NWD123456/NWD123456.recab.cram,CSG=gs://topmed-mapping/results/NWD123456/NWD123456.bcf.csi,FLAGSTAT=gs://topmed-mapping/results/NWD123456/NWD123456.recab.cram.flagstat,METRICS=gs://topmed-mapping/results/NWD123456/NWD123456.dedup_lowmem.metrics --labels build=b38,stage=post-align,year=y3';
  my $pipeline = Topmed::GCE::Pipelines->new(sample_id => 'NWD12346', stage => 'post-align');

  ok($pipeline->add_input(NWDID => $pipeline->sample_id), 'added nwdid param');
  ok($pipeline->add_input(INPUT_PATH => 'gs://topmed-mapping/align/NWD123456/*.cram'), 'added input_path param');

  ok($pipeline->add_output(CRAM => 'gs://topmed-mapping/results/NWD123456/NWD123456.recab.cram'), 'added cram param');
  ok($pipeline->add_output(CRAI => 'gs://topmed-mapping/results/NWD123456/NWD123456.recab.cram.crai'), 'added crai param');
  ok($pipeline->add_output(BCF  => 'gs://topmed-mapping/results/NWD123456/NWD123456.bcf'), 'added bcf param');
  ok($pipeline->add_output(CSG  => 'gs://topmed-mapping/results/NWD123456/NWD123456.bcf.csi'), 'added csi param');
  ok($pipeline->add_output(METRICS => 'gs://topmed-mapping/results/NWD123456/NWD123456.dedup_lowmem.metrics'), 'added metrics param');
  ok($pipeline->add_output(FLAGSTAT => 'gs://topmed-mapping/results/NWD123456/NWD123456.recab.cram.flagstat'), 'added flagstat param');

  ok($pipeline->add_label(STAGE => $pipeline->stage), 'added stage label');
  ok($pipeline->add_label(BUILD => 'b38'),            'added build label');
  ok($pipeline->add_label(YEAR  => 'y3'),             'added year label');

  is($pipeline->_cmd, $cmd, 'command matches');
}
