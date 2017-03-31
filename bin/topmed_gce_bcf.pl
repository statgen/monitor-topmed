#!/usr/bin/env perl

use FindBin;
use lib (
  qq($FindBin::Bin/../local/lib/perl5),
  qq($FindBin::Bin/../lib/perl5),
);

use Topmed::Base qw(cmds);
use Topmed::Constants qw(:all);
use Topmed::DB;

GetOptions(
  'limit|l=i' => \(my $limit   = 1),
  'dry-run|n' => \(my $dryrun  = undef),
  'verbose|v' => \(my $verbose = undef),
  'help|h'    => \&HelpMessage,
) or HelpMessage(1);

my $schema = Topmed::DB->new();
my $samples = $schema->resultset('Bamfile')->find_gce_bcf_samples;

for my $sample ($samples->slice(0, ($limit - 1))) {
  my $nwdid    = $sample->expt_sampleid;
  my $gs_uri   = $GOOGLE_BUCKETS{bcf} . $nwdid;
  my $cram_uri = $gs_uri . $SLASH . $nwdid . '.recab.cram';
  my $bcf_uri  = $gs_uri . $SLASH . $nwdid . '.bcf';
  my $pipeline = qq{$FindBin::Bin/../etc/pipelines/bcf.yaml};
  my @args     = ();

  push @args, q{run};
  push @args, qq{--pipeline-file $pipeline};
  push @args, qq{--logging $gs_uri};
  push @args, qq{--inputs INPUT_FILE=$cram_uri,NWDID=$nwdid};
  push @args, qq{--outputs OUTPUT_FILE=$bcf_uri};
  push @args, sprintf q{--labels nwdid=%s,stage=b38-bcf}, lc($nwdid);
  push @args, '2>&1';

  my $cmd = join($SPACE, ($GOOGLE_PIPELINE_CMD, @args));

  if ($dryrun) {
    say $cmd;
  } else {

    # Example output:
    # Running [operations/EPbpyo6wKxjSu7ndsdvjEyC6tq6izhcqD3Byb2R1Y3Rpb25RdWV1ZQ].
    chomp(my $output = capture($cmd));
    (my $operation_id = $output) =~ s|Running \[operations/([^\]]+)\]\.|$1|g;

    $sample->update(
      {
        state_gce38bcf => $SUBMITTED,
        gce38bcf_opid  => $operation_id,
      }
    );

    say "SAMPLE: $nwdid Operation: $operation_id" if $verbose;
  }
}

__END__

=head1 NAME

topmed_gce_bcf.pl - Submit Google Genomics Pipeline jobs for BCF generation

=head1 SYNOPSIS

topmed_gce_bcf.pl [options]

  Options:

    -l, --limit   Limit the number of jobs submitted to Google
    -n, --dry-run Output what would be run but do not actually run it
    -v, --verbose Output what actions are happening
    -h, --help    Usage info

=head1 DESCRIPTION

This script will submit Google Genomics Pipeline jobs for BCF generation.

Example pipelines command:

    gcloud alpha genomics pipelines run \
      --pipeline-file ./bcf.yaml \
      --logging gs://topmed-bcf/${sample} \
      --inputs INPUT_FILE=${input},NWDID=${sample} \
      --outputs OUTPUT_FILE=${output} \
      --labels nwdid=${sample,,},stage=b38-bcf

=head1 ENVIRONMENT VARIABLES

=head1 VERSION

0.1

=cut
