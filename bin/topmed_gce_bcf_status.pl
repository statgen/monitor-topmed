#!/usr/bin/env perl

use FindBin;
use lib (qq($FindBin::Bin/../local/lib/perl5), qq($FindBin::Bin/../lib/perl5),);

use Topmed::Base qw(cmds);
use Topmed::Constants qw(:all);
use Topmed::DB;

use List::MoreUtils qw(part);
use Parallel::ForkManager;

GetOptions(
  'procs|p=i' => \(my $procs   = 10),
  'dry-run|n' => \(my $dryrun  = undef),
  'verbose|v' => \(my $verbose = undef),
  'help|h'    => \&HelpMessage,
  )
  or HelpMessage(1);

my $i      = 0;
my $pm     = Parallel::ForkManager->new($procs);
my $schema = Topmed::DB->new();
my @parts  = part {$i++ % $procs} $schema->resultset('Bamfile')->check_gce_bcf_status()->all();

for my $part (@parts) {
  $pm->start and next;
  for my $sample (@{$part}) {
    my $cmd    = join($SPACE, ($GOOGLE_OPERATION_CMD, 'describe', $sample->gce38bcf_opid));
    my $status = YAML::Syck::Load(io->pipe($cmd)->chomp->all());

    unless ($status->{done} eq 'true') {
      say $sample->expt_sampleid . ' : NOT COMPLETE' if $verbose;
      next;
    }

    if (exists $status->{error}) {
      say $sample->expt_sampleid . ' : FAILED' if $verbose;
      $sample->update({state_gce38bcf => $FAILED}) unless $dryrun;
      next;
    }

    # XXX - best guess that the operation succeeded:
    #
    #     - done is true
    #     - no error attribute
    #     - has an event list
    #     - last event was an ok state
    #
    if (scalar @{$status->{metadata}->{events}}) {
      if ($status->{metadata}->{events}->[-1]->{description} eq 'ok') {
        say $sample->expt_sampleid . ' : COMPLETED' if $verbose;
        $sample->update({state_gce38bcf => $COMPLETED}) unless $dryrun;
      }
    }
  }
  $pm->finish;
}
$pm->wait_all_children;

__END__

=head1 NAME

topmed_gce_bcf_status.pl - Check on the operation status of bcf generation jobs.

=head1 SYNOPSIS

topmed_gce_bcf_status.pl [options]

  Options:

    -p, --procs   Number of concurrent samples to check (default: 10)
    -v, --verbose Output what actions are happening
    -h, --help    Usage info

=head1 DESCRIPTION

This script will check the status of Google Genomics Pipeline jobs for BCF generation.

Example pipelines command:

    gcloud alpha genomics operations descripbe <operation_id>

=head1 ENVIRONMENT VARIABLES

=head1 VERSION

0.1

=cut

