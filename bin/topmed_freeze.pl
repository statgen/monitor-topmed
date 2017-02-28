#!/usr/bin/env perl

use FindBin;
use lib (qq($FindBin::Bin/../lib/perl5), qq($FindBin::Bin/../local/lib/perl5));

use Topmed::Base;
use Topmed::Constants;
use Topmed::DB;

use Getopt::Long qw(HelpMessage);
use Class::CSV;

GetOptions(
  'build|b=s'  => \(my $build  = 'b38'),
  'output|o=s' => \(my $output = $DASH),
  )
  or HelpMessage(1);

unless ($build =~ /^b37|b38$/) {
  say 'Invalid build specified!';
  exit 1;
}

if ($output ne $DASH and not -w $output) {
  say 'Invalid output file';
  exit 1;
}

my $headers = [qw(nwdid pi center run cram)];
my %csv_params = (fields => $headers);

if ($output eq $DASH) {
  $csv_params{filehandle} = io->stdout->tie;
} else {
  $csv_params{filename} = $output;
}

my $csv     = Class::CSV->new(%csv_params);
my $schema  = Topmed::DB->new();
my $samples = $schema->resultset('Bamfile')->completed_for_build(build => $build);

$csv->add_line({map {$_ => $_} @{$headers}});

for my $sample ($samples->all) {
  my $meth = qq(${build}_mapped_path);
  $csv->add_line(
    {
      nwdid  => $sample->nwdid,
      pi     => $sample->piname,
      center => $sample->center,
      run    => $sample->runname,
      cram   => $sample->$meth,
    }
  );
}

__END__

=head1 NAME

topmed_freeze.pl - Generate a CSV list of remapped crams

=head1 SYNOPSIS

topmed_freeze.pl [optoin]

  Options:

    -b, --build   Build to export in format b# (e.g. b38|b37) [default: b38]
    -o, --output  File to write csv output to or - for STDOUT [default: STDOUT]

=head1 DESCRIPTION

This script will generate a CSV of sample info for a given build.
Columns that are exported:

  NWDID, PI, CENTER, RUN, CRAM

=head1 ENVIRONMENT VARIABLES

=head1 VERSION

0.1

=cut
