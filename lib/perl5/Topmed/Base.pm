package Topmed::Base;

use base qw(Import::Base);

our @IMPORT_MODULES = (
  qw(
    FindBin
    Data::Dumper
    Modern::Perl
    IO::All
  ),
  'autodie' => [qw(:all)],
  'English' => [qw(-no_match_vars)],
);

our %IMPORT_BUNDLES = (
);

1;
