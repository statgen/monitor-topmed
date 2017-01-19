package Topmed::Base;

use base qw(Import::Base);

our @IMPORT_MODULES = (
  qw(
    FindBin
    Data::Dumper
    Modern::Perl
    IO::All
    YASF
    Path::Class
  ),
  'autodie' => [qw(:all)],
  'English' => [qw(-no_match_vars)],
  'Cwd'     => [qw(abs_path)],
);

our %IMPORT_BUNDLES = (
  files => [
    'File::Path' => [qw(make_path remove_tree)],
  ],
);

1;
