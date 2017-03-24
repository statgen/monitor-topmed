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
    YAML::Syck
  ),
  'autodie'      => [qw(:all)],
  'English'      => [qw(-no_match_vars)],
  'Cwd'          => [qw(abs_path)],
  'Getopt::Long' => [qw(HelpMessage)],
);

our %IMPORT_BUNDLES = (
  cmds => [
    'IPC::System::Simple' => [qw(run capture EXIT_ANY $EXITVAL)],
  ],
  files => [
    'File::Path' => [qw(make_path remove_tree)],
  ],
);

1;
