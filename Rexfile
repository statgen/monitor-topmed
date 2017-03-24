# vim: set ft=perl
#
use Rex -base;
use Rex::Commands::SCM;
use Rex::Commands::Sync;

set 'install_dir' => '/var/tmp/monitor';
set 'build_dir'   => '/net/topmed/working/build/monitor-topmed';

set repository => "master", url => 'https://github.com/statgen/monitor-topmed.git', type => 'git';

environment 'dev' => sub {
  set 'install_dir'=> $ENV{PWD},
  set 'build_dir'  => $ENV{PWD},
};

environment 'test' => sub {
  set 'install_dir' => '/var/tmp/monitor-topmed',
  set 'build_dir'   => '/tmp/monitor-topmed',
};

unless (environment eq 'dev') {
  umask 0002;
}

task 'build', sub {
  my $install_dir = get 'install_dir';
  my $build_dir   = get 'build_dir';
  my $local_lib   = "$build_dir/local";
  my $cpanm       = "$local_lib/bin/cpanm";
  my $carton      = "$local_lib/bin/carton";
  my $umask       = umask;

  unless (environment eq 'dev') {
    checkout 'master', path => $build_dir;
  }

  run "mkdir -p $local_lib/bin",
    unless => "test -d $local_lib/bin";

  run "curl -s -L https://cpanmin.us/ > $cpanm",
    unless => "test -x $cpanm", cwd => $build_dir;

  run "chmod 755 $cpanm",
    unless => "test -x $cpanm";

  run "$cpanm --self-contained --local-lib $local_lib Carton",
    unless => "test -e $local_lib/lib/perl5/Carton.pm";

  my $carton_cmd = case environment, {
    dev     => "$carton install",
    default => "$carton install --deployment",
  };

  run "umask $umask ; $carton_cmd",
    cwd     => $build_dir,
    env     => {
      PERL5LIB         => "$local_lib/lib/perl5",
      PERL_CARTON_PATH => "$local_lib",
    };
};

task 'install', sub {
  my $install_dir = get 'install_dir';
  my $build_dir   = get 'build_dir';

  run "rsync -a --exclude-from=.rsync-excludes $build_dir/ $install_dir/";
};

task 'schema', sub {
  run "$ENV{PWD}/scripts/schema.pl";
};
