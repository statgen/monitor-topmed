# vim: set ft=perl
#
# TODO
#   - checkout source     - DONE
#   - get cpanminus       - DONE
#   - set exec on cpanm   - DONE
#   - install App::Carton - DONE
#   - run carton install  - DONE
#   - install changes     - DONE
#
use Rex -base;
use Rex::Commands::SCM;
use Rex::Commands::Sync;

my $install_dir = '/var/tmp/monitor';
my $build_dir   = '/tmp/monitor';
my $local_lib   = "$build_dir/local";
my $cpanm       = "$local_lib/bin/cpanm";
my $carton      = "$local_lib/bin/carton";

set repository => "master", url => 'git@github.com:statgen/monitor.git', type => 'git';

task 'build', sub {
  checkout 'master', path => $build_dir;

  run "mkdir -p $local_lib/bin",                 unless => "test -d $local_lib/bin";
  run "curl -s -L https://cpanmin.us/ > $cpanm", unless => "test -x $cpanm", cwd => $build_dir;
  run "chmod 755 $cpanm",                        unless => "test -x $cpanm";
  run "$cpanm --self-contained --local-lib $local_lib Carton";

  run "$carton install --deployment",
    cwd     => $build_dir,
    env     => {
      PERL5LIB         => "$local_lib/lib/perl5",
      PERL_CARTON_PATH => "$local_lib",
    };
};

task 'install', sub {
  run "rsync -a $build_dir/ $install_dir/";
};
