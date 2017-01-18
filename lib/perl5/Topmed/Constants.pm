package Topmed::Constants;

use base qw(Exporter);
use Readonly;

use Topmed::Base;

our @EXPORT = (
  qw(
    $EMPTY
    $COMMA
    $UNDERSCORE
    $PERIOD
    $TRUE
    $FALSE
    $PIPE
    $DASH
    $SPACE
    $TAB
    $COLON
    $SLASH
    $MAX_DELAY
    $TIMEZONE
    )
);

our @EXPORT_OK = (
  qw(
    $EMPTY
    $COMMA
    $UNDERSCORE
    $PERIOD
    $TRUE
    $FALSE
    $PIPE
    $DASH
    $SPACE
    $TAB
    $COLON
    $SLASH
    $MAX_DELAY
    $TIMEZONE
    %GOOGLE_BUCKETS
    )
);

our %EXPORT_TAGS = (
  all    => \@EXPORT_OK,
  google => [
    qw(
      %GOOGLE_BUCKETS
      )
  ],
);

Readonly::Scalar our $EMPTY      => q{};
Readonly::Scalar our $COMMA      => q{,};
Readonly::Scalar our $UNDERSCORE => q{_};
Readonly::Scalar our $PERIOD     => q{.};
Readonly::Scalar our $TRUE       => q{1};
Readonly::Scalar our $FALSE      => q{0};
Readonly::Scalar our $PIPE       => q{|};
Readonly::Scalar our $DASH       => q{-};
Readonly::Scalar our $SPACE      => q{ };
Readonly::Scalar our $TAB        => qq{\t};
Readonly::Scalar our $COLON      => q{:};
Readonly::Scalar our $SLASH      => q{/};
Readonly::Scalar our $TIMEZONE   => q{America/Detroit};

Readonly::Hash our %GOOGLE_BUCKETS => (
  incoming => 'gs://topmed-incoming/',
  fastqs   => 'gs://topmed-fastqs/',
  crams    => 'gs://topmed-crams/',
  recabs   => 'gs://topmed-recabs/',
  logs     => 'gs://topmed-logs/',
);

1;
