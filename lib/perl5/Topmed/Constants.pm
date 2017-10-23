package Topmed::Constants;

use base qw(Exporter);
use Readonly;

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
    $EQUAL
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
    $EQUAL
    $MAX_DELAY
    $TIMEZONE
    %GOOGLE_BUCKETS
    $GOOGLE_PIPELINE_CMD
    $GOOGLE_OPERATION_CMD
    $NOTSET
    $REQUESTED
    $SUBMITTED
    $STARTED
    $DELIVERED
    $COMPLETED
    $CANCELLED
    $FAILEDCHECKSUM
    $FAILED
    )
);

our %EXPORT_TAGS = (
  all    => \@EXPORT_OK,
  google => [
    qw(
      %GOOGLE_BUCKETS
      $GOOGLE_PIPELINE_CMD
      $GOOGLE_OPERATION_CMD
      )
  ],
  states => [
    qw(
      $NOTSET
      $REQUESTED
      $SUBMITTED
      $STARTED
      $DELIVERED
      $COMPLETED
      $CANCELLED
      $FAILEDCHECKSUM
      $FAILED
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
Readonly::Scalar our $EQUAL      => q{=};
Readonly::Scalar our $TIMEZONE   => q{America/Detroit};

Readonly::Hash our %GOOGLE_BUCKETS => (
  incoming => 'gs://topmed-incoming/',
  fastqs   => 'gs://topmed-fastqs/',
  crams    => 'gs://topmed-crams/',
  recabs   => 'gs://topmed-recabs/',
  logs     => 'gs://topmed-logs/',
  bcf      => 'gs://topmed-bcf/',
  mapping  => 'gs://topmed-mapping',
);

Readonly::Scalar our $GOOGLE_PIPELINE_CMD  => q{gcloud alpha genomics pipelines};
Readonly::Scalar our $GOOGLE_OPERATION_CMD => q{gcloud alpha genomics operations};

Readonly::Scalar our $NOTSET         => 0;     # Not set
Readonly::Scalar our $REQUESTED      => 1;     # Task requested
Readonly::Scalar our $SUBMITTED      => 2;     # Task submitted to be run
Readonly::Scalar our $STARTED        => 3;     # Task started
Readonly::Scalar our $DELIVERED      => 19;    # Data delivered, but not confirmed
Readonly::Scalar our $COMPLETED      => 20;    # Task completed successfully
Readonly::Scalar our $CANCELLED      => 89;    # Task cancelled
Readonly::Scalar our $FAILEDCHECKSUM => 98;    # Task failed, because checksum at NCBI bad
Readonly::Scalar our $FAILED         => 99;    # Task failed

1;
