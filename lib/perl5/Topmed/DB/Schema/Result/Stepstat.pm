use utf8;
package Topmed::DB::Schema::Result::Stepstat;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::Stepstat

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 COMPONENTS LOADED

=over 4

=item * L<DBIx::Class::InflateColumn::DateTime>

=back

=cut

__PACKAGE__->load_components("InflateColumn::DateTime");

=head1 TABLE: C<stepstats>

=cut

__PACKAGE__->table("stepstats");

=head1 ACCESSORS

=head2 yyyymmdd

  data_type: 'char'
  is_nullable: 0
  size: 10

=head2 count_verify

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 avetime_verify

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 count_bai

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 avetime_bai

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 count_qplot

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 avetime_qplot

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 count_cram

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 avetime_cram

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 count_expt

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 avetime_expt

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 ncbicount_expt

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 count_orig

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 avetime_orig

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 ncbicount_orig

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 count_b37

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 avetime_b37

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 ncbicount_b37

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 count_b38

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 avetime_b38

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 ncbicount_b38

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 bamcount

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 errcount

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 errorigcount

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 loadedorigbamcount

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 errckorigcount

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 errb37count

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 loadedb37bamcount

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 errckb37count

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 errb38count

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 loadedb38bamcount

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 errckb38count

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "yyyymmdd",
  { data_type => "char", is_nullable => 0, size => 10 },
  "count_verify",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "avetime_verify",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "count_bai",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "avetime_bai",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "count_qplot",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "avetime_qplot",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "count_cram",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "avetime_cram",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "count_expt",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "avetime_expt",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "ncbicount_expt",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "count_orig",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "avetime_orig",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "ncbicount_orig",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "count_b37",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "avetime_b37",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "ncbicount_b37",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "count_b38",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "avetime_b38",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "ncbicount_b38",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "bamcount",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "errcount",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "errorigcount",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "loadedorigbamcount",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "errckorigcount",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "errb37count",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "loadedb37bamcount",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "errckb37count",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "errb38count",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "loadedb38bamcount",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "errckb38count",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</yyyymmdd>

=back

=cut

__PACKAGE__->set_primary_key("yyyymmdd");


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-06 14:26:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:vjSAgZQNnw5M1/8p7TKziw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
