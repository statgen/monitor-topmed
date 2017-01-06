use utf8;
package Topmed::DB::Schema::Result::SummaryPassingGenome;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::SummaryPassingGenome

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

=head1 TABLE: C<summary_passing_genomes>

=cut

__PACKAGE__->table("summary_passing_genomes");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 col_name

  data_type: 'varchar'
  is_nullable: 1
  size: 100

=head2 mean

  data_type: 'float'
  is_nullable: 1
  size: [10,2]

=head2 sd

  data_type: 'float'
  is_nullable: 1
  size: [10,2]

=head2 min

  data_type: 'float'
  is_nullable: 1
  size: [10,2]

=head2 25-pct

  accessor: '25_pct'
  data_type: 'float'
  is_nullable: 1
  size: [10,2]

=head2 50-pct

  accessor: '50_pct'
  data_type: 'float'
  is_nullable: 1
  size: [10,2]

=head2 75-pct

  accessor: '75_pct'
  data_type: 'float'
  is_nullable: 1
  size: [10,2]

=head2 max

  data_type: 'float'
  is_nullable: 1
  size: [12,2]

=head2 label

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 label_indent

  data_type: 'integer'
  default_value: 0
  is_nullable: 0

=head2 id_compact

  data_type: 'integer'
  is_nullable: 1

=head2 stat_compacts

  data_type: 'integer'
  is_nullable: 1

=head2 number_fixed

  data_type: 'integer'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "col_name",
  { data_type => "varchar", is_nullable => 1, size => 100 },
  "mean",
  { data_type => "float", is_nullable => 1, size => [10, 2] },
  "sd",
  { data_type => "float", is_nullable => 1, size => [10, 2] },
  "min",
  { data_type => "float", is_nullable => 1, size => [10, 2] },
  "25-pct",
  {
    accessor => "25_pct",
    data_type => "float",
    is_nullable => 1,
    size => [10, 2],
  },
  "50-pct",
  {
    accessor => "50_pct",
    data_type => "float",
    is_nullable => 1,
    size => [10, 2],
  },
  "75-pct",
  {
    accessor => "75_pct",
    data_type => "float",
    is_nullable => 1,
    size => [10, 2],
  },
  "max",
  { data_type => "float", is_nullable => 1, size => [12, 2] },
  "label",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "label_indent",
  { data_type => "integer", default_value => 0, is_nullable => 0 },
  "id_compact",
  { data_type => "integer", is_nullable => 1 },
  "stat_compacts",
  { data_type => "integer", is_nullable => 1 },
  "number_fixed",
  { data_type => "integer", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</id>

=back

=cut

__PACKAGE__->set_primary_key("id");


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-06 14:26:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:J2+Zxao7M7tU1K8t7p/36Q


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
