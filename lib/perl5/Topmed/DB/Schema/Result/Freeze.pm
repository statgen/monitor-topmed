use utf8;
package Topmed::DB::Schema::Result::Freeze;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::Freeze

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

=head1 TABLE: C<freezes>

=cut

__PACKAGE__->table("freezes");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 bamid

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 freeze

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 created_at

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 0

=head2 modified_at

  data_type: 'timestamp'
  datetime_undef_if_invalid: 1
  default_value: current_timestamp
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "bamid",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "freeze",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "created_at",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 0,
  },
  "modified_at",
  {
    data_type => "timestamp",
    datetime_undef_if_invalid => 1,
    default_value => \"current_timestamp",
    is_nullable => 0,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</id>

=back

=cut

__PACKAGE__->set_primary_key("id");

=head1 RELATIONS

=head2 bamid

Type: belongs_to

Related object: L<Topmed::DB::Schema::Result::Bamfile>

=cut

__PACKAGE__->belongs_to(
  "bamid",
  "Topmed::DB::Schema::Result::Bamfile",
  { bamid => "bamid" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-09 15:30:58
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:0Hir77gMskfj51TNzfnJjg


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
