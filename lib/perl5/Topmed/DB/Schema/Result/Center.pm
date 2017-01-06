use utf8;
package Topmed::DB::Schema::Result::Center;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::Center

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

=head1 TABLE: C<centers>

=cut

__PACKAGE__->table("centers");

=head1 ACCESSORS

=head2 centerid

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 centername

  data_type: 'varchar'
  is_nullable: 0
  size: 16

=head2 datamethod

  data_type: 'varchar'
  default_value: 'pull'
  is_nullable: 0
  size: 8

=head2 centerdesc

  data_type: 'varchar'
  is_nullable: 1
  size: 96

=head2 designdesc

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "centerid",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "centername",
  { data_type => "varchar", is_nullable => 0, size => 16 },
  "datamethod",
  {
    data_type => "varchar",
    default_value => "pull",
    is_nullable => 0,
    size => 8,
  },
  "centerdesc",
  { data_type => "varchar", is_nullable => 1, size => 96 },
  "designdesc",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</centerid>

=back

=cut

__PACKAGE__->set_primary_key("centerid");


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-06 14:26:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:bxh8swF4cy9kWAUiQYb4BQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
