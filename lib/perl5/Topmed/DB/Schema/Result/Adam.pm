use utf8;
package Topmed::DB::Schema::Result::Adam;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::Adam

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

=head1 TABLE: C<adam>

=cut

__PACKAGE__->table("adam");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 name_1

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 md5sum_1

  data_type: 'varchar'
  is_nullable: 1
  size: 32

=head2 upload_id_1

  data_type: 'integer'
  is_nullable: 1

=head2 track_date_1

  data_type: 'varchar'
  is_nullable: 1
  size: 8

=head2 status1

  data_type: 'varchar'
  is_nullable: 1
  size: 16

=head2 name_2

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 md5sum_2

  data_type: 'varchar'
  is_nullable: 1
  size: 32

=head2 upload_id_2

  data_type: 'integer'
  is_nullable: 1

=head2 track_date_2

  data_type: 'varchar'
  is_nullable: 1
  size: 8

=head2 status_2

  data_type: 'varchar'
  is_nullable: 1
  size: 16

=cut

__PACKAGE__->add_columns(
  "id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "name_1",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "md5sum_1",
  { data_type => "varchar", is_nullable => 1, size => 32 },
  "upload_id_1",
  { data_type => "integer", is_nullable => 1 },
  "track_date_1",
  { data_type => "varchar", is_nullable => 1, size => 8 },
  "status1",
  { data_type => "varchar", is_nullable => 1, size => 16 },
  "name_2",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "md5sum_2",
  { data_type => "varchar", is_nullable => 1, size => 32 },
  "upload_id_2",
  { data_type => "integer", is_nullable => 1 },
  "track_date_2",
  { data_type => "varchar", is_nullable => 1, size => 8 },
  "status_2",
  { data_type => "varchar", is_nullable => 1, size => 16 },
);

=head1 PRIMARY KEY

=over 4

=item * L</id>

=back

=cut

__PACKAGE__->set_primary_key("id");


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-06 14:26:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:rUcYiSDSRB3V5vArVskZ8g


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
