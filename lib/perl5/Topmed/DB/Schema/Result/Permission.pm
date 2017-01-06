use utf8;
package Topmed::DB::Schema::Result::Permission;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::Permission

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

=head1 TABLE: C<permissions>

=cut

__PACKAGE__->table("permissions");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 centerid

  data_type: 'integer'
  is_nullable: 0

=head2 centername

  data_type: 'varchar'
  is_nullable: 0
  size: 16

=head2 runid

  data_type: 'integer'
  is_nullable: 0

=head2 dirname

  data_type: 'varchar'
  is_nullable: 0
  size: 64

=head2 operation

  data_type: 'char'
  is_nullable: 1
  size: 12

=cut

__PACKAGE__->add_columns(
  "id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "centerid",
  { data_type => "integer", is_nullable => 0 },
  "centername",
  { data_type => "varchar", is_nullable => 0, size => 16 },
  "runid",
  { data_type => "integer", is_nullable => 0 },
  "dirname",
  { data_type => "varchar", is_nullable => 0, size => 64 },
  "operation",
  { data_type => "char", is_nullable => 1, size => 12 },
);

=head1 PRIMARY KEY

=over 4

=item * L</id>

=back

=cut

__PACKAGE__->set_primary_key("id");


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-06 14:26:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:rnDCHM8iU80GIRaBcHelpw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
