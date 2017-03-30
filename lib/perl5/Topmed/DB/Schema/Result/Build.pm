use utf8;
package Topmed::DB::Schema::Result::Build;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::Build

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

=head1 TABLE: C<builds>

=cut

__PACKAGE__->table("builds");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=cut

__PACKAGE__->add_columns(
  "id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 45 },
);

=head1 PRIMARY KEY

=over 4

=item * L</id>

=back

=cut

__PACKAGE__->set_primary_key("id");

=head1 RELATIONS

=head2 bamfiles_actions

Type: has_many

Related object: L<Topmed::DB::Schema::Result::BamfilesAction>

=cut

__PACKAGE__->has_many(
  "bamfiles_actions",
  "Topmed::DB::Schema::Result::BamfilesAction",
  { "foreign.build_id" => "self.id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-03-30 09:45:26
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:HXgfz3aSe3FnjPVlLev8tQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
