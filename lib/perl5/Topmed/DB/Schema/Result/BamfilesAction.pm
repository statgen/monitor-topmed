use utf8;
package Topmed::DB::Schema::Result::BamfilesAction;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::BamfilesAction

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

=head1 TABLE: C<bamfiles_actions>

=cut

__PACKAGE__->table("bamfiles_actions");

=head1 ACCESSORS

=head2 bam_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

	

=head2 action_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 state_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 build_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

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
  "bam_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "action_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "state_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "build_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
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

=head1 UNIQUE CONSTRAINTS

=head2 C<index4>

=over 4

=item * L</bam_id>

=item * L</action_id>

=item * L</build_id>

=back

=cut

__PACKAGE__->add_unique_constraint("index4", ["bam_id", "action_id", "build_id"]);

=head1 RELATIONS

=head2 action

Type: belongs_to

Related object: L<Topmed::DB::Schema::Result::Action>

=cut

__PACKAGE__->belongs_to(
  "action",
  "Topmed::DB::Schema::Result::Action",
  { id => "action_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 bam

Type: belongs_to

Related object: L<Topmed::DB::Schema::Result::Bamfile>

=cut

__PACKAGE__->belongs_to(
  "bam",
  "Topmed::DB::Schema::Result::Bamfile",
  { bamid => "bam_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 build

Type: belongs_to

Related object: L<Topmed::DB::Schema::Result::Build>

=cut

__PACKAGE__->belongs_to(
  "build",
  "Topmed::DB::Schema::Result::Build",
  { id => "build_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 state

Type: belongs_to

Related object: L<Topmed::DB::Schema::Result::State>

=cut

__PACKAGE__->belongs_to(
  "state",
  "Topmed::DB::Schema::Result::State",
  { id => "state_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-03-30 09:45:26
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:05vco4XkZyztiFM+0idVQw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->load_components("InflateColumn::CreatedAt");
1;
