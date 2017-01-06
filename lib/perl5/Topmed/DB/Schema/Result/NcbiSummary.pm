use utf8;
package Topmed::DB::Schema::Result::NcbiSummary;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::NcbiSummary

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

=head1 TABLE: C<ncbi_summary>

=cut

__PACKAGE__->table("ncbi_summary");

=head1 ACCESSORS

=head2 realm

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 upload_id

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 upload_date

  data_type: 'varchar'
  is_nullable: 1
  size: 20

=head2 file_name

  data_type: 'varchar'
  default_value: (empty string)
  is_nullable: 0
  size: 64

=head2 file_size

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 file_md5sum

  data_type: 'varchar'
  is_nullable: 1
  size: 32

=head2 upload_name

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 upload_size

  data_type: 'varchar'
  is_nullable: 1
  size: 24

=head2 upload_md5sum

  data_type: 'varchar'
  is_nullable: 1
  size: 32

=head2 file_status

  data_type: 'varchar'
  is_nullable: 1
  size: 20

=head2 file_type

  data_type: 'varchar'
  is_nullable: 1
  size: 20

=head2 load_date

  data_type: 'varchar'
  is_nullable: 1
  size: 24

=head2 file_error

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 submissions

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 loaded_runs

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 unloaded_runs

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 suppressed_runs

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 loaded_analyses

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 unloaded_analyses

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 suppressed_analyses

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=cut

__PACKAGE__->add_columns(
  "realm",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "upload_id",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "upload_date",
  { data_type => "varchar", is_nullable => 1, size => 20 },
  "file_name",
  { data_type => "varchar", default_value => "", is_nullable => 0, size => 64 },
  "file_size",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "file_md5sum",
  { data_type => "varchar", is_nullable => 1, size => 32 },
  "upload_name",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "upload_size",
  { data_type => "varchar", is_nullable => 1, size => 24 },
  "upload_md5sum",
  { data_type => "varchar", is_nullable => 1, size => 32 },
  "file_status",
  { data_type => "varchar", is_nullable => 1, size => 20 },
  "file_type",
  { data_type => "varchar", is_nullable => 1, size => 20 },
  "load_date",
  { data_type => "varchar", is_nullable => 1, size => 24 },
  "file_error",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "submissions",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "loaded_runs",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "unloaded_runs",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "suppressed_runs",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "loaded_analyses",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "unloaded_analyses",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "suppressed_analyses",
  { data_type => "varchar", is_nullable => 1, size => 64 },
);

=head1 PRIMARY KEY

=over 4

=item * L</file_name>

=back

=cut

__PACKAGE__->set_primary_key("file_name");


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-06 14:26:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:aymod6Ybh4sr6eziZc4QSQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
