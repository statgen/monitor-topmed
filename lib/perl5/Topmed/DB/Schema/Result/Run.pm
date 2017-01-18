use utf8;
package Topmed::DB::Schema::Result::Run;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::Run

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

=head1 TABLE: C<runs>

=cut

__PACKAGE__->table("runs");

=head1 ACCESSORS

=head2 runid

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 centerid

  data_type: 'integer'
  is_nullable: 1

=head2 dirname

  data_type: 'varchar'
  is_nullable: 0
  size: 64

=head2 status

  data_type: 'varchar'
  is_nullable: 1
  size: 256

=head2 bamcount

  data_type: 'integer'
  is_nullable: 1

=head2 datayear

  data_type: 'integer'
  default_value: 2
  is_nullable: 1

=head2 offsite

  data_type: 'char'
  default_value: 'N'
  is_nullable: 1
  size: 1

=head2 xmlfound

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 arrived

  data_type: 'char'
  default_value: 'N'
  is_nullable: 1
  size: 1

=head2 dateinit

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datecomplete

  data_type: 'varchar'
  is_nullable: 1
  size: 10

=head2 comments

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "runid",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "centerid",
  { data_type => "integer", is_nullable => 1 },
  "dirname",
  { data_type => "varchar", is_nullable => 0, size => 64 },
  "status",
  { data_type => "varchar", is_nullable => 1, size => 256 },
  "bamcount",
  { data_type => "integer", is_nullable => 1 },
  "datayear",
  { data_type => "integer", default_value => 2, is_nullable => 1 },
  "offsite",
  { data_type => "char", default_value => "N", is_nullable => 1, size => 1 },
  "xmlfound",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "arrived",
  { data_type => "char", default_value => "N", is_nullable => 1, size => 1 },
  "dateinit",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datecomplete",
  { data_type => "varchar", is_nullable => 1, size => 10 },
  "comments",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</runid>

=back

=cut

__PACKAGE__->set_primary_key("runid");


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-06 14:26:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:JVxYT98VHBX50JOe5rkeTg


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->belongs_to(
  center =>
  'Topmed::DB::Schema::Result::Center',
  {'foreign.centerid' => 'self.centerid'}
);

__PACKAGE__->has_many(
  samples =>
  'Topmed::DB::Schema::Result::Bamfile',
  {'foreign.runid' => 'self.runid'}
);

1;
