use utf8;
package Topmed::DB::Schema::Result::BadBamfile;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::BadBamfile

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

=head1 TABLE: C<bad_bamfiles>

=cut

__PACKAGE__->table("bad_bamfiles");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 bamid

  data_type: 'integer'
  is_nullable: 0

=head2 expt_sampleid

  data_type: 'varchar'
  is_nullable: 0
  size: 24

=head2 runid

  data_type: 'integer'
  is_nullable: 0

=head2 dateinit

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datayear

  data_type: 'integer'
  is_nullable: 0

=head2 build

  data_type: 'varchar'
  is_nullable: 1
  size: 4

=head2 bamname_orig

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 bamname

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 cramname

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 bamsize

  data_type: 'varchar'
  is_nullable: 0
  size: 16

=head2 checksum

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 cramchecksum

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 bamflagstat

  data_type: 'bigint'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 cramflagstat

  data_type: 'bigint'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 nominal_length

  data_type: 'integer'
  is_nullable: 1

=head2 nominal_sdev

  data_type: 'integer'
  is_nullable: 1

=head2 base_coord

  data_type: 'integer'
  is_nullable: 1

=head2 library_name

  data_type: 'varchar'
  is_nullable: 1
  size: 96

=head2 studyname

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 piname

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 phs

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 phs_consent_short_name

  data_type: 'varchar'
  is_nullable: 1
  size: 24

=head2 phs_sra_sample_id

  data_type: 'varchar'
  is_nullable: 1
  size: 24

=head2 phs_sra_data_details

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 nwdid_known

  data_type: 'char'
  is_nullable: 1
  size: 1

=head2 state_arrive

  data_type: 'integer'
  is_nullable: 1

=head2 state_md5ver

  data_type: 'integer'
  is_nullable: 1

=head2 state_cram

  data_type: 'integer'
  is_nullable: 1

=head2 state_bai

  data_type: 'integer'
  is_nullable: 1

=head2 state_qplot

  data_type: 'integer'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "bamid",
  { data_type => "integer", is_nullable => 0 },
  "expt_sampleid",
  { data_type => "varchar", is_nullable => 0, size => 24 },
  "runid",
  { data_type => "integer", is_nullable => 0 },
  "dateinit",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datayear",
  { data_type => "integer", is_nullable => 0 },
  "build",
  { data_type => "varchar", is_nullable => 1, size => 4 },
  "bamname_orig",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "bamname",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "cramname",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "bamsize",
  { data_type => "varchar", is_nullable => 0, size => 16 },
  "checksum",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "cramchecksum",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "bamflagstat",
  { data_type => "bigint", extra => { unsigned => 1 }, is_nullable => 1 },
  "cramflagstat",
  { data_type => "bigint", extra => { unsigned => 1 }, is_nullable => 1 },
  "nominal_length",
  { data_type => "integer", is_nullable => 1 },
  "nominal_sdev",
  { data_type => "integer", is_nullable => 1 },
  "base_coord",
  { data_type => "integer", is_nullable => 1 },
  "library_name",
  { data_type => "varchar", is_nullable => 1, size => 96 },
  "studyname",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "piname",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "phs",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "phs_consent_short_name",
  { data_type => "varchar", is_nullable => 1, size => 24 },
  "phs_sra_sample_id",
  { data_type => "varchar", is_nullable => 1, size => 24 },
  "phs_sra_data_details",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "nwdid_known",
  { data_type => "char", is_nullable => 1, size => 1 },
  "state_arrive",
  { data_type => "integer", is_nullable => 1 },
  "state_md5ver",
  { data_type => "integer", is_nullable => 1 },
  "state_cram",
  { data_type => "integer", is_nullable => 1 },
  "state_bai",
  { data_type => "integer", is_nullable => 1 },
  "state_qplot",
  { data_type => "integer", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</id>

=back

=cut

__PACKAGE__->set_primary_key("id");


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-25 10:12:38
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:zUA3IGlYaioomHxhJTU/PQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
