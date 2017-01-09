use utf8;
package Topmed::DB::Schema::Result::QcResult;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::QcResult

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

=head1 TABLE: C<qc_results>

=cut

__PACKAGE__->table("qc_results");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 bam_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 pct_freemix

  data_type: 'float'
  is_nullable: 1

=head2 n_reads_m

  data_type: 'float'
  is_nullable: 1

=head2 pct_mapped

  data_type: 'float'
  is_nullable: 1

=head2 pct_mq0

  data_type: 'float'
  is_nullable: 1

=head2 pct_paired

  data_type: 'float'
  is_nullable: 1

=head2 pct_prop_paired

  data_type: 'float'
  is_nullable: 1

=head2 mapped_gb

  data_type: 'float'
  is_nullable: 1

=head2 q20_gb

  data_type: 'float'
  is_nullable: 1

=head2 pct_q20_base

  data_type: 'float'
  is_nullable: 1

=head2 mean_depth

  data_type: 'float'
  is_nullable: 1

=head2 pct_genome_cov

  data_type: 'float'
  is_nullable: 1

=head2 isize_mode

  data_type: 'float'
  is_nullable: 1

=head2 isize_median

  data_type: 'float'
  is_nullable: 1

=head2 pct_dups

  data_type: 'float'
  is_nullable: 1

=head2 pct_genome_dp5

  data_type: 'float'
  is_nullable: 1

=head2 pct_genome_dp10

  data_type: 'float'
  is_nullable: 1

=head2 pct_genome_dp20

  data_type: 'float'
  is_nullable: 1

=head2 pct_genome_dp30

  data_type: 'float'
  is_nullable: 1

=head2 vmr_depth

  data_type: 'float'
  is_nullable: 1

=head2 depth_q10

  data_type: 'float'
  is_nullable: 1

=head2 depth_q20

  data_type: 'float'
  is_nullable: 1

=head2 depth_q30

  data_type: 'float'
  is_nullable: 1

=head2 raw_base_gb

  data_type: 'float'
  is_nullable: 1

=head2 pct_overlap_reads

  data_type: 'float'
  is_nullable: 1

=head2 pct_overlap_bases

  data_type: 'float'
  is_nullable: 1

=head2 isize_iqr

  data_type: 'float'
  is_nullable: 1

=head2 isize_stdev

  data_type: 'float'
  is_nullable: 1

=head2 gc_depth_0_1

  data_type: 'float'
  is_nullable: 1

=head2 gc_depth_1_5

  data_type: 'float'
  is_nullable: 1

=head2 gc_depth_5_25

  data_type: 'float'
  is_nullable: 1

=head2 gc_depth_25_75

  data_type: 'float'
  is_nullable: 1

=head2 gc_depth_75_95

  data_type: 'float'
  is_nullable: 1

=head2 gc_depth_95_99

  data_type: 'float'
  is_nullable: 1

=head2 gc_depth_99_100

  data_type: 'float'
  is_nullable: 1

=head2 library_size_m

  data_type: 'float'
  is_nullable: 1

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
  "bam_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "pct_freemix",
  { data_type => "float", is_nullable => 1 },
  "n_reads_m",
  { data_type => "float", is_nullable => 1 },
  "pct_mapped",
  { data_type => "float", is_nullable => 1 },
  "pct_mq0",
  { data_type => "float", is_nullable => 1 },
  "pct_paired",
  { data_type => "float", is_nullable => 1 },
  "pct_prop_paired",
  { data_type => "float", is_nullable => 1 },
  "mapped_gb",
  { data_type => "float", is_nullable => 1 },
  "q20_gb",
  { data_type => "float", is_nullable => 1 },
  "pct_q20_base",
  { data_type => "float", is_nullable => 1 },
  "mean_depth",
  { data_type => "float", is_nullable => 1 },
  "pct_genome_cov",
  { data_type => "float", is_nullable => 1 },
  "isize_mode",
  { data_type => "float", is_nullable => 1 },
  "isize_median",
  { data_type => "float", is_nullable => 1 },
  "pct_dups",
  { data_type => "float", is_nullable => 1 },
  "pct_genome_dp5",
  { data_type => "float", is_nullable => 1 },
  "pct_genome_dp10",
  { data_type => "float", is_nullable => 1 },
  "pct_genome_dp20",
  { data_type => "float", is_nullable => 1 },
  "pct_genome_dp30",
  { data_type => "float", is_nullable => 1 },
  "vmr_depth",
  { data_type => "float", is_nullable => 1 },
  "depth_q10",
  { data_type => "float", is_nullable => 1 },
  "depth_q20",
  { data_type => "float", is_nullable => 1 },
  "depth_q30",
  { data_type => "float", is_nullable => 1 },
  "raw_base_gb",
  { data_type => "float", is_nullable => 1 },
  "pct_overlap_reads",
  { data_type => "float", is_nullable => 1 },
  "pct_overlap_bases",
  { data_type => "float", is_nullable => 1 },
  "isize_iqr",
  { data_type => "float", is_nullable => 1 },
  "isize_stdev",
  { data_type => "float", is_nullable => 1 },
  "gc_depth_0_1",
  { data_type => "float", is_nullable => 1 },
  "gc_depth_1_5",
  { data_type => "float", is_nullable => 1 },
  "gc_depth_5_25",
  { data_type => "float", is_nullable => 1 },
  "gc_depth_25_75",
  { data_type => "float", is_nullable => 1 },
  "gc_depth_75_95",
  { data_type => "float", is_nullable => 1 },
  "gc_depth_95_99",
  { data_type => "float", is_nullable => 1 },
  "gc_depth_99_100",
  { data_type => "float", is_nullable => 1 },
  "library_size_m",
  { data_type => "float", is_nullable => 1 },
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


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-01-06 14:26:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:ppH7bFFfs3ss2ax9vtH8tA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
