use utf8;
package Topmed::DB::Schema::Result::NhlbiQcMetric;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::NhlbiQcMetric - VIEW

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
__PACKAGE__->table_class("DBIx::Class::ResultSource::View");

=head1 TABLE: C<nhlbi_qc_metrics>

=cut

__PACKAGE__->table("nhlbi_qc_metrics");
__PACKAGE__->result_source_instance->view_definition("select `b`.`expt_sampleid` AS `sample_id`,`b`.`piname` AS `pi_name`,`b`.`studyname` AS `study`,`c`.`centername` AS `center`,'2016-04-15' AS `seq_date`,from_unixtime(`b`.`dateinit`) AS `bam_date`,from_unixtime(`b`.`dateqplot`) AS `qc_date`,`q`.`pct_freemix` AS `pct_freemix`,`q`.`n_reads_m` AS `n_reads_m`,`q`.`pct_mapped` AS `pct_mapped`,`q`.`pct_mq0` AS `pct_mq0`,`q`.`pct_paired` AS `pct_paired`,`q`.`pct_prop_paired` AS `pct_prop_paired`,`q`.`mapped_gb` AS `mapped_gb`,`q`.`q20_gb` AS `q20_gb`,`q`.`pct_q20_base` AS `pct_q20_base`,`q`.`mean_depth` AS `mean_depth`,`q`.`pct_genome_cov` AS `pct_genome_cov`,`q`.`isize_mode` AS `isize_mode`,`q`.`isize_median` AS `isize_median`,`q`.`pct_dups` AS `pct_dups`,`q`.`pct_genome_dp5` AS `pct_genome_dp5`,`q`.`pct_genome_dp10` AS `pct_genome_dp10`,`q`.`pct_genome_dp20` AS `pct_genome_dp20`,`q`.`pct_genome_dp30` AS `pct_genome_dp30`,`q`.`vmr_depth` AS `vmr_depth`,`q`.`depth_q10` AS `depth_q10`,`q`.`depth_q20` AS `depth_q20`,`q`.`depth_q30` AS `depth_q30`,`q`.`raw_base_gb` AS `raw_base_gb`,`q`.`pct_overlap_reads` AS `pct_overlap_reads`,`q`.`pct_overlap_bases` AS `pct_overlap_bases`,`q`.`isize_iqr` AS `isize_iqr`,`q`.`isize_stdev` AS `isize_stdev`,`q`.`gc_depth_0_1` AS `gc_depth_0_1`,`q`.`gc_depth_1_5` AS `gc_depth_1_5`,`q`.`gc_depth_5_25` AS `gc_depth_5_25`,`q`.`gc_depth_25_75` AS `gc_depth_25_75`,`q`.`gc_depth_75_95` AS `gc_depth_75_95`,`q`.`gc_depth_95_99` AS `gc_depth_95_99`,`q`.`gc_depth_99_100` AS `gc_depth_99_100`,`q`.`library_size_m` AS `library_size_m`,0 AS `qc`,if(((`q`.`pct_freemix` < 3) and (`q`.`pct_genome_dp10` > 95) and (`q`.`mean_depth` > 30)),1,0) AS `qc_pass`,if((`q`.`mean_depth` < 30),1,0) AS `qc_flagged`,if(((`q`.`pct_freemix` > 3) or (`q`.`pct_genome_dp10` < 95)),1,0) AS `qc_fail`,from_unixtime(`b`.`datearrived`) AS `recieved`,`b`.`bamsize` AS `size`,`s_qplot`.`name` AS `status_qplot`,`s_b37`.`name` AS `status_remap_hg37`,`s_b38`.`name` AS `status_remap_hg38`,`s_backup`.`name` AS `status_backup`,`s_ncbi37`.`name` AS `status_ncbi_hg37`,`s_ncbi38`.`name` AS `status_ncbi_hg38`,from_unixtime(`b`.`datemapping`) AS `mapped_b37`,`b`.`datemapping_b38` AS `mapped_b38` from (((((((((`nhlbi`.`bamfiles` `b` join `nhlbi`.`runs` `r` on((`b`.`runid` = `r`.`runid`))) join `nhlbi`.`centers` `c` on((`r`.`centerid` = `c`.`centerid`))) join `nhlbi`.`states` `s_qplot` on((`b`.`state_qplot` = `s_qplot`.`id`))) left join `nhlbi`.`qc_results` `q` on((`q`.`bam_id` = `b`.`bamid`))) join `nhlbi`.`states` `s_b37` on((`b`.`state_b37` = `s_b37`.`id`))) join `nhlbi`.`states` `s_b38` on((`b`.`state_b38` = `s_b38`.`id`))) join `nhlbi`.`states` `s_backup` on((`b`.`state_backup` = `s_backup`.`id`))) join `nhlbi`.`states` `s_ncbi37` on((`b`.`state_ncbib37` = `s_ncbi37`.`id`))) join `nhlbi`.`states` `s_ncbi38` on((`b`.`state_ncbib38` = `s_ncbi38`.`id`)))");

=head1 ACCESSORS

=head2 sample_id

  data_type: 'varchar'
  is_nullable: 1
  size: 24

=head2 pi_name

  data_type: 'varchar'
  is_nullable: 1
  size: 96

=head2 study

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 center

  data_type: 'varchar'
  is_nullable: 0
  size: 16

=head2 seq_date

  data_type: 'varchar'
  default_value: (empty string)
  is_nullable: 0
  size: 10

=head2 bam_date

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=head2 qc_date

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

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

=head2 qc

  data_type: 'integer'
  default_value: 0
  is_nullable: 0

=head2 qc_pass

  data_type: 'integer'
  default_value: 0
  is_nullable: 0

=head2 qc_flagged

  data_type: 'integer'
  default_value: 0
  is_nullable: 0

=head2 qc_fail

  data_type: 'integer'
  default_value: 0
  is_nullable: 0

=head2 recieved

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=head2 size

  data_type: 'varchar'
  default_value: 0
  is_nullable: 1
  size: 16

=head2 status_qplot

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 status_remap_hg37

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 status_remap_hg38

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 status_backup

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 status_ncbi_hg37

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 status_ncbi_hg38

  data_type: 'varchar'
  is_nullable: 0
  size: 45

=head2 mapped_b37

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=head2 mapped_b38

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "sample_id",
  { data_type => "varchar", is_nullable => 1, size => 24 },
  "pi_name",
  { data_type => "varchar", is_nullable => 1, size => 96 },
  "study",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "center",
  { data_type => "varchar", is_nullable => 0, size => 16 },
  "seq_date",
  { data_type => "varchar", default_value => "", is_nullable => 0, size => 10 },
  "bam_date",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
  "qc_date",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
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
  "qc",
  { data_type => "integer", default_value => 0, is_nullable => 0 },
  "qc_pass",
  { data_type => "integer", default_value => 0, is_nullable => 0 },
  "qc_flagged",
  { data_type => "integer", default_value => 0, is_nullable => 0 },
  "qc_fail",
  { data_type => "integer", default_value => 0, is_nullable => 0 },
  "recieved",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
  "size",
  { data_type => "varchar", default_value => 0, is_nullable => 1, size => 16 },
  "status_qplot",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "status_remap_hg37",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "status_remap_hg38",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "status_backup",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "status_ncbi_hg37",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "status_ncbi_hg38",
  { data_type => "varchar", is_nullable => 0, size => 45 },
  "mapped_b37",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
  "mapped_b38",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
);


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-03-30 09:45:26
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:QDhP8MhWfQIuJ417SzFAww


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
