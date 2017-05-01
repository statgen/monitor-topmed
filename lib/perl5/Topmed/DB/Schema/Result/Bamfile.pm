use utf8;
package Topmed::DB::Schema::Result::Bamfile;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Topmed::DB::Schema::Result::Bamfile

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

=head1 TABLE: C<bamfiles>

=cut

__PACKAGE__->table("bamfiles");

=head1 ACCESSORS

=head2 bamid

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 runid

  data_type: 'integer'
  is_nullable: 0

=head2 bamname

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 base_coord

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 library_name

  data_type: 'varchar'
  is_nullable: 1
  size: 96

=head2 nominal_sdev

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 nominal_length

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 cramname

  data_type: 'varchar'
  is_nullable: 1
  size: 96

=head2 cramchecksum

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 b37cramchecksum

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 b38cramchecksum

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

=head2 b37flagstat

  data_type: 'bigint'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 b38flagstat

  data_type: 'bigint'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 datemapping_b38

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=head2 studyname

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 piname

  data_type: 'varchar'
  is_nullable: 1
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

=head2 emsg

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 checksum

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 expt_sampleid

  data_type: 'varchar'
  is_nullable: 1
  size: 24

=head2 nwdid_known

  data_type: 'char'
  default_value: 'N'
  is_nullable: 1
  size: 1

=head2 poorquality

  data_type: 'char'
  default_value: 'N'
  is_nullable: 1
  size: 1

=head2 state_arrive

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_verify

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_gce38backup

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_cram

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_qplot

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_b37

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_b38

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_ncbiexpt

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_ncbiorig

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_ncbib37

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_gce38push

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_gce38pull

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_bcf

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_gce38bcf_push

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_gce38bcf_pull

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_gce38bcf

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 gce38bcf_opid

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 state_38cp2gce

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 datearrived

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datemapping_b37

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=head2 bamsize

  data_type: 'varchar'
  default_value: 0
  is_nullable: 1
  size: 16

=head2 datayear

  data_type: 'integer'
  default_value: 3
  is_nullable: 1

=head2 build

  data_type: 'varchar'
  default_value: 37
  is_nullable: 1
  size: 4

=head2 dateinit

  data_type: 'varchar'
  is_nullable: 1
  size: 10

=head2 bamname_orig

  data_type: 'varchar'
  is_nullable: 1
  size: 96

=cut

__PACKAGE__->add_columns(
  "bamid",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "runid",
  { data_type => "integer", is_nullable => 0 },
  "bamname",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "base_coord",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "library_name",
  { data_type => "varchar", is_nullable => 1, size => 96 },
  "nominal_sdev",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "nominal_length",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "cramname",
  { data_type => "varchar", is_nullable => 1, size => 96 },
  "cramchecksum",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "b37cramchecksum",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "b38cramchecksum",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "bamflagstat",
  { data_type => "bigint", extra => { unsigned => 1 }, is_nullable => 1 },
  "cramflagstat",
  { data_type => "bigint", extra => { unsigned => 1 }, is_nullable => 1 },
  "b37flagstat",
  { data_type => "bigint", extra => { unsigned => 1 }, is_nullable => 1 },
  "b38flagstat",
  { data_type => "bigint", extra => { unsigned => 1 }, is_nullable => 1 },
  "datemapping_b38",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
  "studyname",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "piname",
  { data_type => "varchar", is_nullable => 1, size => 96 },
  "phs",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "phs_consent_short_name",
  { data_type => "varchar", is_nullable => 1, size => 24 },
  "phs_sra_sample_id",
  { data_type => "varchar", is_nullable => 1, size => 24 },
  "phs_sra_data_details",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "emsg",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "checksum",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "expt_sampleid",
  { data_type => "varchar", is_nullable => 1, size => 24 },
  "nwdid_known",
  { data_type => "char", default_value => "N", is_nullable => 1, size => 1 },
  "poorquality",
  { data_type => "char", default_value => "N", is_nullable => 1, size => 1 },
  "state_arrive",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_verify",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_gce38backup",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_cram",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_qplot",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_b37",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_b38",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_ncbiexpt",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_ncbiorig",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_ncbib37",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_gce38push",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_gce38pull",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_bcf",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_gce38bcf_push",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_gce38bcf_pull",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_gce38bcf",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "gce38bcf_opid",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "state_38cp2gce",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "datearrived",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datemapping_b37",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
  "bamsize",
  { data_type => "varchar", default_value => 0, is_nullable => 1, size => 16 },
  "datayear",
  { data_type => "integer", default_value => 3, is_nullable => 1 },
  "build",
  { data_type => "varchar", default_value => 37, is_nullable => 1, size => 4 },
  "dateinit",
  { data_type => "varchar", is_nullable => 1, size => 10 },
  "bamname_orig",
  { data_type => "varchar", is_nullable => 1, size => 96 },
);

=head1 PRIMARY KEY

=over 4

=item * L</bamid>

=back

=cut

__PACKAGE__->set_primary_key("bamid");

=head1 UNIQUE CONSTRAINTS

=head2 C<expt_sampleid>

=over 4

=item * L</expt_sampleid>

=back

=cut

__PACKAGE__->add_unique_constraint("expt_sampleid", ["expt_sampleid"]);

=head1 RELATIONS

=head2 bamfiles_actions

Type: has_many

Related object: L<Topmed::DB::Schema::Result::BamfilesAction>

=cut

__PACKAGE__->has_many(
  "bamfiles_actions",
  "Topmed::DB::Schema::Result::BamfilesAction",
  { "foreign.bam_id" => "self.bamid" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 freezes

Type: has_many

Related object: L<Topmed::DB::Schema::Result::Freeze>

=cut

__PACKAGE__->has_many(
  "freezes",
  "Topmed::DB::Schema::Result::Freeze",
  { "foreign.bamid" => "self.bamid" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 qc_results

Type: has_many

Related object: L<Topmed::DB::Schema::Result::QcResult>

=cut

__PACKAGE__->has_many(
  "qc_results",
  "Topmed::DB::Schema::Result::QcResult",
  { "foreign.bam_id" => "self.bamid" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-05-01 14:56:31
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:5aho/HvYdsqKqbYfUXMEcw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
use Class::Method::Modifiers;

use Topmed::Base qw(files);
use Topmed::Constants qw(:google);

__PACKAGE__->belongs_to(
  run => 'Topmed::DB::Schema::Result::Run',
  {'foreign.runid' => 'self.runid'}
);

around 'piname' => sub {
  my $orig = shift;
  my $self = shift;

  return $self->$orig // 'unknown';
};

sub center {
  return shift->run->center->centername;
}

sub nwdid {
  return shift->expt_sampleid;
}

sub runname {
  return shift->run->dirname;
}

sub host {
  my $self = shift;
  my $ptn  = YASF->new('/net/topmed/working/backups/incoming/topmed/{center}/{run}/{nwdid}.src.cram');
  my $path = $ptn % {center => $self->center, run => $self->run->dirname, nwdid => $self->expt_sampleid};
  my $file = Path::Class->file(abs_path($path));

  return ($file->components)[3];
}

sub b37_mapped_path {
  my $self = shift;
  my $path = YASF->new('/net/{host}/working/schelcj/results/{center}/{pi}/{nwdid}');

  return $path % {
    host   => $self->host,
    center => $self->center,
    pi     => $self->piname,
    nwdid  => $self->expt_sampleid,
  };
}

#   The path returned here is more convoluted than we want
#   because we started allocating b38 files on 9 and 10
#   and then after allocating a large number of samples
#   realized we needed more space.
#   Hence we must support two schemes to determine the path.
sub b38_mapped_path {
  my $self = shift;
  my @host_partialpath = (
    [ qw/topmed10 topmed9 topmed6  topmed7  topmed9  topmed10/ ],
    [ qw/working  working incoming incoming incoming incoming/ ]
  );
  # First determine the old path and if that directory exists, use it
  my $mod = $self->bamid % 2;
  my $path = YASF->new('/net/{host}/{part}/mapping/results/{center}/{pi}/b38/{nwdid}');
  my $outdir = $path % {
    host   => $host_partialpath[0][$mod],
    part => $host_partialpath[1][$mod],
    center => $self->center,
    pi     => $self->piname,
    nwdid  => $self->expt_sampleid,
 };
 if ( -e $outdir) { return $outdir; }       # If path exists, use that path

  # Path does not exist, allocate using the new scheme
  $mod = $self->bamid % 6;
  $path = YASF->new('/net/{host}/{part}/mapping/results/{center}/{pi}/b38/{nwdid}');
  $outdir = $path % {
    host   => $host_partialpath[0][$mod],
    part => $host_partialpath[1][$mod],
    center => $self->center,
    pi     => $self->piname,
    nwdid  => $self->expt_sampleid,
 };
 return $outdir;            # Caller must make directory if he wants it
}

sub gce_recab_uri {
  return $GOOGLE_BUCKETS{recabs} . shift->expt_sampleid;
}

sub incoming_path {
  my ($self) = @_;

  # Original BAM path:
  # /<prefix>/<host>/<project_incoming_dir>/<center>/<run_dir>/<filename>
  my $bam = File::Spec->join('/net/topmed/incoming', $self->center, $self->run->dirname, $self->bamname);

  # Backed up/squeezed path:
  # /<prefix>/<host>/<project_backup_dir>/<project_incoming_dir>/<center>/<run_dir>/<sample_id>.src.cram
  my $cram = File::Spec->join('/net/topmed/working/backups/incoming/topmed', $self->center, $self->run->dirname, $self->nwdid . '.src.cram');

  return $bam if -e $bam;
  return $cram if -e $cram;
  return;
}

sub b37_remapped_path {
  my ($self) = @_;
  return File::Spec->join($self->b37_mapped_path, 'bams', $self->nwdid . '.recal.cram');
}

sub b38_remapped_path {
  my ($self) = @_;
  return File::Spec->join($self->b38_mapped_path, $self->nwdid . '.recab.cram');
}

1;
