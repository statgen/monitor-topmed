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

=head2 studyid

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

=head2 ncbierr

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 checksum

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 refname

  data_type: 'varchar'
  default_value: 'UNKNOWN'
  is_nullable: 1
  size: 96

=head2 expt_refname

  data_type: 'varchar'
  default_value: 'UNKNOWN'
  is_nullable: 1
  size: 96

=head2 expt_sampleid

  data_type: 'varchar'
  default_value: 'UNKNOWN'
  is_nullable: 1
  size: 96

=head2 nwdid_known

  data_type: 'char'
  default_value: 'N'
  is_nullable: 1
  size: 1

=head2 donot_remap

  data_type: 'char'
  default_value: (empty string)
  is_nullable: 1
  size: 4

=head2 state_arrive

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_md5ver

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_backup

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_cram

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_bai

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

=head2 state_ncbib38

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 time_ncbiexpt

  data_type: 'char'
  is_nullable: 1
  size: 19

=head2 time_ncbiorig

  data_type: 'char'
  is_nullable: 1
  size: 19

=head2 time_ncbib37

  data_type: 'char'
  is_nullable: 1
  size: 19

=head2 time_ncbib38

  data_type: 'char'
  is_nullable: 1
  size: 19

=head2 state_gce38push

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_gce38pull

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_gce38post

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 state_bcf

  data_type: 'integer'
  default_value: 0
  is_nullable: 1

=head2 datearrived

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datemd5ver

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 dateqplot

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datemapping

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datereport

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datebackup

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datecram

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datebai

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 datecp2ncbi

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidarrived

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidmd5ver

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidbackup

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidcram

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidbai

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidqplot

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidmapping

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidb37

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidb38

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidnwdid

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidncbiorig

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidncbib37

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 jobidncbib38

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 cramb37sent

  data_type: 'char'
  default_value: 'N'
  is_nullable: 1
  size: 1

=head2 cramb37checksum

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 cramb38checksum

  data_type: 'varchar'
  is_nullable: 0
  size: 96

=head2 jobidcp2ncbi

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 bam_delivered

  data_type: 'varchar'
  is_nullable: 1
  size: 12

=head2 cramorigsent

  data_type: 'char'
  default_value: 'N'
  is_nullable: 1
  size: 1

=head2 bamsize

  data_type: 'varchar'
  default_value: 0
  is_nullable: 1
  size: 16

=head2 datayear

  data_type: 'integer'
  default_value: 2
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
  "studyid",
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
  "ncbierr",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "checksum",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "refname",
  {
    data_type => "varchar",
    default_value => "UNKNOWN",
    is_nullable => 1,
    size => 96,
  },
  "expt_refname",
  {
    data_type => "varchar",
    default_value => "UNKNOWN",
    is_nullable => 1,
    size => 96,
  },
  "expt_sampleid",
  {
    data_type => "varchar",
    default_value => "UNKNOWN",
    is_nullable => 1,
    size => 96,
  },
  "nwdid_known",
  { data_type => "char", default_value => "N", is_nullable => 1, size => 1 },
  "donot_remap",
  { data_type => "char", default_value => "", is_nullable => 1, size => 4 },
  "state_arrive",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_md5ver",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_backup",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_cram",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_bai",
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
  "state_ncbib38",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "time_ncbiexpt",
  { data_type => "char", is_nullable => 1, size => 19 },
  "time_ncbiorig",
  { data_type => "char", is_nullable => 1, size => 19 },
  "time_ncbib37",
  { data_type => "char", is_nullable => 1, size => 19 },
  "time_ncbib38",
  { data_type => "char", is_nullable => 1, size => 19 },
  "state_gce38push",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_gce38pull",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_gce38post",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "state_bcf",
  { data_type => "integer", default_value => 0, is_nullable => 1 },
  "datearrived",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datemd5ver",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "dateqplot",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datemapping",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datereport",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datebackup",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datecram",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datebai",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "datecp2ncbi",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidarrived",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidmd5ver",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidbackup",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidcram",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidbai",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidqplot",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidmapping",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidb37",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidb38",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidnwdid",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidncbiorig",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidncbib37",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "jobidncbib38",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "cramb37sent",
  { data_type => "char", default_value => "N", is_nullable => 1, size => 1 },
  "cramb37checksum",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "cramb38checksum",
  { data_type => "varchar", is_nullable => 0, size => 96 },
  "jobidcp2ncbi",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "bam_delivered",
  { data_type => "varchar", is_nullable => 1, size => 12 },
  "cramorigsent",
  { data_type => "char", default_value => "N", is_nullable => 1, size => 1 },
  "bamsize",
  { data_type => "varchar", default_value => 0, is_nullable => 1, size => 16 },
  "datayear",
  { data_type => "integer", default_value => 2, is_nullable => 1 },
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

=head1 RELATIONS

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


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2017-02-08 15:52:57
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:uauA/Idko0hdTLbFzUULgQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
use Class::Method::Modifiers;

use Topmed::Base qw(files);
use Topmed::Constants qw(:google);

__PACKAGE__->belongs_to(
  run => 'Topmed::DB::Schema::Result::Run',
  {'foreign.runid' => 'self.runid'}
);

__PACKAGE__->belongs_to(
  study => 'Topmed::DB::Schema::Result::Study',
  {'foreign.studyid' => 'self.studyid'}
);

__PACKAGE__->belongs_to(
  mapping => 'Topmed::DB::Schema::Result::Mapping',
  {'foreign.bam_id' => 'self.bamid'}
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
