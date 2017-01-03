#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
###################################################################
#
# Name:	topmed_xml.pl
#
# Description:
#   Use this program to create the XML files needed to submit
#   data to NCBI. Acts as a decision filter to make sure
#   the XML data is looks correct for submission
#
# Test sequence:
#   mkdir -p /tmp/test
#   ~/bin/topmed_xml.pl -xmlpr /tmp/test4/ -type expt 12891
#   ~/bin/topmed_xml.pl -xmlpr /tmp/test4/ -build 37 -type secondary 12891 /c/NWD559667.src.bam 12341234124
#   ~/bin/topmed_xml.pl -xmlpr /tmp/test4/ -build 37 -type remap 12891 /c/NWD559667.remap.bam  567856785678
#
# ChangeLog:
#   $Log: topmed_xml.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/perl5";
use My_DB;
use Getopt::Long;

use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my $xmlns_url = "xmlns:xsi = \"http://www.w3.org/2001/XMLSchema-instance\"\n" .
        "  xsi:noNamespaceSchemaLocation = \"http://www.ncbi.nlm.nih.gov/viewvc/v1/" .
            "trunk/sra/doc/SRA_1-5/SRA";
my $incomingdir = '/net/topmed/incoming';
my $topmed = '/net/topmed';
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    broker_name => "UM-SPH",
    lab_name => "Abecasis",
    master => 'Tom Blackwell',
    master_email => 'unattended@umich.edu',
    inst_model => 'HiSeq X Ten',
    design_description => 'Not set',
    build => 37,                    # 37, 38 etc
    run_processing => 'run_processing.txt',
    verbose => 0,
    xmlprefix => '',
    type => 'run',                  # Can be remap, secondary, or expt
    xmlns_run => "$xmlns_url.run.xsd?view=co\">\n",
    xmlns_experiment => "$xmlns_url.experiment.xsd?view=co\">\n",
    xmlns_submission => "$xmlns_url.submission.xsd?view=co\"\n",
    processingsectiondir => "$topmed/incoming/study.reference/seq.xml.extend",
    experiment_title => "30x whole genome DNA sequence for NHLBI TOPMed sample",
    schema =>  'http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA/SRA.%TYPE%.xsd?view=co',
);

Getopt::Long::GetOptions( \%opts,qw(
    help nolint realm=s verbose build=s inst_model=s design_description=s
    run_processing=s xmlprefix=s type=s bamchecksum=s master_email=s
    )) || die "$Script Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options -type expt] bamid\n" .
        "  or\n" .
        "$Script [options -type run ] bamid filename checksum\n" .
        "\n" .
        "Create XML files necessary to submit data to NCBI'.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my ($bamid, $filename, $checksum) = @ARGV;

if ($opts{type} ne 'secondary' && $opts{type} ne 'remap' && $opts{type} ne 'expt') {
    die "$Script Invalid option 'type' ($opts{type}). Must be 'secondary', 'remap' or 'expt'\n";
}
if ($opts{build} =~ /\D/) {
    die "$Script Invalid 'build' ($opts{build}). Must be a number like '37'\n";
}

my $dbh = DBConnect($opts{realm});

#   Figure out run and center for this bamid
my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE bamid=$bamid";
my $sth = DoSQL($sql);
my $rowsofdata = $sth->rows();
if (! $rowsofdata) { die "$Script Unknown bamid '$bamid'\n"; }
my $href = $sth->fetchrow_hashref;              # Save all info for this bamid

$sql = "SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$href->{runid}";
my $sth2 = DoSQL($sql);
$rowsofdata = $sth->rows();
if (! $rowsofdata) { die "$Script Unknown runid '$href->{runid}'  How'd that happen?\n"; }
my $href2 = $sth2->fetchrow_hashref;
my $run = $href2->{dirname};

$sql = "SELECT * FROM $opts{centers_table} WHERE centerid=$href2->{centerid}";
$sth2 = DoSQL($sql);
$rowsofdata = $sth->rows();
if (! $rowsofdata) { die "$Script Unknown centerid, how'd that happen?\n"; }
$href2 = $sth2->fetchrow_hashref;
my $center = $href2->{centername};
my $centerdesc = $href2->{centerdesc};
if (defined($href2->{designdesc}) && $href2->{designdesc} ne '') {
    $opts{design_description} = $href2->{designdesc};
}

#   CD some place nice and safe, just in case.  I don't think this is necessary
my $d = '/tmp';
if (! chdir($d)) {
    die "$Script Unable to CD to '$d', how can that happen?\n";
}
print "Creating XML for BAMID '$bamid' from $center/$run\n";

#--------------------------------------------------------------
#   Create XML files
#--------------------------------------------------------------
if ($opts{xmlprefix} && $opts{xmlprefix} !~ /[\/\/]$/) { $opts{xmlprefix} .= '.'; }

if ($opts{type} eq 'expt') {
    my $submitfile = $opts{xmlprefix} . "$href->{expt_sampleid}-$opts{type}.submit.xml";
    my $exptfile   = $opts{xmlprefix} . "$href->{expt_sampleid}.expt.xml";
    CreateSubmit($submitfile, "$href->{expt_sampleid}.expt.submit", 'expt', $exptfile);
    RunLINT($submitfile, 'submission');
    CreateExpt($exptfile, $href);
    RunLINT($exptfile, 'experiment');
    exit;
}

if ((! $filename) || (! $checksum)) { die "$Script filename or CHECKSUM not provided\n"; }

my $submitfile = $opts{xmlprefix} . "$href->{expt_sampleid}-$opts{type}.$opts{build}.submit.xml";
my $runfile = $opts{xmlprefix} . "$href->{expt_sampleid}-$opts{type}.$opts{build}.run.xml";
CreateSubmit($submitfile, "$href->{expt_sampleid}.run.$opts{type}.submit", 'run', $runfile);
RunLINT($submitfile, 'submission');
CreateRun($runfile, $href, $filename, $checksum);
RunLINT($runfile, 'run' );
exit;

#==================================================================
# Subroutine:
#   CreateRun($f, $href, $filename, $checksum)
#
#   Create XML file of run information
#   References global data
#==================================================================
sub CreateRun {
    my ($f, $href, $filename, $checksum) = @_;
    my ($processingfile, $title);
    $filename = basename($filename);            # No fully qualified path

    if ($opts{type} eq 'secondary') {           # Original bam
        $processingfile = "$opts{processingsectiondir}/$center.$href->{build}.run_processing.txt";
        $title = "Secondary mapping Build $opts{build} for $href->{expt_sampleid} [$href->{bamid}]";
        if ($filename =~ /(^NWD\d+)/) {         # BAM we send should be called NWDxxxxx.src.bam
            my $n = $1;
            if ($filename =~ /cram/) {          # Unless it is a CRAM from a bam
                $filename = $n . '.src.cram';
            }
            else { $filename = $n . '.src.bam'; }
        }
    }
    else {
        if (! defined($filename)) { die "$Script filename not provided for RUN XML\n"; }
        if (! defined($checksum)) { die "$Script Checksum not provided for file '$filename'\n"; }
        $processingfile =  "$opts{processingsectiondir}/umich.$opts{build}.run_processing.txt";
        $title = 'Undefined';
        if ($opts{type} eq 'remap') { $title = 'Remapped file'; }
        if ($opts{type} eq 'primary') { $title = 'Primary mapping'; }
        $title .= " Build $opts{build} for $href->{expt_sampleid} [$href->{bamid}]";
    }

    #   Read in fixed XML data for all files in this run
    open(IN, $processingfile) ||
        die "$Script Unable to read file '$processingfile': $!\n";
    my @l = <IN>;
    close(IN);
    my $runlines = join('', @l);    # Insert this in every RUN block

    #   Create the run XML file
    open(OUT, '>' . $f) ||
        die "$Script Unable to create file '$f': $!\n";
    print OUT "<?xml\n  version = \"1.0\"\n  encoding = \"UTF-8\"?>\n";
    print OUT "<RUN_SET\n" . $opts{xmlns_run};
    my $refname = "$href->{expt_sampleid}-expt";
    #   What a POS.  NYGC created their own experiments, so we need to know what
    #   they used for the name. Fortunately, this mess only belongs to Burchard
    if ($center eq 'nygc' &&
        $href->{piname} eq 'Burchard' &&
        $href->{datayear} eq '1') {
        $refname = $href->{expt_sampleid};
    }
    my $xml = GenRUNXML($checksum, $title, $refname, $runlines, $filename, $checksum);
    print OUT $xml . "</RUN_SET>\n";
    close(OUT);
    print "  Created Run XML $f\n";
}

#==================================================================
# Subroutine:
#   $xml = GenRUNXML($alias, $title, $refname, $runlines, $filename, $checksum)
#
#   Create XML fragment for <RUN> clause
#==================================================================
sub GenRUNXML {
    my ($alias, $title, $refname, $runlines, $filename, $checksum) = @_;

    my $ext = 'bam';
    if ($filename =~ /\.(\w+)$/) { $ext = $1; }  # Should be cram or bam
    my $provider = 'IRC Harmonized';
    if ($filename =~ /\.src\./) { $provider = 'Sequencing Center'; }
    my $xml .= "<RUN\n" .
        "  alias = \"$alias\"\n" .
        "  center_name = \"$centerdesc\"\n" .
        "  broker_name = \"$opts{broker_name}\"\n" .
        "  run_center = \"$centerdesc\">\n" .
        "<TITLE>$title</TITLE>\n" .
        "<EXPERIMENT_REF\n" .
        "  refname = \"$refname\"\n" .
        "  refcenter = \"$centerdesc\"/>\n" .
        $runlines .
        "<DATA_BLOCK>\n" .
        "  <FILES>\n" .
        "    <FILE\n" .
        "      filename = \"$filename\"\n" .
        "      filetype = \"$ext\"\n" .
        "      checksum_method = \"MD5\"\n" .
        "      checksum = \"$checksum\"></FILE>\n" .
        "  </FILES>\n" .
        "</DATA_BLOCK>\n" .
        "<RUN_ATTRIBUTES>\n" .
        "  <RUN_ATTRIBUTE>\n" .
        "    <TAG>assembly</TAG>\n" .
        "    <VALUE>GRCh$opts{build}</VALUE>\n" .
        "  </RUN_ATTRIBUTE>\n" .

        "  <RUN_ATTRIBUTE>\n" .                     #   Add another run_attribute per Adam Stine
        "    <TAG>Alignment Provider</TAG>\n" .
        "    <VALUE>$provider</VALUE>\n" .          # Sequencing Center for src or IRC Harmonized
        "  </RUN_ATTRIBUTE>\n" .

        "</RUN_ATTRIBUTES>\n" .
        "</RUN>\n";
    return $xml;
}

#==================================================================
# Subroutine:
#   $aref = CreateExpt($f, $href)
#
#   Create XML file of experiment information
#   References global data
#==================================================================
sub CreateExpt {
    my ($f, $href) = @_;

    #   Create the experiment XML file
    open(OUT, '>' . $f) ||
        die "Script Unable to create file '$f': $!\n";

    print OUT "<?xml\n  version = \"1.0\"\n  encoding = \"UTF-8\"?>\n";
    print OUT "<EXPERIMENT_SET\n" . $opts{xmlns_experiment};
    my $x = Experiment($href);          # Make the EXPERIMENT section
    print OUT $x;
    print OUT "</EXPERIMENT_SET>\n";
    close(OUT);
    print "  Created Experiment XML $f\n";
}

#==================================================================
# Subroutine:
#   $xmlstring = Experiment($href)
#
#   Create XML EXPERIMENT clause for this file
#
#   Returns:  string of XML or null if there was an error
#==================================================================
sub Experiment {
    my ($href) = @_;

    #   Build the EXPERIMENT XML for this BAM - should never see uninitialized field
    my $xml = "<EXPERIMENT\n" .
        "  alias = \"$href->{expt_sampleid}-expt\"\n" .
        "  center_name = \"$centerdesc\"\n" .
        "  broker_name = \"$opts{broker_name}\">\n" .
        "<TITLE>$opts{experiment_title} $href->{expt_sampleid}</TITLE>\n" .
        "<STUDY_REF\n" .
        "  accession = \"$href->{phs}\"/>\n" .
        "<DESIGN>\n" .
        "  <DESIGN_DESCRIPTION>$opts{design_description}</DESIGN_DESCRIPTION>\n" .
        "  <SAMPLE_DESCRIPTOR\n" .
        "    refname = \"$href->{expt_sampleid}\"\n" .
        "    refcenter = \"$href->{phs}\"/>\n";
    $xml .= 
        "  <LIBRARY_DESCRIPTOR>\n" .
        "    <LIBRARY_NAME>$href->{library_name}</LIBRARY_NAME>\n" .
        "    <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>\n" .
        "    <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>\n" .
        "    <LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>\n" .
        "    <LIBRARY_LAYOUT>\n";
    $xml .=
        "    <PAIRED\n" .
	    "      NOMINAL_LENGTH = \"$href->{nominal_length}\"\n" .
	    "      NOMINAL_SDEV = \"$href->{nominal_sdev}\"/>\n" .
	    "    </LIBRARY_LAYOUT>\n  </LIBRARY_DESCRIPTOR>\n";
	$xml .= 
        "  <SPOT_DESCRIPTOR>\n" .
        "    <SPOT_DECODE_SPEC>\n" .
        "      <READ_SPEC>\n" .
        "        <READ_INDEX>0</READ_INDEX>\n" .
        "        <READ_LABEL>forward</READ_LABEL>\n" .
        "        <READ_CLASS>Application Read</READ_CLASS>\n" .
        "        <READ_TYPE>Forward</READ_TYPE>\n" .
        "        <BASE_COORD>1</BASE_COORD>\n" .
        "      </READ_SPEC>\n" .
        "      <READ_SPEC>\n" .
        "        <READ_INDEX>1</READ_INDEX>\n" .
        "        <READ_LABEL>reverse</READ_LABEL>\n" .
        "        <READ_CLASS>Application Read</READ_CLASS>\n" .
        "        <READ_TYPE>Reverse</READ_TYPE>\n" .
        "        <BASE_COORD>$href->{base_coord}</BASE_COORD>\n" .
        "      </READ_SPEC>\n" .
        "    </SPOT_DECODE_SPEC>\n" .
        "  </SPOT_DESCRIPTOR>\n";
    $xml .=
        "</DESIGN>\n" .
        "<PLATFORM>\n" .
        "  <ILLUMINA>\n" .
        "    <INSTRUMENT_MODEL>$opts{inst_model}</INSTRUMENT_MODEL>\n" .
        "  </ILLUMINA>\n" .
        "</PLATFORM>\n" .
        "</EXPERIMENT>\n";
    return $xml;
}

#==================================================================
# Subroutine:
#   CreateSubmit($f, $alias, $t, $f2)
#
#   Create XML file of submission information
#   References global data
#==================================================================
sub CreateSubmit {
    my ($f, $alias, $t, $f2) = @_;
    my $action = 'MODIFY';
    if ($t eq 'expt') { $action = 'ADD'; $t = 'experiment'; }
    $f2 = basename($f2);

    #   Create the submit XML file
    open(OUT, '>' . $f) ||
        die "$Script Unable to create file '$f': $!\n";

    print OUT "<?xml\n  version = \"1.0\"\n  encoding = \"UTF-8\"?>\n";
    print OUT "<SUBMISSION\n" . $opts{xmlns_submission} .
        "  alias = \"$alias\"\n" .
        "  center_name = \"$centerdesc\"\n" .
        "  broker_name = \"$opts{broker_name}\"\n" .
        "  lab_name = \"$opts{lab_name}\"\n" .
        "  submission_comment = \"\">\n";

    print OUT "<TITLE>$center $run $t</TITLE>\n";

    print OUT "<CONTACTS>\n" .
        "<CONTACT\n" .
        "  name = \"$opts{master}\"\n";
    if ($opts{master_email}) {
        print OUT "  inform_on_status = \"$opts{master_email}\"\n" .
        "  inform_on_error = \"$opts{master_email}\"";
    }
    print OUT "/>\n" .
        "</CONTACTS>\n";
 
    print OUT "<ACTIONS>\n" .
        "  <ACTION>\n" .
        "    <$action\n" .
        "      source = \"$f2\"\n" .
        "      schema = \"$t\"/>\n" .
        "  </ACTION>\n" .
        "  <ACTION>\n" .
        "    <PROTECT/>\n" .
        "  </ACTION>\n" .
        "</ACTIONS>\n" .
        "</SUBMISSION>\n";
    close(OUT);
    print "  Created Submit XML $f\n";
}

#==================================================================
# Subroutine:
#   RunLINT($f, $schema)
#
#   Run XMLLINT on an XML file. Dies on error
#==================================================================
sub RunLINT {
    my ($f, $schema) = @_;
    if ($opts{nolint}) { return; }

    print "  Checking XML syntax for $f ... ";
    my $s = $opts{schema};
    $s =~ s/%TYPE%/$schema/;
    my $cmd = "xmllint --schema $s $f > /dev/null";     # Only show errors
    my $rc = system($cmd);
    if (! $rc) { print "OK\n"; return; }
    print "ERROR\n";
    my $emsg = "$Script XML file '$f' was mal-formed and renamed to '$f.failed'\n";
    rename($f,"$f.failed");
    die $emsg . "\n";
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_xml.pl - Create XML files so we can sent data to NCBI

=head1 SYNOPSIS

  topmed_xml.pl  broad  2015jun22

=head1 DESCRIPTION

This program generates the three XML files necessary so we can send
data to NCBI.  It verifies the fields in the database.
This is expected to be used in a shell script that will drive the
process of preparing the files to be sent and then copy them there.

=head1 OPTIONS

=over 4

=item B<-build nn|orig>

Specifies the build number for the BAMs to be sent
The default is B<37>.
This field controls XML included which describes the programs used
to create the BAM file.
These fragments of XML live in B<$opts{processingsectiondir}>.

=item B<-design_description string>

Specifies the DESIGN_DESCRIPTION value in the experiment file.
This defaults to B<Equivalent to Illumina TruSeq PCR-free DNA sample prep>.

=item B<-help>

Generates this output.

=item B<-inst_model>

Specifies the INSTRUMENT_MODEL value in the experiment file.
This defaults to B<HiSeq X Ten>.

=item B<-master_email emailaddress>

Specifies the error and notification Email address to be used when submitting
XML. The default is B<unattended@umich.edu>.
Specifying a null string will remove the Email address from the XML.

=item B<-nolint>

Specifies that xmllint should not be run on the XML files.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<topmed>.

=item B<-run_processing string>

Specifies the name of a file containing XML lines to be included in the run file.
These lines describe the processes used to create the BAM file, as well as
the names of the programs and their versions.
This differs for each run and will over time change.
This file must exist in the run directory and defaults to B<run_processing.txt>.

=item B<-type run|expt>

Specifies the type of XML to be sent, XML for a run or an experiment.

=item B<-verbose N>

Provided for developers to see additional information.

=item B<-xmlprefix string>

Specifies a prefix for the three XML files to be created.
This allows the caller to make these unique.
The XML files will normally begin with the NWDID.


=back

=head1 PARAMETERS

=over 4

=item B<bamid>

Specifies the bam database identifier, e.g. B<5467>.

=item B<filename>

Specifies the filename of the BAM or CRAM file.
Be sure you have specified the correct B<build> for this file.
This is not verified.

=item B<checksum>

Specifies the checksum of the bam file.
This is not verified.

=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

