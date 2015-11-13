#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name:	topmed_xml.pl
#
# Description:
#   Use this program to create the XML files needed to submit
#   data to NCBI. Acts as a decision filter to make sure
#   the XML data is looks correct for submission
#
# ChangeLog:
#   $Log: nhlbi_xml.pl,v $
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

use IO::File;
use XML::Writer;
use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    topdir => '/net/topmed/incoming/topmed',
    broker_name => "UM-SPH",
    lab_name => "Abecasis",
    master => 'Tom Blackwell',
    master_email => 'tblackw@umich.edu',
    exptalias => 'not_set_yet',
    inst_model => 'HiSeq X Ten',
    design_description => 'Not set',
    build => 37,
    run_processing => 'run_processing.txt',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help nolint realm=s verbose=i build=i inst_model=s design_description=s run_processing=s
    )) || die "$Script Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 1 || $opts{help}) {
    warn "$Script [options] center run\n" .
        "Create XML files necessary to submit data to NCBI'.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $center = shift @ARGV;
my $run = shift @ARGV;

my $dbh = DBConnect($opts{realm});
my $now = time();

my $sql = "SELECT * FROM $opts{centers_table} WHERE centername='$center'";
my $sth = DoSQL($sql);
my $rowsofdata = $sth->rows();
if (! $rowsofdata) { die "$Script Unknown center '$center'\n"; }
my $href = $sth->fetchrow_hashref;
my $centerdesc = $href->{centerdesc};
if (defined($href->{designdesc}) && $href->{designdesc} ne '') {
    $opts{design_description} = $href->{designdesc};
}

my $d = $opts{topdir} . "/$center/$run";
if (! chdir($d)) {
    die "$Script Unable to CD to '$d', likely unknown run '$run'\n";
}

#--------------------------------------------------------------
#   Figure out names of XML files to be created
#   If they already exist, make another since data CAN be
#   submitted to NCBI more than once for the same run
#--------------------------------------------------------------
my %files = (
    submit => "nhlbi.$center.$run.submit",
    expt =>   "nhlbi.$center.$run.expt",
    run =>    "nhlbi.$center.$run.run"
);
foreach my $f (keys %files) {
    if (! -f $files{$f}) { next; }              # File does not exist, use this
    foreach my $i (1 .. 10) {                   # Else try another name
        if (! -f "$files{$f}$i") { $files{$f} = "$files{$f}$i"; last; } 
    }
    if (-f $files{$f}) { die "$Script Unable to figure out XML name for '$files{$f}'\n"; }
}
foreach my $f (keys %files) { $files{$f} = $files{$f} . '.xml'; }   # Full name of files

$opts{exptalias} = "$center.$run.$now";         # This ties expt with submit XML

#--------------------------------------------------------------
#   Create XML files or die   $center and $run are referenced globally
#--------------------------------------------------------------
CreateSubmit($files{submit});
RunLINT($files{submit});

#   Create XML of the experiment. Get array of runids to send
my $runidref = CreateExpt($files{expt});
RunLINT($files{expt});

#   Create XML of BAM details
CreateRun($files{run}, $runidref, $opts{run_processing});
RunLINT($files{run});

exit;

#==================================================================
# Subroutine:
#   CreateRun($f, $runidsref, $processingfile)
#
#   Create XML file of run information
#==================================================================
sub CreateRun {
    my ($f, $runidsref, $processingfile) = @_;
    
    #   Read in fixed XML data for all files in this run
    open(IN, $processingfile) ||
        die "$Script Unable to read 'run_processing' file '$processingfile': $!\n";
    my @l = <IN>;
    close(IN);
    my $runlines = join('', @l);        # Insert this in every RUN block

    #   Create the run XML file
    open(OUT, '>' . $f) ||
        die "$Script Unable to create file '$f': $!\n";

    print OUT "<?xml\n  version = \"1.0\"\n  encoding = \"UTF-8\"?>\n";
    print OUT "<RUN_SET\n" .
        "  xmlns:xsi = \"http://www.w3.org/2001/XMLSchema-instance\"\n" .
        "  xsi:noNamespaceSchemaLocation = \"http://www.ncbi.nlm.nih.gov/viewvc/v1/" .
            "trunk/sra/doc/SRA_1-5/SRA.run.xsd?view=co\">\n";

    #   Now walk through all bams for this run
    my $count = 0;
    foreach my $runid (@{$runidsref}) {
        my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE runid=$runid";
        my $sth = DoSQL($sql);
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script Run '$run' [$runid] has no BAM files\n"; }
        $href = $sth->fetchrow_hashref;
        $count++;
        print OUT "<RUN\n" .
            "  alias = \"$href->{checksum}\"\n" .
            "  center_name = \"$centerdesc\"\n" .
            "  broker_name = \"$opts{broker_name}\"\n" .
            "  run_center = \"$centerdesc\">\n" .
            "<TITLE>Secondary mapping for Build $opts{build}</TITLE>\n" .
            $runlines .
            "<DATA_BLOCK>\n" .
            "  <FILES>\n" .
            "    <FILE\n" .
            "      filename = \"$href->{bamname}\"\n" .
            "      filetype = \"bam\"\n" .
            "      checksum_method = \"MD5\"\n" .
            "      checksum = \"$href->{checksum}\">\n" .
            "    </FILE>\n" .
            "  </FILES>\n" .
            "</DATA_BLOCK>\n" .
            "</RUN>\n";
    }
    print OUT "</RUN_SET>\n";
    close(OUT);
    if (! $count) {
        unlink($f);
        die "$Script No eligible BAMs were found for '$center' '$run'. File '$f' not created.\n";
    }
    print "Created '$f' for $count BAMs\n";
}

#==================================================================
# Subroutine:
#   $aref = CreateExpt($f)
#
#   Create XML file of experiment information
#
#   Returns:  Reference to array of runids to be sent
#==================================================================
sub CreateExpt {
    my ($f) = @_;

    #   Create the experiment XML file
    open(OUT, '>' . $f) ||
        die "Script Unable to create file '$f': $!\n";

    print OUT "<?xml\n  version = \"1.0\"\n  encoding = \"UTF-8\"?>\n";
    print OUT "<EXPERIMENT_SET\n" .
        "  xmlns:xsi = \"http://www.w3.org/2001/XMLSchema-instance\"\n" .
        "  xsi:noNamespaceSchemaLocation = \"http://www.ncbi.nlm.nih.gov/viewvc/v1/" .
            "trunk/sra/doc/SRA_1-5/SRA.experiment.xsd?view=co\">\n";

    my $sql = "SELECT * FROM $opts{runs_table} WHERE dirname='$run'";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script Unknown center '$center'\n"; }
    my $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};                 # Uniq id for this run

    #   Now walk through all bams for this run and make an EXPERIMENT section
    $sql = "SELECT * FROM $opts{bamfiles_table} WHERE runid=$runid";
    $sth = DoSQL($sql);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script Run '$run' [$runid] has no BAM files\n"; }
    my @runids = ();                    # Keep track of all runs we found
    for (my $i=1; $i<=$rowsofdata; $i++) {
        $href = $sth->fetchrow_hashref;
        my $x = Experiment($href);
        if (! $x) { next; }             # Nope, did not like this BAM
        print OUT $x;
        push @runids, $href->{runid};   # Save for later
    }

    print OUT "</EXPERIMENT_SET>\n";
    close(OUT);
    if (! @runids) {
        unlink($f);
        die "$Script No eligible BAMs were found for '$center' '$run'. File '$f' not created.\n";
    }
    print "Created '$f' for " . (scalar(@runids)+1) . " BAMs\n";
    return \@runids;                    # Return list of BAMs
}

#==================================================================
# Subroutine:
#   $xmlstring = Experiment($href)
#
#   Create XML EXPERIMENT clause for this BAM
#
#   Returns:  string of XML or null if there was an error
#==================================================================
sub Experiment {
    my ($href) = @_;

    #   Check if this BAM can be sent
    foreach my $col (qw(checksum phs nominal_length nominal_sdev base_coord library_name expt_sampleid)) {
        if (exists($href->{$col}) && $href->{$col}) { next; }
        print "  Skipping BAM '$href->{bamname}' which has no value for '$col'\n";
        return '';
    }

    if (exists($href->{bam_delivered}) && $href->{bam_delivered} && ($href->{bam_delivered} > 10)) {
        print "  Skipping '$href->{bamname}' which has already been sent\n";
        return '';
    }

    #   Build the EXPERIMENT XML for this BAM - should never see uninitialized field
    my $xml = "<EXPERIMENT\n" .
        "  alias = \"$href->{checksum}\"\n" .
        "  center_name = \"$centerdesc\"\n" .
        "  broker_name = \"$opts{broker_name}\">\n" .
        "<TITLE>30x whole genome DNA sequence for NHLBI TOPMed sample</TITLE>\n" .
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
#   CreateSubmit($f)
#
#   Create XML file of submission information
#==================================================================
sub CreateSubmit {
    my ($f) = @_;

    #   Create the submit XML file
    open(OUT, '>' . $f) ||
        die "$Script Unable to create file '$f': $!\n";

    print OUT "<?xml\n  version = \"1.0\"\n  encoding = \"UTF-8\"?>\n";
    print OUT "<SUBMISSION\n" .
        "  xmlns:xsi = \"http://www.w3.org/2001/XMLSchema-instance\"\n" .
        "  xsi:noNamespaceSchemaLocation = \"http://www.ncbi.nlm.nih.gov/viewvc/v1/" .
            "trunk/sra/doc/SRA_1-5/SRA.submission.xsd?view=co\"\n" .
        "  alias = \"$opts{exptalias}\"\n" .
        "  center_name = \"$centerdesc\"\n" .
        "  broker_name = \"$opts{broker_name}\"\n" .
        "  lab_name = \"$opts{lab_name}\"\n" .
        "  submission_comment = \"\">\n";

    print OUT "<TITLE>$center $run $now</TITLE>\n";

    print OUT "<CONTACTS>\n" .
        "<CONTACT\n" .
        "  name = \"$opts{master}\"\n" .
        "  inform_on_status = \"$opts{master_email}\"\n" .
        "  inform_on_error = \"$opts{master_email}\"/>\n" .
        "</CONTACTS>\n";
 
    print OUT "<ACTIONS>\n  <ACTION>\n" .
        "    <ADD\n" .
        "      source = \"$files{expt}\"\n" .
        "      schema = \"experiment\"/>\n" .
        "  </ACTION>\n  <ACTION>\n" .
        "  <ADD\n" .
        "      source = \"$files{run}\"\n" .
        "      schema = \"run\"/>\n" .
        "  </ACTION>\n";
    print OUT "  <ACTION>\n    <PROTECT/>\n  </ACTION>\n</ACTIONS>\n";

    print OUT "</SUBMISSION>\n";
    close(OUT);
    print "Created '$f'\n";
}

#==================================================================
# Subroutine:
#   RunLINT($f)
#
#   Run XMLLINT on an XML file. Dies on error
#==================================================================
sub RunLINT {
    my ($f) = @_;
    if ($opts{nolint}) { return; }

    print "  Checking XML syntax for $f ... ";
    my $cmd = "xmllint $f > /dev/null";     # Only show errors
    my $rc = system($cmd);
    if (! $rc) { print "OK\n"; return; }
    print "ERROR\n";
    unlink($f);
    die "$Script XML file '$f' was mal-formated and has been removed\n";
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

=item B<-build nn>

Specifies the build number for the BAMs to be sent
The default is B<37>.

=item B<-design_description string>

Specifies the DESIGN_DESCRIPTION value in the experiment file.
This defaults to B<Equivalent to Illumina TruSeq PCR-free DNA sample prep>.

=item B<-help>

Generates this output.

=item B<-nolint>

Specifies that xmllint should not be run on the XML files.

=item B<-inst_model>

Specifies the INSTRUMENT_MODEL value in the experiment file.
This defaults to B<HiSeq X Ten>.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<topmed>.

=item B<-run_processing string>

Specifies the name of a file containing XML lines to be included in the run file.
These lines describe the processes used to create the BAM file, as well as
the names of the programs and their versions.
This differs for each run and will over time change.
This file must exist in the run directory and defaults to B<run_processing.txt>.

=item B<-verbose N>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<center>

Specifies the name of the center, e.g. B<broad>, B<uw> etc.

=item B<run>

Specifies the name of the directory (run) for the center, e.g. B<2015jun22> etc.

=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

