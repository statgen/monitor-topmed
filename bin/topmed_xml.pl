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

use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my $xmlns_url = "xmlns:xsi = \"http://www.w3.org/2001/XMLSchema-instance\"\n" .
        "  xsi:noNamespaceSchemaLocation = \"http://www.ncbi.nlm.nih.gov/viewvc/v1/" .
            "trunk/sra/doc/SRA_1-5/SRA";
my $incomingdir = '/net/topmed/incoming';
my $topmed = '/net/topmed';
my $topmed2 = '/net/topmed2';
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    topdir      => "$topmed/incoming/topmed",
    topdir2     => "$topmed2/incoming/topmed",
    backupdir   => "$topmed/working/backups/incoming/topmed",
    backupdir2  => "$topmed2/working/backups/incoming/topmed",
    resultsdir  => "$topmed/working/schelcj/results",
    resultsdir2 => "$topmed2/working/schelcj/results",
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
    xmlns_run => "$xmlns_url.run.xsd?view=co\">\n",
    xmlns_experiment => "$xmlns_url.experiment.xsd?view=co\">\n",
    xmlns_submission => "$xmlns_url.submission.xsd?view=co\">\n",
    processingsectiondir => "$topmed/incoming/study.reference/send2ncbi",
    experiment_title => "30x whole genome DNA sequence for NHLBI TOPMed sample",
);

Getopt::Long::GetOptions( \%opts,qw(
    help nolint realm=s verbose build=i inst_model=s design_description=s run_processing=s
    bamlist=s xmlprefix=s sendlist=s 
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
#   Prepare input to create XML
#--------------------------------------------------------------
#   Figure out names of XML files to be created
#   If they already exist, they will be overwritten
if (! $opts{xmlprefix}) { $opts{xmlprefix} = "nhlbi.$center.$run"; }
my %files = (
    submit => $opts{xmlprefix} . '.submit.xml',
    expt =>   $opts{xmlprefix} . '.expt.xml',
    run =>    $opts{xmlprefix} . '.run.xml'
);
$opts{exptalias} = "$center.$run.$now";         # This ties expt with submit XML

#   Get list of all BAMIDs
my @bamids = ();
if ($opts{bamlist}) {
    open(IN,$opts{bamlist}) ||
        die "$Script Unable to open list of BAM files '$opts{bamlist}': $!\n";
    @bamids = <IN>;
    close(IN);
}
else {                          # List not provided, get all bams for this center/RUN
    my $sql = "SELECT * FROM $opts{runs_table} WHERE dirname='$run'";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script Unknown run '$run'\n"; }
    my $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};                 # Uniq id for this run
    $sql = "SELECT * FROM $opts{bamfiles_table} WHERE runid=$runid";
    $sth = DoSQL($sql);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script Run '$run' [$runid] has no BAM files\n"; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        $href = $sth->fetchrow_hashref;
        push @bamids,$href->{bamid};
    }
    print "Using all BAMs for $center/$run. This is almost certainly not what you want\n";
    print "Do you really want to continue? (Y/N): ";
    $_ = <STDIN>;
    if (! /^y/i) { die "Nothing done\n"; }
}
if (! @bamids) { die "$Script No BAMIDs were provided\n"; }
print "Creating XML for a total of " . scalar(@bamids) . " BAMS in $center/$run\n";

#--------------------------------------------------------------
#   Create XML files or die
#--------------------------------------------------------------
#   Initialize the file of files to send
if ($opts{sendlist}) {
    open(APP, '>' . $opts{sendlist}) ||
        die "$Script Unable to create to file '$opts{sendlist}': $!\n";
    print APP "# XML This is a list of files to be sent to NCBI\n";
    close(APP);
}

CreateSubmit($files{submit}, $files{run}, $files{expt});
RunLINT($files{submit});

#   Create experiment XML from the list of bamids
CreateExpt($files{expt}, \@bamids);
RunLINT($files{expt});
my $firstbamid = $bamids[0];            # Use this as alias for failure if XML files
chomp($firstbamid);

#   Create XML of BAM details
#   Append fully qualified path of files to send in $opts{sendlist} on one line
CreateRun($files{run}, \@bamids, "$opts{processingsectiondir}/$center.run_processing.txt",
    "$opts{processingsectiondir}/um$opts{build}.run_processing.txt", $opts{sendlist});
RunLINT($files{run});

#   Last thing to be sent are the XML files themselves
if ($opts{sendlist}) {
    open(APP, '>>' . $opts{sendlist}) ||
        die "$Script Unable to append to file '$opts{sendlist}': $!\n";
    print APP "$firstbamid XML $files{submit} $files{expt} $files{run}\n";
    close(APP);
    print "Created list of files to send in '$opts{sendlist}'\n";
}

exit;

#==================================================================
# Subroutine:
#   CreateRun($f, $bamidlistref, $centerprocessingfile, $remapprocessingfile, $path_files2send)
#
#   Create XML file of run information
#   For each file in the set of bams, append this information to $path_files2send
#     fully qualified path to file  + database column to be set
#   The contents of these lines become part of the command line to topmed_ncbi.sh
#==================================================================
sub CreateRun {
    my ($f, $bamidlistref, $centerprocessingfile, $buildprocessingfile, $path_files2send) = @_;

    #   Read in fixed XML data for all files in this run
    open(IN, $centerprocessingfile) ||
        die "$Script Unable to read file '$centerprocessingfile': $!\n";
    my @l = <IN>;
    close(IN);
    my $centerrunlines = join('', @l);  # Insert this in every RUN block

    #   Read in fixed XML data for all remapped files in this run
    open(IN, $buildprocessingfile) ||
        die "$Script Unable to read file '$buildprocessingfile': $!\n";
    @l = <IN>;
    close(IN);
    my $buildrunlines = join('', @l);   # Insert this in every remapped RUN block

    #   Create the run XML file
    open(OUT, '>' . $f) ||
        die "$Script Unable to create file '$f': $!\n";

    print OUT "<?xml\n  version = \"1.0\"\n  encoding = \"UTF-8\"?>\n";
    print OUT "<RUN_SET\n" . $opts{xmlns_run};

    #   For all bams in our list, create a RUN section
    my $count = 0;
    my $xml = '';
    my $files2send = '';
    foreach my $bamid (@{$bamidlistref}) {
        my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE bamid=$bamid";
        my $sth = DoSQL($sql);
        my $rowsofdata = $sth->rows();
        $href = $sth->fetchrow_hashref;
        my $s = '';
        #   If the original BAM (cram) has not been sent, include that
        if ($href->{cramorigsent} ne 'Y') {
            $count++;
            $xml .= GenRUNXML($href->{cramchecksum},
                "Secondary mapping Build $opts{build} for $href->{expt_sampleid} [$href->{bamid}]",
                $href->{checksum},
                $centerrunlines,
                $href->{cramname},
                $href->{cramchecksum});
            #   Path of file to send
            $s .= "$opts{backupdir}/$center/$run/$href->{cramname} cramorigsent";
        }
        if ($href->{cramb37sent} ne 'Y') {
            $count++;
            $xml .= GenRUNXML($href->{cramb37checksum},
                "Primary mapping Build $opts{build} for $href->{expt_sampleid} [$href->{bamid}]",
                $href->{checksum},
                $buildrunlines,
                $href->{expt_sampleid} . ".recal.cram",
                $href->{cramb37checksum});
            #   Path of file to send
            $s .= " $opts{resultsdir}/$center/$href->{piname}/" .
                "$href->{expt_sampleid}/bams/$href->{expt_sampleid}.recal.cram cramb37sent";
        }
        if ($s) { $files2send .= $href->{bamid} . " $s\n"; }
    }
    print OUT $xml . "</RUN_SET>\n";
    close(OUT);
    print "  Created '$f' for $count BAMs\n";

    #   Append files to be sent to $path_files2send
    if ($files2send) {
        open(APP, '>>' . $path_files2send) ||
            die "$Script Unable to append to file '$path_files2send': $!\n";
        print APP $files2send;
        close(APP);
    }
}

#==================================================================
# Subroutine:
#   $xml = GenRUNXML($alias, $title, $refname, $runlines, $filename, $checksum)
#
#   Create XML fragment for <RUN> clause
#==================================================================
sub GenRUNXML {
    my ($alias, $title, $refname, $runlines, $filename, $checksum) = @_;

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
        "      filetype = \"cram\"\n" .
        "      checksum_method = \"MD5\"\n" .
        "      checksum = \"$checksum\">\n" .
        "      <READ_LABEL>forward</READ_LABEL>\n" .
        "      <READ_LABEL>reverse</READ_LABEL>\n" .
        "      <READ_LABEL>unmapped</READ_LABEL>\n" .
        "    </FILE>\n" .
        "  </FILES>\n" .
        "</DATA_BLOCK>\n" .
        "<RUN_ATTRIBUTES>\n" .
        "  <RUN_ATTRIBUTE>\n" .
        "    <TAG>assembly</TAG>\n" .
        "    <VALUE>GRCh37</VALUE>\n" .
        "  </RUN_ATTRIBUTE>\n" .
        "</RUN_ATTRIBUTES>\n" .
        "</RUN>\n";
    return $xml;
}

#==================================================================
# Subroutine:
#   $aref = CreateExpt($f, $bamidlistref)
#
#   Create XML file of experiment information
#==================================================================
sub CreateExpt {
    my ($f, $bamidlistref, $path_files2send) = @_;

    #   Create the experiment XML file
    open(OUT, '>' . $f) ||
        die "Script Unable to create file '$f': $!\n";

    print OUT "<?xml\n  version = \"1.0\"\n  encoding = \"UTF-8\"?>\n";
    print OUT "<EXPERIMENT_SET\n" . $opts{xmlns_experiment};
            
    #   Make an EXPERIMENT section for each BAMID
    my $count = 0;
    foreach my $bamid (@{$bamidlistref}) {
        my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE bamid=$bamid";
        my $sth = DoSQL($sql);
        $href = $sth->fetchrow_hashref;
        my $x = Experiment($href);
        if (! $x) { die "$Script Failed to create experiment for $href->{bamname} [$href->{bamid}]\n"; }
        print OUT $x;
        $count++;
    }

    print OUT "</EXPERIMENT_SET>\n";
    close(OUT);
    print "  Created '$f' for $count BAMs\n";
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

    if (exists($href->{bam_delivered}) && $href->{bam_delivered} && ($href->{bam_delivered} > 10)) {
        print "  '$href->{bamname}' has already been sent\n";
        return '';
    }

    #   Build the EXPERIMENT XML for this BAM - should never see uninitialized field
    my $xml = "<EXPERIMENT\n" .
        "  alias = \"$href->{checksum}\"\n" .
        "  center_name = \"$centerdesc\"\n" .
        "  broker_name = \"$opts{broker_name}\">\n" .
        "<TITLE>$opts{experiment_title} for $href->{expt_sampleid}</TITLE>\n" .
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
#   CreateSubmit($f, $runf, $exptf)
#
#   Create XML file of submission information
#==================================================================
sub CreateSubmit {
    my ($f, $runf, $exptf) = @_;

    #   Create the submit XML file
    open(OUT, '>' . $f) ||
        die "$Script Unable to create file '$f': $!\n";

    print OUT "<?xml\n  version = \"1.0\"\n  encoding = \"UTF-8\"?>\n";
    print OUT "<SUBMISSION\n" . $opts{xmlns_submission} .
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
        "      source = \"$exptf\"\n" .
        "      schema = \"experiment\"/>\n" .
        "  </ACTION>\n  <ACTION>\n" .
        "  <ADD\n" .
        "      source = \"$runf\"\n" .
        "      schema = \"run\"/>\n" .
        "  </ACTION>\n";
    print OUT "  <ACTION>\n    <PROTECT/>\n  </ACTION>\n</ACTIONS>\n";

    print OUT "</SUBMISSION>\n";
    close(OUT);
    print "  Created '$f'\n";
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
    my $emsg = "$Script XML file '$f' was mal-formated\n";
    if (! $opts{verbose}) { unlink($f); }
    else { $emsg .= " and has been removed.  Specify -v to keep the file around."; }
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

=item B<-bamlist path>

Specifies the path to a file which contains a list of all the BAMIDs to be
sent to NCBI.  Each line must contain a BAMID beginning in column one.
If none is provided, the list of bams is taken from the database
directly. This is almost never the option you want.

=item B<-build nn>

Specifies the build number for the BAMs to be sent
The default is B<37>.

=item B<-design_description string>

Specifies the DESIGN_DESCRIPTION value in the experiment file.
This defaults to B<Equivalent to Illumina TruSeq PCR-free DNA sample prep>.

=item B<-help>

Generates this output.

=item B<-inst_model>

Specifies the INSTRUMENT_MODEL value in the experiment file.
This defaults to B<HiSeq X Ten>.

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

=item B<-sendlist path>

Specifies the path of a file of data to be sent to NCBI.
Each line of the file contains the parameters passed to the shell script B<topmed_ncbi.sh>.
Two formats for lines in this file are:

bamid + multiple pairs of filetosend databasecolumn 
  or
bamid XML filetosend1 filetosend2 ...

=item B<-verbose N>

Provided for developers to see additional information.

=item B<-xmlprefix string>

Specifies a prefix for the three XML files to be created.
This allows the caller to make these unique.
This defaults to B<nhlbi.$center.$run>.


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

