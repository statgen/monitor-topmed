#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
###################################################################
#
# Name: topmedqplot.pl
#
# Description:
#   Use this program to update the NHLBI database with data created by qplot
#
# ChangeLog:
#   $Log: topmedqplot.pl,v $
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
use Parse_QC_Files;
use My_DB;
use Getopt::Long;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    qcdata_table => 'qc_results',
    bamfiles_table => 'bamfiles',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 1 || $opts{help}) {
    warn "$Script [options] dirname nwdid\n" .
        "\nUpdate the topmed database with QPLOT data\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $dirname = shift @ARGV;
my $nwdid   = shift @ARGV;

#--------------------------------------------------------------
#   Extract qplot valies or die trying
#--------------------------------------------------------------
my $f = $dirname . '/' . $nwdid;
my $self  = `ls $f*.vb.selfSM`;
my $stats = `ls $f*.qp.stats`;
my $r     = `ls $f*.qp.R`;
my @metrics = parseQCFiles($self, $stats, $r);

if ($opts{verbose}) { print "$nwdid parsed:\n"; foreach my $v (@metrics) { print "$v\n"; } }

my @colnames = parseQCColumns();
if ($#colnames != $#metrics) {
    die "$Script Number of columns [$#colnames] returned by parseQCFiles [$#metrics] did not match\n";
}

#--------------------------------------------------------------
#   Update database with values we obtained
#--------------------------------------------------------------
#   Our network can get flakey and the DBConnect can fail
#   Catch when this happens and wait a bit and try again
my $dbh;
my $sleeptime = 10;
for (1 .. 10) {
    eval { $dbh = DBConnect($opts{realm}); };
    if ($@) {                           # Failed, wait a bit and try again
        print "Datbase connection failed, wait and retry\n";
        sleep($sleeptime);
        $sleeptime += 10;
    }
    else { last; }
}
if ($@) { die $@ . "\n"; }

#   See if the NWDID we were provided is valid
my $sth = DoSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'", 0);
my $rowsofdata = $sth->rows();
if (! $rowsofdata) { die "$Script - NWDID '$nwdid' is unknown\n"; }
my $href = $sth->fetchrow_hashref;
my $bamid = $href->{bamid};

#   First value returned is a string, rest are floats
my $date = shift(@metrics);
$_ = shift(@colnames);              # Ignore name of this first column
my $sql = "INSERT INTO $opts{qcdata_table} (bam_id,created_at," . lc(join(',',@colnames)) .
    ") VALUES($bamid,'$date',";
for (my $i=0; $i<=$#metrics; $i++) {
    $sql .= $metrics[$i] . ',';
}
chop($sql);
$sql .= ')';
if ($opts{verbose}) { print "SQL=$sql\n"; }
$sth = DoSQL($sql, 0);
if (! $sth) { die "$Script Failed to update database. SQL=$sql\n"; }

print "Updated QC database data for '$nwdid'\n";
exit;

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedqplot.pl - Update the database with qplot values for a NWDID sampleid

=head1 SYNOPSIS

  topmedqplot.pl /net/topmed/incoming/qc.results/broad/2015may01.single.fram NWD321439

=head1 DESCRIPTION

This program extracts values of interest from qplot output and saves them
in the database.

See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<dirname>

The directory where the qplot output resides

=item B<nwdid>

The NWDID for the data to extract. This serves as the prefix for each of the qplot files.

=item B<-verbose>

Provided for developers to see additional information.

=back


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

