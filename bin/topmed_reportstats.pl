#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name: topmed_reportstats.pl   - completely unused
#
# Description:
#   Use this program to automatically request actions on data
#   from NHLBI TopMed centers.
#   This program is expected to run as a crontab job.
#   If the center contains the file 'donotmonitor', we ignore it.
#
# ChangeLog:
#   $Log: topmed_monitor.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use File::Basename;
use Cwd;
use Getopt::Long;
use DBIx::Connector;

die "This is dead code. Irina chose to generate report data using her own scripts and methods   June 2015\n";

my($me, $mepath, $mesuffix) = fileparse($0, '\.pl');
(my $version = '$Revision: 1.0 $ ') =~ tr/[0-9].//cd;

#   Map qlpot file key to database column name
my %key2col = (
    'Stats\BAM' => 'bamfile_remapped',
    'TotalReads(e6)' => 'totalreads',
    'MappingRate(%)' => 'mappingrate',
    'MapRate_MQpass(%)' => 'maprate_mqpass',
    'TargetMapping(%)' => 'targetmapping',
    'ZeroMapQual(%)' => 'zeromapqual',
    'MapQual<10(%)' => 'mapqual',
    'PairedReads(%)' => 'pairedreads',
    'ProperPaired(%)' => 'properpair',
    'MappedBases(e9)' => 'mappedbases',
    'Q20Bases(e9)' => 'q20bases',
    'Q20BasesPct(%)' => 'q20basespct',
    'MeanDepth' => 'meandepth',
    'GenomeCover(%)' => 'genomecover',
    'EPS_MSE' => 'eps_mse',
    'EPS_Cycle_Mean' => 'eps_cycle_mean',
    'GCBiasMSE' => 'gcbiasmse',
    'ISize_mode' => 'isize_mode',
    'ISize_medium' => 'isize_medium',
    'SecondaryRate(%)' => '',
    'DupRate(%)' => 'duptate',
    'QCFailRate(%)' => 'qcfailrate',
    'BaseComp_A(%)' => 'base_comp_a',
    'BaseComp_C(%)' => 'base_comp_c',
    'BaseComp_G(%)' => 'base_comp_g',
    'BaseComp_T(%)' => 'base_comp_t',
    'BaseComp_O(%)' => 'base_comp_o',
    'Depth>=1(%)' => '',
    'Depth>=5(%)' => '',
    'Depth>=10(%)' => '',
    'Depth>=15(%)' => '',
    'Depth>=25(%)' => '',
    'Depth>=30(%)' => ''
);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    report_table => 'reportdata',
    verbose => 0,
);
my ($dbc);                      # For access to DB
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 1 || $opts{help}) {
    warn "$me$mesuffix [options] bamid statsfile\n" .
        "\nVersion $version\n" .
        "Collect qplot stats for a bam and save in the table $opts{report_table}.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $bamid = shift(@ARGV);
my $statsfile = shift(@ARGV);

my $dbh = DBConnect($opts{realm});
if ($opts{verbose}) { print "$me$mesuffix Version $version, realm=$opts{realm}\n"; }

#--------------------------------------------------------------
#   Read the stats file and extract data we want for database
#--------------------------------------------------------------
my %data = (qplotfile => $statsfile, bamid => $bamid);   # Collect all data here
open(IN, $statsfile) ||
    die "Unable to open file '$statsfile': $!\n";
while (<IN>) {
    s/\t/ /g;
    my @c = split(' ', $_);
    if (! exists($key2col{$c[0]})) {
        warn "Ignoring unknown qplot stats field: $c[0]\n";
        next;
    }
    if ($key2col{$c[0]} eq '') {
        #warn "Not interested in qplot stats field: $c[0]\n";
        next;
    }
    $data{$key2col{$c[0]}} = $c[1];     # Save value for this column
}
close(IN);

#   If necessary, create a new row, then update it
my $reportid = -1;
my $sql = "SELECT reportid from $opts{report_table} WHERE bamid=$bamid";
my $sth = DoSQL($sql);
my $rowsofdata = $sth->rows();
if (! $rowsofdata) {                # Does not exist, create it
    $sql = "INSERT INTO $opts{report_table} SET bamid=$bamid";
    $sth = DoSQL($sql);
    $reportid = $sth->{'mysql_insertid'};
}
else {
    my $href = $sth->fetchrow_hashref;
    $reportid = $href->{reportid};
}

#   Now generate an UPDATE for the entry
$sql = "UPDATE $opts{report_table} SET ";
foreach my $c (keys %data) { $sql .= "$c='$data{$c}',"; }
chop($sql);                         # Drop trailing comma
$sql .= " WHERE reportid=$reportid";
$sth = DoSQL($sql);
print "Updated QPLOT stats for '$bamid' [id=$reportid]\n";
exit;

#==================================================================
# Subroutine:
#   DBConnect($realm)
#
#   Connect to our database using realm '$realm'. Return a DB handle.
#   Get the connection information from DBIx::Connector
#   Fully qualified realm file may be provided
#==================================================================
sub DBConnect {
    my ($realm) = @_;
    if (! $realm) { return 0; }
    #   Get the connection information from DBIx::Connector
    #   Fully qualified realm file may be provided
    if ($realm =~ /^(\/.+)\/([^\/]+)$/) {
        my $d = $1;
        $realm = $2;
        $dbc = new DBIx::Connector(-realm => $realm, -connection_dir => $d,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    else {
        $dbc = new DBIx::Connector(-realm => $realm,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    return $dbc->connect();
}

#==================================================================
# Subroutine:
#   DoSQL - Execute an SQL command
#
# Arguments:
#   sql - String of SQL to run
#   die - boolean if we should die on error
#
# Returns:
#   SQL handle for subsequent MySQL actions
#   Does not return if error detected
#==================================================================
sub DoSQL {
    my ($sql, $die) = @_;
    if (! defined($die)) { $die = 1; }
    if ($opts{verbose}) { warn "DEBUG: SQL=$sql\n"; }
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    if ($DBI::err) {
        if (! $die) { return 0; }
        die "$me$mesuffix SQL failure: $DBI::errstr\n  SQL=$sql\n";
    }
    return $sth;
}
#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_reportstats.pl - Collect QPLOT stats for a bam

=head1 SYNOPSIS

  topmed_reportstats.pl bamid  statsfile

=head1 DESCRIPTION

Use program as to extract QPLOT data of interest for BAM file and 
so that a report can be generated of the results.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<bamid>

The database identifier for the BAM file of interest.

=item B<statsfile>

A file of QPLOT output for the BAM file.

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

