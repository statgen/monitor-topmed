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
    password => 'please',               # Password to be entered for --remove  (yes, not secure)
    verbose => 0,
);
my $dbh;

Getopt::Long::GetOptions( \%opts,qw(
    help finddups realm=s remove=s verbose
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
#   Special admin hook - show all duplicates in the database
#--------------------------------------------------------------
#   Show all the duplicates. There should be none, but ...
if ($opts{finddups}) {
    my $aref = Dups();
    if ($#$aref < 0) { print "$Script - No duplicates found!\n"; exit; }
    print "Found " . ($#$aref+1) . " duplicates\n";
    if ($opts{verbose}) { print join(' ', @$aref) . "\n"; }
    exit;
}

#   Generate the SQL to remove duplicate records
if ($opts{remove}) {
    if ($opts{remove} ne $opts{password}) {
        die "$Script - password '$opts{remove}' invalid\n";
    }

    my $aref = Dups();
    if ($#$aref < 0) { print "$Script - No duplicates found!\n"; exit; }
    foreach my $bamid (@$aref) {
        my $sql = "SELECT id from $opts{qcdata_table} WHERE bam_id=$bamid ORDER by id";
        my $sth = DoSQL($sql);
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) { print "$Script - Oops, now there is no duplicate '$bamid'?\n"; next; }
        my @dupbamids = ();
        for (my $i=0; $i<$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref();
            push @dupbamids,$href->{id};            
        }
        my $keepbamid = pop(@dupbamids);
        if ($opts{verbose}) { print "Keep $keepbamid, remove " . join(',',@dupbamids) . "\n"; }
        else {
            foreach my $id (@dupbamids) {
                print "DELETE FROM $opts{qcdata_table} WHERE id=$id;\n";
            }
        }
    }
    exit;
}

#--------------------------------------------------------------
#   Extract qplot values or die trying
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
#   See if the NWDID we were provided is valid
my $sth = PersistDoSQL($opts{realm},"SELECT bamid FROM $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'");
my $rowsofdata = $sth->rows();
if (! $rowsofdata) { die "$Script - NWDID '$nwdid' is unknown\n"; }
my $href = $sth->fetchrow_hashref;
my $bamid = $href->{bamid};

my $date = shift(@metrics);
$_ = shift(@colnames);              # Ignore name of this first column
my $sql;

#   See if there already is a qc result for this NWDID
$sth = PersistDoSQL($opts{realm},"SELECT bam_id FROM $opts{qcdata_table} WHERE bam_id=$bamid");
$rowsofdata = $sth->rows();
if ($rowsofdata) {                  # This already exists, do UPDATE
    $sql = "UPDATE $opts{qcdata_table} SET created_at='$date',";
    for (my $i=0; $i<=$#metrics; $i++) {        # No need to specify bam_id
        $sql .= "$colnames[$i]=$metrics[$i],";
    }
    chop($sql);
    $sql .= " WHERE bam_id=$bamid";
}
else {                              # First time, do INSERT
    #   First value is a string, rest are floats
    $sql = "INSERT INTO $opts{qcdata_table} (bam_id,created_at," . lc(join(',',@colnames)) .
        ") VALUES($bamid,'$date',";
    for (my $i=0; $i<=$#metrics; $i++) {
        $sql .= $metrics[$i] . ',';
    }
    chop($sql);
    $sql .= ')';
}
if ($opts{verbose}) { print "SQL=$sql\n"; }
$sth = PersistDoSQL($opts{realm},$sql);
if (! $sth) { die "$Script Failed to update database. SQL=$sql\n"; }

print "Updated QC database data for '$nwdid'\n";
exit;

#==================================================================
# Subroutine:
#   $aref = Dups()
#
#   Return a reference to an array of duplicate qc_result entries
#==================================================================
sub Dups {
    #my ($href, $ncbihref) = @_;
    my @dups = ();
    
    $dbh = DBConnect($opts{realm}) ||
        die "$Script - Unable to connect to database, realm=$opts{realm}\n";
    my $sth = DoSQL("SELECT count(bam_id) as count,bam_id FROM $opts{qcdata_table} group by bam_id having count > 1");
    my $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        my $href;
        for (my $i=0; $i<$rowsofdata; $i++) {
            $href = $sth->fetchrow_hashref();
            push @dups, $href->{bam_id};
        }
    }
    return \@dups;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedqplot.pl - Update the database with qplot values for a NWDID sampleid

=head1 SYNOPSIS

  topmedqplot.pl /net/topmed/incoming/qc.results/broad/2015may01.single.fram NWD321439

  topmedqplot.pl -finddups  ignore ignore       # Debugging, show number qc_results duplicates
  topmedqplot.pl -finddups -verbose ignore ignore 
  topmedqplot.pl -remove=please ignore ignore   # Remove show qc_results duplicates

=head1 DESCRIPTION

This program extracts values of interest from qplot output and saves them
in the database.

See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-finddups>

Scan the database of qc results and show instances of multiple bam_id 
This is only for admins to fing bugs.
No processing of normal parameters is done (e.g. add/modify database).

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-remove=password>

Remove instances of multiple bam_id records from the database.
This generates the SQL to remove all duplicates (all but the last record).
This is only for admins to fing bugs.
No processing of normal parameters is done.

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

