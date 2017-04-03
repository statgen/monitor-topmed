#!/usr/bin/perl
###################################################################
#
# Name: topmed_status.pl
#
# Description:
#   Use this program to automatically update the states for runs
#   from NHLBI TopMed centers.
#   This program is expected to run as a crontab job.
#
# ChangeLog:
#   $Log: topmed_status.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use FindBin qw($Bin $Script);
use lib (
  qq($FindBin::Bin),
  qq($FindBin::Bin/../lib),
  qq($FindBin::Bin/../lib/perl5),
  qq($FindBin::Bin/../local/lib/perl5),
  qq(/usr/cluster/topmed/lib/perl5),
  qq(/usr/cluster/topmed/local/lib/perl5),
);
use My_DB;
use TopMed_Get;
use Getopt::Long;

#   How to add a column to this program
#       Add db column => 'letter'

#   These represent states the run can be in. column name to letter
my %attributes2letter = (
    'state_arrive'   => 'A',
    'state_verify'   => '5',
    'state_qplot'    => 'Q',
    'state_cram'     => 'C',
    'state_b37'      => '7',
    'state_gce38push'=> 's',
    'state_gce38pull'=> 'r',
    'state_gce38post'=> 'Z',
    'state_b38'      => '8',
    'state_gce38bcf_push'=> 'x',
    'state_gce38bcf_pull'=> 'y',
    'state_bcf'      => 'B',
    'state_ncbiexpt' => 'X',
    'state_ncbiorig' => 'S',
    'state_ncbib37'  => 'P',
);

my $NOTSET    = 0;            # Not set
my $REQUESTED = 1;            # Task requested
my $SUBMITTED = 2;            # Task submitted to be run
my $STARTED   = 3;            # Task started
my $DELIVERED = 19;           # Data delivered, but not confirmed
my $COMPLETED = 20;           # Task completed successfully
my $CANCELLED = 89;           # Task cancelled
my $FAILEDCHECKSUM = 98;      # Task failed, because checksum at NCBI bad
my $FAILED    = 99;           # Task failed

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] runstatus\n" .
        "Update the status fields in the database.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift(@ARGV);

my $dbh = DBConnect($opts{realm});

#--------------------------------------------------------------
#   Calculate overall status for each run 
#--------------------------------------------------------------
if ($fcn eq 'runstatus') {
    #   Get all the known centers in the database
    my @attributes = keys %attributes2letter;
    my $attributes_str = join(',', @attributes);
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        #   For each run, calculate a status for all the bam files in the run
        my $runsref = GetRunsCount($cid);
        if (! %$runsref) { next; }          # No directories we want
        foreach my $runid (keys %{$runsref}) {
            my $numberbams = $runsref->{$runid};
            #   Get list of all bams for this run
            my $sql = "SELECT bamid,$attributes_str FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            my %counts = ();
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;      # Data for this BAM
                #   Get counts of states we are interested for each attribute
                foreach my $a (@attributes) {
                    if (! defined($href->{$a})) { next; }
                    if ($href->{$a} == $NOTSET) { next; }
                    if ($href->{$a} == $CANCELLED)  { next; }    # cancelled
                    #   These go into done, processing or neither states
                    if ($href->{$a} ==  $COMPLETED)  { $counts{$a}{done}++;   next; } # yea!
                    if ($href->{$a} >= $FAILEDCHECKSUM) { $counts{$a}{failed}++; next; } # boo!
                    $counts{$a}{processing}++;
                }
            }

            #   Figure out the status for this run
            #   Status consists of a letter and a summary value
            #   The letter represents the step (e.g. A for arrived)
            #   the summary value is one of done processing or unknown
            my $s = '';
            foreach my $a (@attributes) {
                if (! exists($counts{$a}{failed}))     { $counts{$a}{failed} = 0; }
                if (! exists($counts{$a}{done}))       { $counts{$a}{done} = 0; }
                if (! exists($counts{$a}{processing})) { $counts{$a}{processing} = 0; }
                #   If we have done all these, it's done
                if ($counts{$a}{failed} > 0) {
                    $s .= "$attributes2letter{$a}=failed,";
                    next;
                }
                if ($counts{$a}{processing} > 0) {
                    $s .= "$attributes2letter{$a}=processing,";
                    next;
                }
                if ($counts{$a}{done} > 0 && $counts{$a}{done} < $numberbams) {
                    $s .= "$attributes2letter{$a}=processing,";
                    next;
                }
                if ($numberbams && $counts{$a}{done} >= $numberbams) {
                    $s .= "$attributes2letter{$a}=done,";
                    next;
                }
                $s .= "$attributes2letter{$a}=unknown,";
            }
            if ($opts{verbose}) { warn "Status runid $runid=$s\n"; }
            #   Update status for this run 
            $sql = "UPDATE $opts{runs_table} SET status='$s' WHERE runid=$runid";
            $sth = DoSQL($sql);
        }
    }
    exit;
}

die "Invalid request '$fcn'. Try '$Script --help'\n";

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_status.pl - automatically update the states for runs

=head1 SYNOPSIS

  topmed_status.pl runstatus

=head1 DESCRIPTION

Use this program as a crontab job to update the database for runs
with an overall status of the BAMs for a particular run.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-runs NAME[,NAME,...]>

Specifies a specific set of runs,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<runstatus>

Directs this program to update the status for all centers.

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

