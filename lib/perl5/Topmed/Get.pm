package TopMed_Get;
#==================================================================
# TopMed_Get.pm
#   Common functions for the various Topmed Perl programs
#==================================================================
use strict;
use warnings;

use vars qw(@EXPORT_OK @EXPORT @ISA);
use Exporter;

@ISA = qw(Exporter);
@EXPORT_OK = qw(GetCenters GetRuns);
@EXPORT    = qw(GetCenters GetRuns);

#==================================================================
# Subroutine:
#   GetCenters - Get list of centers. If $opts{center} is defined
#     the list returned will only contain that center.
#
# Arguments:
#   none
#
# Returns:
#   Reference to hash of center ids  to center names
#==================================================================
sub GetCenters {
    my %centers = ();

    #   Get all the known centers in the database
    my $sql = "SELECT centerid,centername FROM $main::opts{centers_table}";
    my $sth = My_DB::DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { warn "$main::Script No centers in '$main::opts{topdir}'\n"; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if (defined($main::opts{center}) &&
            $main::opts{center} eq $href->{centername}) {
            my %c = ($href->{centerid} => $href->{centername});
            print "$main::Script - Running on center '$main::opts{center}' only\n";
            return \%c;
        }
        $centers{$href->{centerid}} = $href->{centername};
    }
    if ($main::opts{center}) { die "$main::Script '$main::opts{center}' unknown\n"; }
    return \%centers;
}

#==================================================================
# Subroutine:
#   GetRuns - Get list of runs for a center.  If $opts{runs} is defined
#     the list returned will only those runs.
#
# Arguments:
#   cid = center id
#
# Returns:
#   Reference to hash of run ids to run dirnames
#==================================================================
sub GetRuns {
    my ($cid) = @_;
    my %run2dir = ();

    my $sql = "SELECT runid,dirname FROM $main::opts{runs_table} WHERE centerid=$cid";
    my $sth = My_DB::DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) {
        if ($main::opts{verbose}) { warn "$main::Script Found no runs for center '$cid'\n"; }
        return \%run2dir;
    }
    my %theseruns = ();
    if ($main::opts{runs}) {              # We only want some runs
        $main::opts{runs} =~ s/,/ /g;
        foreach my $r (split(' ', $main::opts{runs})) {
            $theseruns{$r} = 1;
            if (! $main::opts{onlyrun}) { 
                print "$main::Script - Using only run '$r'\n";
                $main::opts{onlyrun}++;
            }
        }
    }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($main::opts{runs} && (! exists($theseruns{$href->{dirname}}))) { next; }
        $run2dir{$href->{runid}} = $href->{dirname};
    }
    return \%run2dir;
}

1;

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

TopMed_Get.pm

=head1 SYNOPSIS

  use TopMed_Get;
  
=head1 DESCRIPTION

This provides common functions used by many TopMed programs.

=head1 Functions

=over 4

=item B<GetCenters>

Returns a reference to an array of the TopMed Centers.

=item B<GetRuns($cid)>

Returns a reference to an array of run directories for a particular center-id.

=back

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2010 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut
