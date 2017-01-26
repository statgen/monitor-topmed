#!/usr/bin/env perl
###################################################################
#
# Name: topmed_gce_pull.pl
#
# Description:
#   Use this program to update the database to request that remapped
#   data in Google Cloud be pulled to local store.
#
# ChangeLog:
#   $Log: topmed_gce_pull.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;

use FindBin qw($Script);
use lib (
  qq($FindBin::Bin),
  qq($FindBin::Bin/../lib),
  qq($FindBin::Bin/../lib/perl5),
  qq($FindBin::Bin/../local/lib/perl5),
  qq(/usr/cluster/topmed/lib/perl5),
  qq(/usr/cluster/topmed/local/lib/perl5),
);

use Topmed::Base qw(cmds);
use Topmed::Constants qw(:states);
use Topmed::DB;
use Getopt::Long;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    verbose => 0,
    max => 100,
);

Getopt::Long::GetOptions( \%opts,qw(
    help dry-run max=i verbose 
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] set|verify\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;
my $schema = Topmed::DB->new();

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'set')      { SetState(@ARGV); exit; }
if ($fcn eq 'verify')   { VerifyState(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   SetState()
#
#   Find all cases of remapped samples at Google and mark
#   them to be pulled
#==================================================================
sub SetState {
    #my ($x) = @_;
    
    my $sample_rs    = $schema->resultset('Bamfile')->find_gce_uploads();
    my $flagstat_ptn = YASF->new('gsutil ls {uri}/{nwdid}.recab.cram.flagstat 2> /dev/null');

    for my $sample ($sample_rs->all) {
        my $nwdid = $sample->expt_sampleid;
        my $cmd = $flagstat_ptn % {uri => $sample->gce_recab_uri, nwdid => $nwdid};
        capture([0..1], $cmd);

        #   Only set state_gce38pull if it is NOTSET
        unless ($EXITVAL) {
            if (! $opts{'dry-run'}) {
                my $state = $sample->state_gce38pull;

                if ($state == $NOTSET) {
                    $sample->update({state_gce38pull => $REQUESTED});
                    if ($opts{verbose}) { say "REQUESTing $nwdid"; }
                }
            }
            else {
                if ($opts{verbose}) { say "Dry-Run: Would REQUEST $nwdid"; }
            }
            $opts{max}--;               # Don't do too much
            if ($opts{max} <= 0) {
                if ($opts{verbose}) { say "Maximum pulls reached"; }
                last;
            }
        }
    }
}

#==================================================================
# Subroutine:
#   VerifyState()
#
#   Verify the database state reflects all the state of remapped
#   samples at Google 
#==================================================================
sub VerifyState {
    #my ($x) = @_;

    my $cmd = 'gsutil ls gs://topmed-recabs/\*/\*.flagstat';    # Get all flagstats for recabs
    if ($opts{verbose}) { print "This could be slow, patience ...\n"; }
    my $in;
    open($in, "$cmd |") ||
        die "$Script - Unable to get list of flagstat files.  CMD=$cmd\n";
    my %counts = ();
    my @didnotdelete = ();
    my @postfailed = ();
    my @pullfailed = ();
    while (my $l = <$in>) {
        if ($l !~ /^gs.*\/(NWD\d+)\/NWD\d*\.recab\.cram\.flagstat/) {
            if ($opts{verbose}) { print "$Script - cannot parse $l"; }
            next;
        }
        $counts{'Completed recabs'}++;
        my $nwdid = $1;
        my $sample = $schema->resultset('Bamfile')->find_by_nwdid($nwdid);
        if (! $sample) {
            if ($opts{verbose}) { print "Ignoring $nwdid\n"; }
            next;
        }
        my $pullstate = $sample->state_gce38pull;
        my $poststate = $sample->state_gce38post;

        #   If flagstat exists, sample is ready to be pulled
        #   See what the database thinks has/will happened
        if ($pullstate == $COMPLETED) {
            if ($poststate == $COMPLETED) {
                $counts{'POST did not delete files'}++;
                push @didnotdelete,$nwdid;
            }
            if ($poststate == $FAILED) {
                $counts{'POST failed'}++;
                push @postfailed,$nwdid;
            }
            next;
        }
        if ($pullstate == $REQUESTED || $pullstate == $SUBMITTED || $pullstate == $STARTED) {
            $counts{'PULL in progress'}++;
        }
        if ($pullstate == $NOTSET) {
            $counts{'PULL yet to be requested'}++;
        }
        if ($pullstate == $FAILED || $pullstate == $CANCELLED) {
            $counts{'PULL failed'}++;
            push @pullfailed,$nwdid;
        }
        #$opts{max}--; if ($opts{max} < 0) { last; }
    }
    close($in);

    # Give summary of what was found
    foreach my $k (sort keys %counts) {
        print "$k $counts{$k} times\n";
    }
    print "\n";
    if (@didnotdelete) {
        print "POST failed to delete for these: " . join(' ',@didnotdelete) . "\n";
    }
    if (@postfailed) {
        print "POST failed for these: " . join(' ',@postfailed) . "\n";
    }
    if (@pullfailed) {
        print "PULL failed for these: " . join(' ',@pullfailed) . "\n";
    }

}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_gce_pull.pl - Manage the state_gce38pull flag and files at Google

=head1 SYNOPSIS
 
  topmed_gce_pull.pl set                    # Sets state_gce38pull in database
  topmed_gce_pull.pl verify                 # Compares files at Google and database state

=head1 DESCRIPTION

This program is used to manage the state_gce38flag in the database.

=head1 OPTIONS

=over 4

=item B<-dry-run>

Does not actually change values in the database.

=item B<-help>

Generates this output.

=item B<-max N>

When B<set> is used, only sets this many state flags.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program can be:

B<set>
Scans the files at Google and for every remapped sample, sets the database
flag so the sample will be pulled to local storage.


B<verify>
Scans the files at Google and verifies the database flag is in the correct state.
=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Chris Scheller I<E<lt>schelcj@umich.eduE<gt>> and
Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut
