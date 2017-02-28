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
use IO::File;
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
use My_DB;

use POSIX qw(strftime tmpnam);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    verbose => 0,
    max => 100,
    localcachefile => "/tmp/topmed_gce_status.localcache.txt",
    gcecachefileprefix => "gce_status",
    gsuri => 'gs://topmed-recabs/\*/\*.flagstat',
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
);

Getopt::Long::GetOptions( \%opts,qw(
    help dry-run max=i verbose some list cache nocache slow
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] check|mark|verify\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;
my $schema = Topmed::DB->new();

my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'check')    { CheckState(@ARGV); exit; }
if ($fcn eq 'mark')     { MarkState(@ARGV); exit; }
if ($fcn eq 'verify')   { VerifyState(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   MarkState()
#
#   Find all cases of remapped samples at Google and mark
#   them to be pulled.  This creates a cache file of flagstats.
#==================================================================
sub MarkState {
    #my ($x) = @_;
    print "$nowdate MarkState\n";

    my $nwdref = CacheGSData();
    my $count = 0;
    my $countcannot = 0;
    my $cannotlist = '';
    foreach my $nwdid (keys %{$nwdref}) {
        my $sample = $schema->resultset('Bamfile')->find_by_nwdid($nwdid);
        #   Nwdid is set, push should be set. If not, force it
        #   since we do so many manual pushes of samples to GCE
        my $push = $sample->state_gce38push;
        if ($push != $COMPLETED) {
            $sample->update({state_gce38push => $COMPLETED});
            print "Forced push for $nwdid\n";
        }

        #   Only set state_gce38pull if it is NOTSET or has failed
        my $pull = $sample->state_gce38pull;
        my $post = $sample->state_gce38post;
        if (! defined($pull)) {
            print "For some reason getting gce38pull failed for '$nwdid'\n";
            next;
        }
        if ($pull == $NOTSET || $pull == $FAILED) {
            if (! $opts{'dry-run'}) {
                $sample->update({state_gce38pull => $REQUESTED});
                if ($opts{verbose}) { print "Requesting $nwdid\n"; }
            }
            else { print "Dry-Run: Requesting $nwdid\n"; }
            $count++;
            if ($count >= $opts{max}) {
                if ($opts{verbose}) { print "Maximum requests reached\n"; }
                last;
            }
        }
        else {
            $countcannot++;
            $cannotlist .= "$nwdid ($push/$pull/$post) ";
        }
    }
    print "Requested $count pulls to be done, $countcannot NWDIDs could not do pull\n";
    if ($opts{verbose}) { print "CANNOTPULL = $cannotlist\n"; }
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

    print "$nowdate VerifyState\n";

    my $nwdref = CacheGSData();
    my $ahref = GetSQLData();        # Get array of nwdid -> hash
    print "Checking " . scalar(keys %{$ahref}) . " database samples\n";

    my %counts = ();
    my @didnotdelete = ();
    my @postfailed = ();
    my @pullfailed = ();
    my @unknownnwdid = ();
    my @nopushdone = ();
    my @rmfailed = ();
    foreach my $nwdid (keys %{$nwdref}) {
        $counts{'Completed recabs'}++;
        my $href = $ahref->{$nwdid};
        if (! defined($href)) {
            my $sample = $schema->resultset('Bamfile')->find_by_nwdid($nwdid);
            if (! $sample) {
                print "Unable to find '$nwdid' in database\n";
                next;
            }
            #   Database had successful push/pull/post
            #   See if remapped sample is here
            my $f = $sample->b38_mapped_path . "/$nwdid.recab.cram";
            my $fflag = $sample->b38_mapped_path . "/$nwdid.recab.cram.flagstat";
            my @s = stat($f);
            my $cramsize = $s[7];
            @s = stat($f);
            my $flagsize = $s[7];
            if ($flagsize < 100 || $cramsize < 1000000) {   # We have local GCE data   Huh??
                $counts{'Unknown NWDID'}++;
                push @unknownnwdid, $nwdid;
            }
            else {                                  # We have local file, RM at GCE failed
                $counts{'RM_at_GCE_failed'}++;
                push @rmfailed, $nwdid;
            }
            next;
        }
        my $pushstate = $href->{state_gce38push};
        my $pullstate = $href->{state_gce38pull};
        my $poststate = $href->{state_gce38post};
        #   If flagstat exists, sample is ready to be pulled
        #   See what the database thinks has/will happened
        if ($pushstate != $COMPLETED) {
            $counts{'PUSH not done'}++;
            push @nopushdone,$nwdid;
            next;
        }

        if ($pullstate == $COMPLETED) {
            $counts{'PULL finished'}++;
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

        if ($poststate == $COMPLETED) {
            $counts{'POST did not delete files'}++;
            push @didnotdelete,$nwdid;
        }
        if ($poststate == $FAILED || $poststate == $CANCELLED) {
            $counts{'POST failed'}++;
            push @postfailed,$nwdid;
        }
        if ($poststate == $REQUESTED || $poststate == $SUBMITTED || $poststate == $STARTED) {
            $counts{'POST in progress'}++;
        }
        if ($poststate == $NOTSET) {
            $counts{'POST yet to be requested'}++;
        }
    }

    # Give summary of what was found
    foreach my $k (sort keys %counts) {
        print "$k $counts{$k} times\n";
    }
    print "\n";
    if (@didnotdelete) {
        print scalar(@didnotdelete) . " POST failed to delete for these: " .
            substr(join(' ',@didnotdelete),0,50) . " ...\n";
        if ($opts{verbose}) { print '  ' . join("\n  ", @didnotdelete) . "\n\n "; }
    }
    if (@postfailed) {
        print scalar(@postfailed) . " POST failed for these: " . 
            substr(join(' ',@postfailed),0,50) . " ...\n";
        if ($opts{verbose}) {
            foreach my $nwd (@postfailed) {
                print $nwd . ' ';
                my $href = $ahref->{$nwd};
                my $s = `tail -1 ~/output/$href->{bamid}-gcepost.out 2>/dev/null`;
                chomp($s);
                print "$nwd $s\n";
            }   
            print "\n ";
        }
    }
    if (@pullfailed) {
        print scalar(@pullfailed) . " PULL failed for these: " .
            substr(join(' ',@pullfailed),0,50) . " ...\n";
        if ($opts{verbose}) {
            foreach my $nwd (@pullfailed) {
                print $nwd . ' ';
                my $href = $ahref->{$nwd};
                my $s = `tail -1 ~/output/$href->{bamid}-gcepull.out 2>/dev/null`;
                chomp($s);
                print "$nwd $s\n";
            }   
            print "\n ";
        }
    }
    if (@nopushdone) {
        print scalar(@nopushdone) . " No PUSH was done, yet file at GCE: " .
            substr(join(' ',@nopushdone),0,50) . " ...\n";
        if ($opts{verbose}) { print '  ' . join("\n  ", @nopushdone) . "\n\n "; }
    }
    if (@rmfailed) {
        print scalar(@rmfailed) . " Might have failed to remove GCE files for these: " .
            substr(join(' ',@rmfailed),0,50) . " ...\n";
        if ($opts{verbose}) {
            foreach my $nwd (@rmfailed) {
                print $nwd . "   ls -l `/usr/cluster/topmedpath.pl wherefile $nwd b38`" .
                    " && gsutil rm -rf gs://topmed-recabs/$nwd/\n";
            }   
            print "\n ";
        }
    }
    if (@unknownnwdid) {
        print scalar(@unknownnwdid) . " Database marked finished for these: " .
            substr(join(' ',@unknownnwdid),0,50) . " ...\n" .
            "   Perhaps POST failed to delete this at GCE?\n";
        if ($opts{verbose}) { print '  ' . join("\n  ", @unknownnwdid) . "\n\n "; }
    }

}

#==================================================================
# Subroutine:
#   CheckState()
#
#   Checks the consistency of remapped files (local and at Google)
#   and the state of the database flags.
#==================================================================
sub CheckState {
    #my ($x) = @_;
    my @hosts = qw{topmed9 topmed10};   # Where local recab files live
    my %counts = ();                    # Used to generate counts of everything
    my %countslist = ();                # Collect lists of NWDIDs here

    print "$nowdate Getting local remapped files\n";
    my %localfiles = ();
    foreach my $h (@hosts) {
        my $n = -1;
        if ($opts{some}) { $n = 100; }
        my $cmd = "find /net/$h/working/mapping/results -name *.recab.cram -print";
        my $in = new IO::File;
        if ($opts{cache} && -f $opts{localcachefile}) {
            $cmd = $opts{localcachefile};
            print "Reading GCE samples from $opts{localcachefile}\n";
        }
        else { $cmd = "$cmd |"; }
        $in->open($cmd) ||
            die "$Script - Unable to get list of Google remapped files.  CMD=$cmd\n";
        while (<$in>) {
            if (/\/(NWD\d*)\.recab\.cram$/) { $localfiles{$1} = 1; }
            else { print "$Script - cannot parse $_"; }
            $n-- || last;
        }
        $in->close();
    }
    print "Found " . keys(%localfiles) . " local remapped files\n";
    if ($opts{cache} && (! -f $opts{localcachefile})) {  # Create cache file one time only
        my $out = IO::File->new($opts{localcachefile}, 'w');
        foreach (keys %localfiles) { print $out '/' . $_ . ".recab.cram\n"; }
        close($out);
        print "Created cache file of local samples in $opts{localcachefile}\n";
    }

    print "Getting remapped files at GCE\n";
    my $cmd = 'gsutil ls gs://topmed-recabs/\*/\*.recab.cram';
    if ($opts{verbose}) { print "This could be slow, patience ...\n"; }
    my $n = -1;
    if ($opts{some}) { $n = 100; }
    my $in = new IO::File;
    if ($opts{cache} && -f $opts{gcecachefile}) {
        $cmd = $opts{gcecachefile};
        print "Reading GCE samples from $opts{gcecachefile}\n";
    }
    else { $cmd = "$cmd |"; }
    $in->open($cmd) ||
        die "$Script - Unable to get list of Google remapped files.  CMD=$cmd\n";
    my %gcefiles = ();
    while (<$in>) {
        if (/\/(NWD\d*)\.recab\.cram$/) { $gcefiles{$1} = 1; }
        else { print "$Script - cannot parse $_"; }
        $n-- || last;
    }
    $in->close();
    print "Found " . keys(%gcefiles) . " GCE remapped files\n";
    if ($opts{cache} && (! -f $opts{gcecachefile})) {  # Create cache file one time only
        my $out = IO::File->new($opts{gcecachefile}, 'w');
        foreach (keys %gcefiles) { print $out '/' . $_ . ".recab.cram\n"; }
        close($out);
        print "Created cache file of GCE samples in $opts{gcecachefile}\n";
    }

    #   Check database status for all GCE files
    foreach my $nwdid (keys %gcefiles) {
        my $sample = $schema->resultset('Bamfile')->find_by_nwdid($nwdid);
        if (! $sample) {
            print "Ignoring GCE sample $nwdid\n";
            next;
        }
        $counts{GCE_samples}++;

        #   Sample is at GCE
        my $pushstate = $sample->state_gce38push;
        my $pullstate = $sample->state_gce38pull;
        my $poststate = $sample->state_gce38post;
        if ($pushstate == $COMPLETED) {
            #   Sample is at GCE and was pushed
            $counts{GCE_samples_pushed}++;
            $countslist{GCE_samples_pushed} .= $nwdid . ' ';
            if ($pullstate == $NOTSET || $pullstate == $FAILED) {
                $counts{GCE_samples_pushed_and_no_pull}++;
                $countslist{GCE_samples_pushed_and_no_pull} .= $nwdid . ' ';
            }
            if ($pullstate == $COMPLETED) {
                $counts{GCE_samples_pushed_and_pulled}++;
                $countslist{GCE_samples_pushed_and_pulled} .= $nwdid . ' ';
            }
            if ($poststate == $NOTSET || $pullstate == $FAILED) {
                $counts{GCE_samples_pushed_and_no_post}++;
                $countslist{GCE_samples_pushed_and_no_post} .= $nwdid . ' ';
            }
            if ($poststate == $COMPLETED) {
                $counts{GCE_samples_pushed_and_post_done_but_not_removed}++;
                $countslist{GCE_samples_pushed_and_post_done_but_not_removed} .= $nwdid . ' ';
            }
        }
        else {
            #   Sample exists, but we did not push it
            $counts{GCE_samples_NOT_pushed}++;
            $countslist{GCE_samples_NOT_pushed} .= $nwdid . ' ';
        }       
    }

    #   Check database status for all local files
    foreach my $nwdid (keys %localfiles) {
        my $sample = $schema->resultset('Bamfile')->find_by_nwdid($nwdid);
        if (! $sample) {
            print "Ignoring local sample $nwdid\n";
            next;
        }
        $counts{LOCAL_samples}++;

    }

    #   Show results
    my $m;
    print "\n";
    foreach my $msg (sort keys %counts) {
        $m = $msg;
        $m =~ s/_/ /g;
        print "$m = $counts{$msg}\n";
        if ($opts{list}) {
            if (exists($countslist{$msg})) {
                print '   ' . $countslist{$msg} . "\n\n";
            }
        }
    }
}

#==================================================================
# Subroutine:
#   GetSQLData()
#
#   Execute the SQL
#   Return a reference to an array of hashes of the SQL data
#==================================================================
sub GetSQLData {
    #my ($x) = @_;

    DBConnect($opts{realm});
    my $sql = "SELECT bamid,expt_sampleid,state_gce38push," .
        "state_gce38pull,state_gce38post FROM $opts{bamfiles_table}" .
        " WHERE state_gce38push!=$COMPLETED " .
        " OR state_gce38pull!=$COMPLETED OR " .
        " state_gce38post!=$COMPLETED";
    my $sth = DoSQL($sql);
    my $ahref = $sth->fetchall_hashref('expt_sampleid');     # Returns array of hash
    return $ahref;

    #   Get rows one by one when we don't trust fetchall_hashref to work
    $sql = "SELECT bamid,expt_sampleid,state_gce38push," .
        "state_gce38pull,state_gce38post FROM $opts{bamfiles_table}";
    $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %data = ();
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if (($href->{state_gce38push} != 20) ||
            ($href->{state_gce38pull} != 20) ||
            ($href->{state_gce38post} != 20)) {
            $data{$href->{expt_sampleid}} = $href;
        }
    }
    $ahref = \%data;
    return $ahref;

    #   Dump some data
    my $n = 10;
    foreach my $nwdid (keys %{$ahref}) {
        my $href = $ahref->{$nwdid};
        print "bamid=$href->{bamid} nwdid=$href->{expt_sampleid} ($nwdid) " .
            " push=$href->{state_gce38push} pull=$href->{state_gce38pull} " .
            " post=$href->{state_gce38post}\n";
        $n--;
        if ($n <= 0) { last; }
    }
    return $ahref;
}

#==================================================================
# Subroutine:
#   CacheGSData()
#
#   Get all flagstat files at GCE. Cache them in a local file
#   and return a reference to a hash of NWDID and the GCE file
#==================================================================
sub CacheGSData {
    #my ($x) = @_;
    my %gcedata = ();

    #   Look for an existing cache file that isn't too old
    #   and either use it or remove it and create a new one
    my $now = time();
    my ($out, $in);
    opendir($in, '/tmp') ||
        die "$Script - cannot read directory '/tmp': $!\n";
    while (readdir $in) {
        if (! /$opts{gcecachefileprefix}\.(\d+)/) { next; }
        my $t = $1;
        if (($now - $t) < 60*60*24) { $now = $t; last; }
        else {              # Old cache file, remove it
            my $f = "/tmp/$opts{gcecachefileprefix}.$t"; 
            unlink($f);
            print "Removing old cache file '$f'\n";
        }
    }
    closedir($in);

    #   If a cache exists, read it and return the data.  $now is time for file
    my $count = 0;
    my $flagstatscache = "/tmp/$opts{gcecachefileprefix}." . $now;
    if ((! $opts{nocache}) && -r $flagstatscache) {
        print "Reading GCE data from '$flagstatscache'\n";
        open($in, $flagstatscache) ||
            die "$Script - Unable to read file '$flagstatscache'; $!\n";
        while (<$in>) {
            my @c = split(' ', $_);
            if ($c[0] eq 'NWD000') { next; }    # GCE test sample
            $gcedata{$c[0]} = $c[1];
            $count++;
        }
        close($in);
        print "Returned $count GCE samples from cache file\n";
        return \%gcedata;
    }

    #   No cache file, get list of files from GCE    
    open($out, '>' . $flagstatscache) ||
        die "$Script - Unable to create file '$flagstatscache'; $!\n";
    my $cmd = "gsutil ls $opts{gsuri}";    
    print "Very slowly creating cache file '$flagstatscache'\n";
    open($in, "$cmd |") ||
        die "$Script - Unable to get list of flagstat files.  CMD=$cmd\n";
    while (my $l = <$in>) {
        chomp($l);
        if ($l !~ /^gs.*\/(NWD\d+)\/NWD\d*\.recab\.cram\.flagstat/) {
            if ($opts{verbose}) { print "$Script - cannot parse $l\n"; }
            next;
        }
        if ($1 eq 'NWD000') { next; }
        $gcedata{$1} = $l;
        print $out "$1 $l\n";           # Save data in cache file
        $count++;
    }
    close($in);
    close($out);
    print "Created cache file '$flagstatscache' of $count samples\n";
    return \%gcedata;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_gce_pull.pl - Manage the state_gce38* flags and remapped files

=head1 SYNOPSIS
 
  topmed_gce_pull.pl check                  # Check files and data for consistency
  topmed_gce_pull.pl mark                   # Sets state_gce38pull in database
  topmed_gce_pull.pl verify                 # Compares files at Google and database state

=head1 DESCRIPTION

This program is used to manage the state_gce38* flags in the database.

=head1 OPTIONS

=over 4

=item B<-cache>

A list of GCE files is saved in B</tmp/topmde_gce_status.cache.txt> the first time.
Subsequent calls will load the GCE file list from this.
This is useful for testing.
Don't forget to remove the cache file.

=item B<-dry-run>

Does not actually change values in the database.

=item B<-help>

Generates this output.

=item B<-list>

Generates a list of all samples for each count.

=item B<-max N>

When B<set> is used, only sets this many state flags.

=item B<-some>

Only read some of the files at Google or in local store. Useful for testing.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program can be:

B<check>
Checks the database flags to see if they properly reflect the state of files
at Google and in the local file system for remapped files.
=back

B<mark>
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
