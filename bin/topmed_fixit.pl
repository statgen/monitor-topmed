#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
###################################################################
#
# Name: topmed_fixit.pl
#
# Description:
#   Use this program to correct problems in a controlled way
#
# ChangeLog:
#   $Log: topmed_fixit.pl,v $
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
use TopMed_Get;
use Getopt::Long;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my $topmedbin = '/usr/cluster/monitor/bin';
our %opts = (
    topmedcmd => "$topmedbin/topmedcmd.pl",
    topmedcmdnwd => "$topmedbin/topmed_nwdid.pl",
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    studies_table => 'studies',
    bamfiles_table => 'bamfiles',
    topdir  => '/net/topmed/incoming/topmed',
    backupdir => '/working/backups/incoming/topmed',
    resultsdir => '/incoming/qc.results',
    count => 1000,
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s count=i
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] nwd|rename\n" .
        "Fix problems with various hacks.\n";
    exit 1;
}
my $fcn = shift(@ARGV);
my $dbh = DBConnect($opts{realm});

#--------------------------------------------------------------
#   Update the database with all sorts of missing data
#--------------------------------------------------------------
if ($fcn eq 'nwd') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles that arrived
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            print "Doing $dirname\n";
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" .
                    $href->{bamname};
                my @stats = stat($f);
                if (! @stats) {
                    print "BAM '$f' not found, continuing\n";
                    next;
                }
                #   Maybe this was already done?
                if (defined($href->{cramname}) && $href->{cramname} ne '') { next; }
                #   Set database fields
                my $cmd = "$opts{topmedcmdnwd} -bamid $href->{bamid} -nonwdid $f";
                system($cmd) &&
                    print "  Failed: CMD=$cmd\n";
                print "  $href->{bamname} updated\n";
            }
        }
    }
    exit;
}

#--------------------------------------------------------------
#   Rename the BAM file to use NWDID, not whatever the center wanted
#--------------------------------------------------------------
if ($fcn eq 'rename') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles that arrived
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            print "Doing $dirname\n";
            #   Get list of all bams we need to rename
            my $sql = "SELECT bamid,bamname,bamname_orig,expt_sampleid,datemapping FROM $opts{bamfiles_table} WHERE runid=$runid AND bamname NOT LIKE 'NWD%'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) {
                print "  All BAMS start with NWD. Continuing\n";
                next;
            }

            #   Rename MD5 files so they do not get reprocessed
            my $d = $opts{topdir} . "/$centername/$dirname";
            if (opendir(DIR, $d)) {
                my @renamelist = ();
                while (readdir(DIR)) {
                    if (/\.md5$/) { push @renamelist,"$d/$_"; next; }
                    if (/\.MD5$/) { push @renamelist,"$d/$_"; next; }
                    if ($_ eq 'Manifest.txt') { push @renamelist,"$d/$_"; next; }
                }
                closedir DIR;
                foreach my $f (@renamelist) {
                    rename($f, "$f.old") ||
                }
            }

            for (my $i=1; $i<=$rowsofdata; $i++) {
                #   Convenience variables
                my $href = $sth->fetchrow_hashref;
                my $oldbamname = $href->{bamname};
                my $bamid = $href->{bamid};
                my $nwdid = $href->{expt_sampleid};
                my $newbamname = '';
                my $oldbamnamepath = '';
                my $newbamnamepath = '';
                #------------------------------------------------------------------
                #   Check that everything is ready for the change
                #------------------------------------------------------------------
                #   Do not play unless NWDID is as we expect and bam exists
                if ($nwdid !~ /^NWD/) {
                    print "  Skipping bam '$oldbamname' [$bamid] Malformed NWDID [$nwdid]\n";
                    next;
                }
                $oldbamnamepath = $opts{topdir} . "/$centername/$dirname/" . $oldbamname;
                if (! -r $oldbamnamepath) {
                    print "  Skipping bam '$oldbamname' [$bamid] File not found\n";
                    next;
                }
                if ($oldbamname !~ /^([^.]+)\.(.+)/) {
                    print "  Skipping bam '$oldbamname [$bamid] Unable to parse bam name\n";
                    next;
                }
                $newbamname = $nwdid . '.' . $2;
                $newbamnamepath = $opts{topdir} . "/$centername/$dirname/" . $newbamname;
                if ($oldbamname =~ /^NWD/) {
                    print "  Skipping bam '$oldbamname [$bamid] No rename necessary\n";
                    next;
                }

                #   Final hook - remapping must be done, otherwise I mess up Chris
                #if ((! defined($href->{datemapping})) ||
                #    ($href->{datemapping} eq '') ||
                #    ($href->{datemapping} < 10)) {
                #    print "  Skipping bam '$oldbamname [$bamid] Remapping not complete\n";
                #    next;
                #}

                #------------------------------------------------------------------
                #   Everything we need is in place, make changes
                #------------------------------------------------------------------
                #   Rename BAM file, bai and anything else remotely like it
                if (RunCMD("mv $oldbamnamepath $newbamnamepath")) { next; }
                foreach my $suff (qw(bai header md5)) {     # Ignore errors
                    if (! -f "$oldbamnamepath.$suff") { next; }
                    RunCMD("mv $oldbamnamepath.$suff $newbamnamepath.$suff");
                }
                
                #   Save original bamname once in database
                if ((! defined($href->{bamname_orig})) || $href->{bamname_orig} eq '') {
                    RunSQL("UPDATE $opts{bamfiles_table} SET bamname_orig='$oldbamname' WHERE bamid=$bamid");
                }

                #   Update bamname in database now that rename is done
                print "  Updated $oldbamname to $newbamname\n";
                RunSQL("UPDATE $opts{bamfiles_table} SET bamname='$newbamname' WHERE bamid=$bamid");

                #   Rename BACKUP data
                my $bakbamnamedir = $opts{backupdir} . "/$centername/$dirname";
                if (-f "$bakbamnamedir/$oldbamname") {
                    RunCMD("mv $bakbamnamedir/$oldbamname $bakbamnamedir/$newbamname");
                }

                #   Rename MD5 files so they do not get reprocessed
                if (! opendir(DIR, $qcdir)) {
                    print "  Unable to read directory '$qcdir'. Continuing\n";
                    next;
                }
                my @filelist = grep { /^$partname/ } readdir(DIR);
                closedir DIR;
                if (! @filelist) {
                    print "  SURPRISE: Found no files in '$qcdir' to rename! Continuing\n";
                    next;
                }
                foreach my $f (@filelist) {
                    my $oldf = "$qcdir/$f";
                    if ($f !~ /^([^.]+)\.(.+)/) {
                        print "  Unable to parse filename '$f' in $qcdir\n";
                        next;
                    }
                    RunCMD("mv $qcdir/$f $qcdir/$nwdid.$2");
                }




                
                $opts{topdir} . "/$centername/$dirname/" . $newbamname;
                
    if (! open(IN, "cat $d/*.md5 $d/*.MD5 *$d/Manifest.txt 2>/dev/null |")) {
         warn "No MD5 files found in directory '$d'. Maybe later, eh?\n";
         return 0;
    }


                #   Rename QPLOT output
                my $qcdir = $opts{resultsdir} . "/$centername/$dirname";
                my $partname = $oldbamname;
                $partname =~ s/.bam//;
                if (! opendir(DIR, $qcdir)) {
                    print "  Unable to read directory '$qcdir'. Continuing\n";
                    next;
                }
                my @filelist = grep { /^$partname/ } readdir(DIR);
                closedir DIR;
                if (! @filelist) {
                    print "  SURPRISE: Found no files in '$qcdir' to rename! Continuing\n";
                    next;
                }
                foreach my $f (@filelist) {
                    my $oldf = "$qcdir/$f";
                    if ($f !~ /^([^.]+)\.(.+)/) {
                        print "  Unable to parse filename '$f' in $qcdir\n";
                        next;
                    }
                    RunCMD("mv $qcdir/$f $qcdir/$nwdid.$2");
                }
                $opts{count}--;
                if ($opts{count} < 1) { print "Count exhausted, stopping early\n"; exit; }
            }
        }
    }
    exit;
}


die "Invalid request '$fcn'. Try '$Script --help'\n";

#==================================================================
# Subroutine:
#   RunCMD($cmd)
#
#   Execute command or just show what would be done
#==================================================================
sub RunCMD {
    my ($cmd) = @_;
    if ($opts{verbose}) {
        print "  CMD: $cmd\n";
        return(0);
    }
    if (system($cmd)) {
        print "  CMD FAILED:  $cmd\n";
        return(1);
    }
    return(0);
}

#==================================================================
# Subroutine:
#   RunSQL($sql)
#
#   Update SQL or just show what would be done
#==================================================================
sub RunSQL {
    my ($sql) = @_;
    if (! $sql) { return; }
    if ($opts{verbose}) { print "  DoSQL: $sql\n"; }
    else { DoSQL($sql); }
    return;
}

