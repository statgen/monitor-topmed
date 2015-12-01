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
my $NOTSET    = 0;            # Not set
my $REQUESTED = 1;            # Task requested
my $SUBMITTED = 2;            # Task submitted to be run
my $STARTED   = 3;            # Task started
my $DELIVERED = 19;           # Data delivered, but not confirmed
my $COMPLETED = 20;           # Task completed successfully
my $CANCELLED = 89;           # Task cancelled
my $FAILED    = 99;           # Task failed

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
    backupdir => '/net/topmed/working/backups/incoming/topmed',
    backupdir2 => '/net/topmed2/working/backups/incoming/topmed',
    resultsdir => '/incoming/qc.results',
    remappeddir => '/net/topmed/working/schelcj/results',
    remappeddir2 => '/net/topmed2/working/schelcj/results',
    count => 1000,
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s count=i
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] nwd|rename|crammd5|fixmapping\n" .
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
#   Calculate MD5 for all our CRAMs
#--------------------------------------------------------------
if ($fcn eq 'crammd5') {
    my $cons = '/net/topmed/working/topmed-output';
    my $submitcmd = "/usr/cluster/bin/sbatch --mem=2G -Q --workdir=$cons";
                    
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles that arrived
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            print "Doing $dirname\n";
            #   Get list of all bams
            my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE runid='$runid' AND piname!=NULL";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            my $b37fmissing = 0;
            my $bfmissing = 0;
            my $alreadydone = 0;
            my $submitted = 0;
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $bamid = $href->{bamid};

                #   Maybe this was already done?
                if ($href->{cramchecksum} && $href->{cramb37checksum}) { $alreadydone++; next; }

                #   Deal with b37 CRAM
                if (! $href->{cramb37checksum}) {
                    my $topmedhost = '';             # This is where the job should be run
                    #   Figure out where b37 cram lives
                    my $recalfile = $href->{expt_sampleid} . '.recal.cram';
                    my $part = "/$centername/" . $href->{piname} . '/' .
                        $href->{expt_sampleid} . '/bams/' . $recalfile;
                    my $f = $opts{remappeddir} . $part;
                    if (-f $f) { $topmedhost = 'topmed'; }
                    if (! $topmedhost) {
                        $f = $opts{remappeddir2} . $part;
                        if (-f $f) { $topmedhost = 'topmed2'; }
                        if (! $topmedhost) {
                            if ($opts{verbose}) { print "Remapped CRAM '$recalfile' not found\n"; }
                            $b37fmissing++;
                            next;
                        }
                    }
                    #   Submit jobs to calculate md5s
                    my $cmd = "$submitcmd -p $topmedhost-incoming --qos=$topmedhost-fixit " .
                        "-J $bamid-fixit --output=$cons/$bamid-fixit.$topmedhost.out " .
                        "/home/topmed/makemd5.sh $bamid $f cramb37checksum";
                    system($cmd) &&
                        print "  Failed: CMD=$cmd\n";
                    $submitted++;
                    if ($opts{verbose}) { print "  $recalfile to be calculated\n"; }
                }

                #   Deal with backup CRAM
                if (! $href->{cramchecksum}) {
                    my $topmedhost = '';             # This is where the job should be run
                    #   Figure out where backup cram lives
                    my $backupfile = $href->{expt_sampleid} . '.src.cram';
                    my $part = "/$centername/$dirname/$backupfile";
                    my $f = $opts{backupdir} . $part;
                    if (-f $f) { $topmedhost = 'topmed'; }
                    if (! $topmedhost) {
                        $f = $opts{backupdir2} . $part;
                        if (-f $f) { $topmedhost = 'topmed2'; }
                        if (! $topmedhost) {
                            if ($opts{verbose}) { print "Backup CRAM '$backupfile' not found\n"; }
                            $bfmissing++;
                            next;
                        }
                    }
                    #   Submit jobs to calculate md5s
                    my $cmd = "$submitcmd -p $topmedhost-incoming --qos=$topmedhost-fixit " .
                        "-J $bamid-fixit --output=$cons/$bamid-fixit.$topmedhost.out " .
                        "/home/topmed/makemd5.sh $bamid $f cramchecksum";
                    system($cmd) &&
                        print "  Failed: CMD=$cmd\n";
                    $submitted++;
                    if ($opts{verbose}) { print "  $backupfile to be calculated\n"; }
                }
                $opts{count}--;
                if ($opts{count} < 1) { last; }
            }
            print "  Submitted $submitted jobs for $centername/$dirname  " .
                "$alreadydone MD5s already done, " .
                "$b37fmissing B37 CRAMs missing, $bfmissing backup CRAMs missing\n";
            if ($opts{count} < 1) { print "Count exhausted, stopping early\n"; exit; }
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
                        die "$Script Failed to rename $f $f.old\n";
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
my $qcdir = 'this got broken somehow';
my $partname = 'this got broken somehow';

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

                #   Rename QPLOT output
                $qcdir = $opts{resultsdir} . "/$centername/$dirname";
                $partname = $oldbamname;
                $partname =~ s/.bam//;
                if (! opendir(DIR, $qcdir)) {
                    print "  Unable to read directory '$qcdir'. Continuing\n";
                    next;
                }
        @filelist = grep { /^$partname/ } readdir(DIR);
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

#--------------------------------------------------------------
#   Update state_b37 from datamapping until Chris' remapping finally completes
#   Remove comments to correct all those date* fields
#--------------------------------------------------------------
if ($fcn eq 'fixmapping') {
    my %old2new = (
#        datearrived => 'state_arrive',
#        datemd5ver => 'state_md5ver',
#        datebackup => 'state_backup',
#        datecram => 'state_cram',
#        datebai => 'state_bai',
#        dateqplot => 'state_qplot',
        datemapping => 'state_b37',
    );

    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles that arrived
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            print "Doing $dirname\n";
            #   Get list of all bams
            my $sql = "SELECT * FROM $opts{bamfiles_table} WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $bamid = $href->{bamid};
                #   Create UPDATE statements for each field to be changed
                my $sets = '';
                my $st = '';
                foreach my $oldcol (keys %old2new) {
                    my $newcol = $old2new{$oldcol};
                    if ($href->{$newcol} == $COMPLETED) { next; }
                    if (! defined($href->{$oldcol})) { $href->{$oldcol} = ' '; }
                    $st .= "newcol=$newcol [$href->{$newcol}] oldcol=$oldcol [$href->{$oldcol}]";
                    if (defined($href->{$oldcol}) && $href->{$oldcol} ne ' ') {
                        if ($href->{$oldcol} == -1) { $sets .= "$newcol=$FAILED,";    next; }
                        if ($href->{$oldcol} < 0)   { $sets .= "$newcol=$STARTED,";   next; }
                        if ($href->{$oldcol} == 0)  { $sets .= "$newcol=$REQUESTED,"; next; }
                        if ($href->{$oldcol} == 2)  { $sets .= "$newcol=$SUBMITTED,"; next; }
                        #if ($href->{$oldcol} == 3)  { $sets .= "$newcol=$DELIVERED,"; next; }
                        if ($href->{$oldcol} == 3)  { $sets .= "$newcol=$COMPLETED,"; next; }
                        if ($href->{$oldcol} == 1)  { $sets .= "$newcol=$CANCELLED,"; next; }
                        if ($href->{$oldcol} > 10)  { $sets .= "$newcol=$COMPLETED,"; next; }
                        print "Unable to handle $oldcol='$href->{$oldcol}'  bamid=$href->{bamid}\n";
                    }
                }
                if (! $sets) {
                    if ($opts{verbose}) { print "Nothing to do for BAMID=$href->{bamid}\n"; }
                    next;
                }
                chop($sets);
                $sql = "UPDATE $opts{bamfiles_table} SET $sets WHERE bamid=$href->{bamid}";
                #print "$st\n  $sql\n";
                my $sth2 = DoSQL($sql);
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

