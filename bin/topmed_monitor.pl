#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name: topmed_monitor.pl
#
# Description:
#   Use this program to automatically request actions on data
#   from NHLBI TopMed centers.
#   This program is expected to run as a crontab job.
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
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../lib/perl5";
use My_DB;
use TopMed_Get;
use Getopt::Long;

use POSIX qw(strftime);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my $topmedbin = '/usr/cluster/monitor/bin';
our %opts = (
    topmedcmd => "$topmedbin/topmedcmd.pl",
    topmedarrive => "$topmedbin/topmed_arrive.sh",
    topmedverify => "$topmedbin/topmed_verify.sh",
    topmedbackup => "$topmedbin/topmed_backup.sh",
    topmedcram   => "$topmedbin/topmed_cram.sh",
    topmedbai    => "$topmedbin/topmed_bai.sh",
    topmedqplot  => "$topmedbin/topmed_qplot.sh",
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    studies_table => 'studies',
    bamfiles_table => 'bamfiles',
    topdir  => '/net/topmed/incoming/topmed',
    topdir2 => '/net/topmed2/incoming/topmed',
    topdir3 => '/net/topmed3/incoming/topmed',
    backupdir => '/working/backups/incoming/topmed',
    resultsdir => '/incoming/qc.results',    
    dryrun => 0,
    verbose => 0,
    maxjobs => 100,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose topdir=s center=s runs=s maxjobs=i
    memory=s partition=s qos=s dryrun suberr
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] arrive|verify|backup|bai|qplot|cram|cp2ncbi\n" .
        "Find runs which need some action and queue a request to do it.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift(@ARGV);

my $dbh = DBConnect($opts{realm});

my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

#   Set environment variables for shell scripts that do -submit
if ($opts{memory})    { $ENV{TOPMED_MEMORY} = $opts{memory}; }
if ($opts{partition}) { $ENV{TOPMED_PARTITION} = $opts{partition}; }
if ($opts{qos})       { $ENV{TOPMED_QOS} = $opts{qos}; }

#--------------------------------------------------------------
#   dateXXXX columns are either a time when something was done
#   or a flag indicating varying states. For example with verify:
#   datemd5ver not defined   - nothing every happened
#   datemd5ver > 10    verify successfully
#   datemd5ver < 0     started to verify
#   datemd5ver = -1    verify failed
#   datemd5ver = 0     verify requested
#   datemd5ver = 2     verify submitted to be done (not done for arrived)
#   datemd5ver = 1     verify cancelled
#--------------------------------------------------------------

#--------------------------------------------------------------
#   Get a list of BAMs we have not noticed they arrived yet
#--------------------------------------------------------------
if ($fcn eq 'arrive') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles that arrived
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,datearrived FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" .
                    $href->{bamname};
                my @stats = stat($f);
                if (! @stats) { next; }
                #   See if we should mark this as arrived. Few states possible
                if (defined($href->{datearrived}) && ($href->{datearrived} ne '')) {
                    if ($href->{datearrived} > 10) { next; }     # Already arrived
                }
                #   If the mtime on the file is very recent, it might still be coming
                my $nowtime = time();
                if ((time() - $stats[9]) < 3600) { next; }
                #   Run the command
                BatchSubmit("$opts{topmedarrive} $href->{bamid} $f");
            }
        }
    }
    if ($opts{jobcount}) { print "$nowdate  $opts{jobcount} BAMs arrived\n"; }
    exit;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have not been verified
#--------------------------------------------------------------
if ($fcn eq 'verify') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    $opts{jobcount} = 0;
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to be verified
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,datearrived,datemd5ver,checksum,bamsize FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                #   If needed, update bamsize
                if (! $href->{bamsize}) {
                    my @stats = stat($f);
                    if ($stats[7]) {
                        $sql = "UPDATE $opts{bamfiles_table} SET bamsize=$stats[7] " .
                            "WHERE bamid=$href->{bamid}";
                        my $sth2 = DoSQL($sql);
                    }
                }
                #   Only do verify if file has arrived and ready to be run
                if (! defined($href->{datearrived})) { next; }
                if ($href->{datearrived} eq '') { next; }    # Not here
                if ($href->{datearrived} =~ /\D/) { next; }  # Not numeric
                if ($href->{datearrived} < 10) { next; }     # Not arrived, avoid verify
                if (defined($href->{datemd5ver}) && ($href->{datemd5ver} ne '')) {
                    if ($opts{suberr} && $href->{datemd5ver} == -1) { $href->{datemd5ver} = 0; }
                    if ($href->{datemd5ver} =~ /\D/) { next; }  # Not numeric
                    if ($href->{datemd5ver} < 0) { next; }  # Already started
                    if ($href->{datemd5ver} == 2) { next; } # Already queued to be run
                    if ($href->{datemd5ver} > 10) { next; } # Already done
                }
                #   Run the command
                BatchSubmit("$opts{topmedverify} -submit $href->{bamid} $href->{checksum} $f");
            }
        }
    }
    if ($opts{jobcount}) { print "$nowdate  Submitted $opts{jobcount} jobs\n"; }
    exit;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have no BAI file
#--------------------------------------------------------------
if ($fcn eq 'bai') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    $opts{jobcount} = 0;
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there is a missing BAI file
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,datemd5ver,datebai FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                #   Only do bai if file has arrived
                if (! defined($href->{datemd5ver})) { next; }
                if ($href->{datemd5ver} eq '') { next; }    # Not here
                if ($href->{datemd5ver} =~ /\D/) { next; }  # Not numeric
                if ($href->{datemd5ver} < 10) { next; }     # Not arrived, avoid bai
                if (defined($href->{datebai}) && ($href->{datebai} ne '')) {
                    if ($opts{suberr} && $href->{datebai} == -1) { $href->{datebai} = 0; }
                    if ($href->{datebai} =~ /\D/) { next; }  # Not numeric
                    if ($href->{datebai} < 0) { next; }  # Already started
                    if ($href->{datebai} == 2) { next; } # Already queued to be run
                    if ($href->{datebai} > 10) { next; } # Already done
                }
                #   Run the command
                BatchSubmit("$opts{topmedbai} -submit $href->{bamid} $f");
            }
        }
    }
    if ($opts{jobcount}) { print "$nowdate  Submitted $opts{jobcount} jobs\n"; }
    exit;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have not been backed up
#--------------------------------------------------------------
if ($fcn eq 'backup') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    $opts{jobcount} = 0;
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to be backed up
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,datearrived,datebackup FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                #   Only do backup if file has arrived and ready to be run
                if (! defined($href->{datearrived})) { next; }
                if ($href->{datearrived} eq '') { next; }    # Not here
                if ($href->{datearrived} =~ /\D/) { next; }  # Not numeric
                if ($href->{datearrived} < 10) { next; }     # Not arrived, avoid backup
                if (defined($href->{datebackup}) && ($href->{datebackup} ne '')) {
                    if ($opts{suberr} && $href->{datebackup} == -1) { $href->{datebackup} = 0; }
                    if ($href->{datebackup} =~ /\D/) { next; }  # Not numeric
                    if ($href->{datebackup} < 0) { next; }  # Already started
                    if ($href->{datebackup} == 2) { next; } # Already queued to be run
                    if ($href->{datebackup} > 10) { next; } # Already done
                }
                #   Run the command
                BatchSubmit("$opts{topmedbackup} -submit $href->{bamid} $f");
            }
        }
    }
    if ($opts{jobcount}) { print "$nowdate  Submitted $opts{jobcount} jobs\n"; }
    exit;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have not been converted to CRAMs
#--------------------------------------------------------------
if ($fcn eq 'cram') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    $opts{jobcount} = 0;
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to be backed up
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,datearrived,datebackup,datecram FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                #   Only do cram if file has arrived and been backed up
                if (! defined($href->{datearrived})) { next; }
                if ($href->{datearrived} eq '') { next; }    # Not here
                if ($href->{datearrived} =~ /\D/) { next; }  # Not numeric
                if ($href->{datearrived} < 10) { next; }     # Not arrived, avoid backup
                if (defined($href->{datebackup}) && ($href->{datebackup} ne '')) {
                    if ($href->{datebackup} =~ /\D/) { next; }  # Not numeric
                    if ($href->{datebackup} <= 10) { next; } # Not done
                }
                if (defined($href->{datecram}) && ($href->{datecram} ne '')) {
                    if ($opts{suberr} && $href->{datecram} == -1) { $href->{datecram} = 0; }
                    if ($href->{datecram} =~ /\D/) { next; }  # Not numeric
                    if ($href->{datecram} < 0) { next; }  # Already started
                    if ($href->{datecram} == 2) { next; } # Already queued to be run
                    if ($href->{datecram} > 10) { next; } # Already done
                }
                #   Run the command
                BatchSubmit("$opts{topmedcram} -submit $href->{bamid} $f");
            }
        }
    }
    if ($opts{jobcount}) { print "$nowdate  Submitted $opts{jobcount} jobs\n"; }
    exit;
}

#--------------------------------------------------------------
#   Run QPLOT on BAMs
#--------------------------------------------------------------
if ($fcn eq 'qplot') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    $opts{jobcount} = 0;
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to run qplot on
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,datebai,dateqplot FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                if (! -f $f) { next; }      # If BAM not there, do not submit
                #   Only do qplot if BAI finished
                if (! defined($href->{datebai})) { next; }
                if ($href->{datebai} eq '') { next; }
                if ($href->{datebai} < 10) { next; }
                if (defined($href->{dateqplot}) && ($href->{dateqplot} ne '')) {
                    if ($opts{suberr} && $href->{dateqplot} == -1) { $href->{dateqplot} = 0; }
                    if ($href->{dateqplot} =~ /\D/) { next; }  # Not numeric
                    if ($href->{dateqplot} < 0) { next; }  # Already started
                    if ($href->{dateqplot} == 2) { next; } # Already queued to be run
                    if ($href->{dateqplot} > 10) { next; } # Already done
                }
                #   Run the command
                BatchSubmit("$opts{topmedqplot} -submit $href->{bamid} $f");
            }
        }
    }
    if ($opts{jobcount}) { print "$nowdate  Submitted $opts{jobcount} jobs\n"; }
    exit;
}

#--------------------------------------------------------------
#   Check all the BAMs and see if the files match what the DB says
#--------------------------------------------------------------
if ($fcn eq 'check') {
    $_ = "#" . '-' x 70;
    print "$_\n#    Found these inconsistencies between the database and file system\n" .
         "#    BEWARE OF FALSE ERRORS  when jobs are being processed\n$_\n";
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    $opts{jobcount} = 0;
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        print "Checking center=$centername\n";
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to run qplot on
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT * FROM " . $opts{bamfiles_table} .
                " WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();

            #   Check if the directory even exists anymore
            my $bd = $opts{topdir} . "/$centername/$dirname";
            if (! -d $bd) {
                print "$runid $dirname no longer exists\n";
                if ($rowsofdata) {      # There is stuff in database
                    print "We have bamfiles in database for $runid ($dirname)\n";
                    print "DELETE FROM $opts{bamfiles_table} WHERE runid=$runid;\n";
                }
                print "DELETE FROM $opts{runs_table} WHERE runid=$runid;\n";
                next;
            }

            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                #   datearrived means the BAM exists
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                if (defined($href->{datearrived}) &&
                    ($href->{datearrived} =~ /^d+/) &&
                    ($href->{datearrived} > 10) &&
                    (! -f $f)) {
                    print "$href->{bamid} $href->{bamname} arrived: File '$f' does not exist\n";
                    next;
                }
                #   Watch for files that are world writable
                my @stats = stat($f);
                my $mode = sprintf('%04o', $stats[2] & 07777);  # Normal looking UNIX permissions
                my $worldpermissions = substr($mode,3,1);
                if ($worldpermissions eq '6' || $worldpermissions eq '7') {
                    print "chmod 0640 $f    # Was $mode\n"; 
                }
                #   Make sure file size is correct
                if (! $href->{bamsize}) {
                    print "$href->{bamid} bamsize not set\n";
                    if ($stats[7]) {
                        print "UPDATE $opts{bamfiles_table} SET bamsize=$stats[7] WHERE bamid=$href->{bamid};\n";
                    }
                }

                #   See if the nwdid file exists
                my $bf = basename($href->{bamname}, '.bam');
                $bf = $opts{resultsdir} . "/$centername/$dirname/$bf.nwdid";
                if (! -f $bf) {
                    print "$href->{bamid} $href->{bamname} arrive: File '$bf' does not exist\n";
                }

                #   dateverify cannot be checked in any easy manner

                #   datebackup means the backup file should exist with the nwdid
                if (! defined($href->{nwdid})) { $href->{nwdid} = ''; }
                if ($href->{nwdid} !~ /^NWD/) {
                    print "$href->{bamid} $href->{bamname} backup: NWDID is not set\n";
                }
                if (defined($href->{datebackup}) &&
                    ($href->{datebackup} =~ /^\d+/) &&
                    ($href->{datebackup} > 10)) {
                    my $bamf = "/$centername/$dirname/" . $href->{bamname};
                    my $bamprefix = FindPrefix($bamf);
                    if (! $bamprefix) {
                        print "$href->{bamid} $href->{bamname} backup: File '$bamf' does not exist\n";
                    }
                    else {
                        $bf = $opts{backupdir} . "/$centername/$dirname/" .
                            $href->{nwdid} . '.src.cram';
                        my $prefix = FindPrefix($bf);
                        if (! $prefix) {
                            print "$href->{bamid} $href->{bamname} backup: Unable to find BAM file '$bf'\n";
                        }
                        else {                          # Poor way to figure out where backup should be
                            my $backdir = '';
                            if ($prefix eq $opts{topdir})  { $backdir = $opts{topdir2}; }
                            if ($prefix eq $opts{topdir2}) { $backdir = $opts{topdir}; }
                            if (! $backdir) {
                                print "$href->{bamid} $href->{bamname} backup: Unable to figure out " .
                                    "where BAM backup file '$bf' should be\n";
                            }
                            else {
                                $bf = $backdir . '/' . $bf;     # Real backup path
                                if (! -f $bf) {
                                    print "$href->{bamid} $href->{bamname} backup: " .
                                        "Unable to find backup file '$bf'\n";
                                }
                            }
                        }
                    }
                }
                #   datebai means the BAI file exists
                $bf = $f . '.bai';
                if (defined($href->{datebai}) &&
                    ($href->{datebai} =~ /^d+/) &&
                    ($href->{datebai} > 10) &&
                    (! -f $bf)) {
                    print "$href->{bamid} $href->{bamname} bai: File '$bf' does not exist\n";
                }
                #   If the bai file exists, then datebai should be a date
                if (! defined($href->{datebai})) { $href->{datebai} = 0; }
                if (-f $bf && $href->{datebai} < 10) {
                    print "$href->{datebai} $href->{bamname} bai: ".
                        "File '$bf' exists, but state=$href->{datebai}\n";
                    print "UPDATE $opts{bamfiles_table} SET datebai=" .
                        (-($href->{datebai})) . " WHERE bamid=$href->{bamid};\n";
                }
                #   dateqplot means the qplot and verifybamid results exist
                $bf = basename($href->{bamname}, '.bam');
                foreach my $suff ('qp.R', 'vb.selfSM') {
                    my $bff = $opts{resultsdir} . "/$centername/$dirname/$bf.$suff";
                    if (defined($href->{dateqplot}) &&
                        ($href->{dateqplot} =~ /^d+/) &&
                        ($href->{dateqplot} > 10) &&
                        (! -f $bff)) {
                        print "$href->{bamid} $href->{bamname} bai: File '$bff' does not exist\n";
                    }
                    #   If the qplot files exist, then dateqplot should be a date
                    if (-f $bff && $href->{dateqplot} < 10) {
                        print "$href->{bamid} $href->{bamname} qplot: " .
                            "File '$bff' exists, but state=$href->{dateqplot}\n";
                        print "UPDATE $opts{bamfiles_table} SET dateqplot=" . (-($href->{dateqplot})) .
                            " WHERE bamid=$href->{bamid};\n";
                    }
                    #   If the qplot files are zero length, then qplot failed
                    if (-z $bff) {
                        print "$href->{bamid} $href->{bamname} qplot: " .
                            "File '$bff' is zero length\n";
                        print "UPDATE $opts{bamfiles_table} SET dateqplot=-1" .
                            " WHERE bamid=$href->{bamid};\n";
                    }
                }
           }
        }
    }
    exit;
}

die "Invalid request '$fcn'. Try '$Script --help'\n";


#==================================================================
# Subroutine:
#   BatchSubmit - Run a command which submits the command to batch
#
# Arguments:
#   cmd - command
#==================================================================
sub BatchSubmit {
    my ($cmd) = @_;
    $opts{maxjobs}--;
    if ($opts{maxjobs} < 0) { return; }
    if ($opts{maxjobs} == 0) { print "Maximum limit of jobs that can be submitted has been reached\n"; }
    if ($opts{dryrun}) { print "dryrun => $cmd\n"; return; }
    system("$cmd 2>&1") || $opts{jobcount}++;
} 

#==================================================================
# Subroutine:
#   FindPrefix - Return topmed prefix for a file
#
# Arguments:
#   f = partial path to file in $opts{topmed*} directories
#
# Returns:
#   prefix for file
#==================================================================
sub FindPrefix {
    my ($f) = @_;
    foreach my $pfx ('topdir','topdir2','topdir3') {                    
        if (-f "$opts{$pfx}/$f") { return $opts{$pfx}; }
    }
    return '';
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_monitor.pl - Find runs that need some action

=head1 SYNOPSIS

  topmed_monitor.pl verify

=head1 DESCRIPTION

Use program as a crontab job to find runs that have not had certain
steps completed. When one is found a request is queued to
complete the step.

This process will be most successful if this program is run
after one might expect the previous step has completed.
For instance for verifying a run, run this in the early morning when you
might expect all the data has arrived.

=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-dryrun>

Do not submit any jobs, just show the command to be executed.

=item B<-help>

Generates this output.

=item B<-maxjobs N>

Do not submit more than N jobs for this invocation.
The default for B<-maxjobs> is B<100>.

=item B<-memory nG>

Force the sbatch --mem setting when submitting a job.

=item B<-partition name>

Force the sbatch --partition setting when submitting a job.

=item B<-qos name>

Force the sbatch --qos setting when submitting a job.

=item B<-realm NAME>

Specifies the database realm to read data from. This defaults to B<topmed>;

=item B<-runs NAME[,NAME,...]>

Specifies a specific set of runs on which to run the action,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-suberr>

Submit the job if the state is B<error>. Normally tasks in this state
are not submitted to be run.

=item B<-topdir PATH>

Specifies the path to where the tree of BAMs exists. This defaults to  B</incoming/topmed>;

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<arrive | verify | backup | bai | qplot | cram | cp2ncbi>

Directs this program to look for runs that have not been through the process name
you provided and to queue a request they be verified.

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

