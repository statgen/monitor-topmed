#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
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

use POSIX qw(strftime tmpnam);

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
my $FAILEDCHECKSUM = 98;      # Task failed, because checksum at NCBI bad
my $FAILED    = 99;           # Task failed

my $topmedbin = '/usr/cluster/monitor/bin';
our %opts = (
    topmedcmd => "$topmedbin/topmedcmd.pl",
    topmedarrive => "$topmedbin/topmed_arrive.sh",
    topmedverify => "$topmedbin/topmed_verify.sh",
    topmedcram   => "$topmedbin/topmed_cram.sh",
    topmedbai    => "$topmedbin/topmed_bai.sh",
    topmedqplot  => "$topmedbin/topmed_qplot.sh",
    topmedexpt  => "$topmedbin/topmed_ncbiexpt.sh",
    topmedncbiorig => "$topmedbin/topmed_ncbiorig.sh",
    topmedncbib37 => "$topmedbin/topmed_ncbib37.sh",
    topmedncbib38 => "$topmedbin/topmed_ncbib38.sh",
    topmedxml    => "$topmedbin/topmed_xml.pl",
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    studies_table => 'studies',
    bamfiles_table => 'bamfiles',
    topdir  => '/net/topmed/incoming/topmed',
    topdir2 => '/net/topmed2/incoming/topmed',
    topdir3 => '/net/topmed3/incoming/topmed',
    topdir4 => '/net/topmed4/incoming/topmed',
    topdir5 => '/net/topmed5/incoming/topmed',
    topdir6 => '/net/topmed6/incoming/topmed',
    backupdir => '/working/backups/incoming/topmed',
    resultsdir => '/incoming/qc.results',    
    dryrun => 0,
    verbose => 0,
    maxjobs => 100,
    jobcount => 0,              # Not actually an option, but stats
    jobsnotpermitted => 0,
    jobsfailedsubmission => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose topdir=s center=s runs=s maxjobs=i
    dryrun suberr datayear=i
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] arrive|verify|bai|qplot|cram|sexpt|sorig|sb37|sb38\n" .
        "Find runs which need some action and queue a request to do it.\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift(@ARGV);

my $dbh = DBConnect($opts{realm});

my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

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
            my $sql = "SELECT bamid,bamname,state_arrive FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" .
                    $href->{bamname};
                my @stats = stat($f);
                if (! @stats) { next; }
                #   See if we should mark this as arrived
                #   All BAMs must start with NWD. Only skip this if the BAM name
                #   is NWD and it is marked as completed
                if ($href->{bamname} =~ /^NWD/ && $href->{state_arrive} == $COMPLETED) { next; }
                #   If the mtime on the file is very recent, it might still be coming
                my $nowtime = time();
                if ((time() - $stats[9]) < 3600) { next; }
                #   Run the command
                BatchSubmit("$opts{topmedarrive} $href->{bamid} $f");
            }
        }
    }
    ShowSummary('BAMs arrived');
    exit;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have not been verified
#--------------------------------------------------------------
if ($fcn eq 'verify') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to be verified
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,state_arrive,state_md5ver,checksum,bamsize FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
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
                if ($href->{state_arrive} != $COMPLETED) { next; }
                if ($opts{suberr} && $href->{state_md5ver} >= $FAILEDCHECKSUM) { $href->{state_md5ver} = $REQUESTED; }
                if ($href->{state_md5ver} != $NOTSET && $href->{state_md5ver} != $REQUESTED) { next; }

                #   Run the command
                BatchSubmit("$opts{topmedverify} -submit $href->{bamid} $href->{checksum} $f");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have no BAI file
#--------------------------------------------------------------
if ($fcn eq 'bai') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there is a missing BAI file
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,state_md5ver,state_bai FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                #   Only do bai if file has been verified
                if ($href->{state_md5ver} != $COMPLETED) { next; }
                if ($opts{suberr} && $href->{state_bai} >= $FAILEDCHECKSUM) { $href->{state_bai} = $REQUESTED; }
                if ($href->{state_bai} != $NOTSET && $href->{state_bai} != $REQUESTED) { next; }
                #   Run the command
                BatchSubmit("$opts{topmedbai} -submit $href->{bamid} $f");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of BAMs that have not been converted to CRAMs
#--------------------------------------------------------------
if ($fcn eq 'cram') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to be backed up
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,state_md5ver,state_cram,state_bai FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                #   Only create cram if file has been verified
                if ($href->{state_md5ver} != $COMPLETED) { next; }
                if ($href->{state_bai} != $COMPLETED) { next; }
                if ($opts{suberr} && $href->{state_cram} >= $FAILEDCHECKSUM) { $href->{state_cram} = $REQUESTED; }
                if ($href->{state_cram} != $NOTSET && $href->{state_cram} != $REQUESTED) { next; }
                #   Run the command
                BatchSubmit("$opts{topmedcram} -submit $href->{bamid} $f");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Run QPLOT on BAMs
#--------------------------------------------------------------
if ($fcn eq 'qplot') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to run qplot on
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,state_bai,state_qplot FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                if (! -f $f) { next; }      # If BAM not there, do not submit
                #   Only do qplot if BAI finished
                if ($href->{state_bai} != $COMPLETED) { next; }
                #   Only do qplot if cram has been created
                ####if ($href->{state_cram} != $COMPLETED) { next; }
                $f = $opts{topdir} . "/$centername/$dirname/" . $href->{bamname};
                if (! -f $f) { next; }      # If BAM not there, do not submit

                if ($opts{suberr} && $href->{state_qplot} >= $FAILEDCHECKSUM) { $href->{state_qplot} = $REQUESTED; }
                if ($href->{state_qplot} != $NOTSET && $href->{state_qplot} != $REQUESTED) { next; }
                #   Run the command
                BatchSubmit("$opts{topmedqplot} -submit $href->{bamid} $f");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}
#--------------------------------------------------------------
#   Create experiment XML and sent it to NCBI for this NWDID
#--------------------------------------------------------------
if ($fcn eq 'sexpt') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to be delivered to NCBI
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams known at NCBI and that have not been sent
            my $sql = "SELECT * FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid' AND nwdid_known='Y'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) {
                if ($opts{runs}) { print "No BAMs from '$dirname' have nwdid_known='Y'\n"; }
                next;
            }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                
                #   Only do if this NWDID if all the local steps have completed
                if ($href->{state_qplot} != $COMPLETED) { next; }
                if ($href->{state_cram} != $COMPLETED) { next; }

                #   Check important fields for this BAM are possibly correct
                my $skip = '';
                foreach my $col (qw(expt_sampleid phs library_name nominal_length nominal_sdev base_coord)) {
                    if (exists($href->{$col}) && $href->{$col}) { next; }
                    $skip .= "$col ";
                }
                #   Do better checking than just is the field non-blank
                if ($href->{expt_sampleid} !~ /^NWD/) {
                    $skip .= ' Invalid_NWDID_' . $href->{expt_sampleid};
                }
                if ($skip) {
                    print "  BAM '$href->{bamname}' [$href->{bamid}] is ignored because of incomplete data for: $skip\n";
                    next;
                }

                if ($opts{suberr} && $href->{state_ncbiexpt} >= $FAILEDCHECKSUM) { $href->{state_ncbiexpt} = $REQUESTED; }
                if ($href->{state_ncbiexpt} != $NOTSET && $href->{state_ncbiexpt} != $REQUESTED) { next; }
                #   Tell NCBI about this NWDID experiment
                BatchSubmit("$opts{topmedexpt} -submit $href->{bamid}");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of secondary BAMs to be sent to NCBI
#--------------------------------------------------------------
if ($fcn eq 'sorig') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to be delivered to NCBI
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams known at NCBI and that have not been sent
            my $sql = "SELECT * FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid' AND nwdid_known='Y'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) {
                if ($opts{runs}) { print "No BAMs from '$dirname' have nwdid_known='Y'\n"; }
                next;
            }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;

                #   Only send the original BAM if the experiment was accepted at NCBI
                if ($href->{state_ncbiexpt} != $COMPLETED) { next; }

                #   Check important fields for this BAM are possibly correct
                my $skip = '';
                foreach my $col (qw(checksum expt_sampleid)) {
                    if (exists($href->{$col}) && $href->{$col}) { next; }
                    if ($opts{verbose}) { print "  No value for '$col'\n"; }
                    $skip .= "$col ";
                }
                if ($skip) {
                    print "  BAM '$href->{bamname}' [$href->{bamid}] is ignored because of incomplete data for: $skip\n";
                    next;
                }
                if ($opts{suberr} && $href->{state_ncbiorig} >= $FAILEDCHECKSUM) { $href->{state_ncbiorig} = $REQUESTED; }
                if ($href->{state_ncbiorig} != $NOTSET && $href->{state_ncbiorig} != $REQUESTED) { next; }
                #   Send the secondary BAM to NCBI
                BatchSubmit("$opts{topmedncbiorig} -submit $href->{bamid}");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of remapped primary BAMs to be sent to NCBI
#--------------------------------------------------------------
if ($fcn eq 'sb37') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to be delivered to NCBI
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams known at NCBI and that have not been sent
            my $sql = "SELECT * FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid' AND nwdid_known='Y'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;

                #   Only send the remapped b37 file for this year's data
                if ($href->{datayear} ne '1') {
                    if ($opts{verbose}) {
                        print "  BAM '$href->{bamname}' [$href->{bamid}] ignored because it is not year 1\n";
                    }
                    next;
                }

                #   Only send the remapped file if the experiment was accepted at NCBI
                if ($href->{state_ncbiexpt} != $COMPLETED) { next; }

                #   Check important fields for this BAM are possibly correct
                my $skip = '';
                foreach my $col (qw(cramname expt_sampleid)) {
                    if (exists($href->{$col}) && $href->{$col}) { next; }
                    if ($opts{verbose}) { print "  No value for '$col'\n"; }
                    $skip .= "$col ";
                }
                if ($skip) {
                    print "  BAM '$href->{bamname}' [$href->{bamid}] is ignored because of incomplete data for: $skip\n";
                    next;
                }
                if ($opts{suberr} && $href->{state_ncbib37} >= $FAILEDCHECKSUM) { $href->{state_ncbib37} = $REQUESTED; }
                if ($href->{state_ncbib37} != $NOTSET && $href->{state_ncbib37} != $REQUESTED) { next; }
                #   Send the remapped CRAM to NCBI
                BatchSubmit("$opts{topmedncbib37} -submit $href->{bamid}");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Get a list of remapped primary BAMs to be sent to NCBI
#--------------------------------------------------------------
if ($fcn eq 'sb38') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles to be delivered to NCBI
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams known at NCBI and that have not been sent
            my $sql = "SELECT * FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid' AND nwdid_known='Y'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;

                #   Only send the remapped b38 file for this year's data
                if ($href->{datayear} ne '2') { next; }

                #   Only send the remapped file if the experiment was accepted at NCBI
                if ($href->{state_ncbiexpt} != $COMPLETED) { next; }

                #   Check important fields for this BAM are possibly correct
                my $skip = '';
                foreach my $col (qw(cramname expt_sampleid)) {
                    if (exists($href->{$col}) && $href->{$col}) { next; }
                    if ($opts{verbose}) { print "  No value for '$col'\n"; }
                    $skip .= "$col ";
                }
                if ($skip) {
                    print "  BAM '$href->{bamname}' [$href->{bamid}] is ignored because of incomplete data for: $skip\n";
                    next;
                }
                if ($opts{suberr} && $href->{state_ncbib38} >= $FAILEDCHECKSUM) { $href->{state_ncbib38} = $REQUESTED; }
                if ($href->{state_ncbib38} != $NOTSET && $href->{state_ncbib38} != $REQUESTED) { next; }
                #   Send the remapped CRAM to NCBI
                BatchSubmit("$opts{topmedncbib38} -submit $href->{bamid}");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Check all the BAMs and see if the files match what the DB says
#--------------------------------------------------------------
if ($fcn eq 'redo-check') {
    $_ = "#" . '-' x 70;
    print "$_\n#    Found these inconsistencies between the database and file system\n" .
         "#    BEWARE OF FALSE ERRORS  when jobs are being processed\n$_\n";
    #   Get all the known centers in the database
    my $centersref = GetCenters();
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
    if ($opts{maxjobs} == 0 && $opts{verbose}) { print "Limit of jobs to be submitted has been reached\n"; }
    if ($opts{dryrun}) { print "dryrun => $cmd\n"; return; }
    my $rc = system("$cmd 2>&1");
    $rc = $rc >> 8;
    if ($rc == 0) {
        $opts{jobcount}++;
        if ($opts{verbose}) { print "submitted => $cmd\n"; }
        return;
    }
    $opts{maxjobs}++;                   # Submit failed, keep trying
    if ($rc == 4) { $opts{jobsnotpermitted}++; }
    else { $opts{jobsfailedsubmission}++; }
    return;
}

#==================================================================
# Subroutine:
#   ShowSummary - Print summary of jobs activity
#
# Arguments:
#   type - type of job, verify, backup etc
#==================================================================
sub ShowSummary {
    my ($type) = @_;

    my $s = '';
    if ($opts{jobcount})            { $s .= "$opts{jobcount} jobs submitted"; }
    if ($opts{jobsnotpermitted})    { $s .=  "  $opts{jobsnotpermitted} job submissions not permitted"; }
    if ($opts{jobsfailedsubmission}) { $s .=  "  $opts{jobsfailedsubmission} job submissions failed"; }
    if (! $s) { return; }
    print "$nowdate $type: $s\n";
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
    foreach my $pfx ('topdir','topdir2','topdir3', 'topdir4', 'topmed5', 'topmed6','topmed7', 'topmed8') {                    
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
  topmed_monitor.pl -run 20150604 verify   # Select only samples from one run
  topmed_monitor.pl -center nygc verify    # Select only samples from a center
  topmed_monitor.pl -maxjobs 5 sorig       # Only submit a few jobs
  topmed_monitor.pl -datayear 2 sexpt      # Send only year 2 samples

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

=item B<-datayear N>

Submit only jobs for samples in a specific year.

=item B<-dryrun>

Do not submit any jobs, just show the command to be executed.

=item B<-help>

Generates this output.

=item B<-maxjobs N>

Do not submit more than N jobs for this invocation.
The default for B<-maxjobs> is B<100>.

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

=item B<arrive | verify | bai | qplot | cram | sexpt | sorig | sb37 | sb38>

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

