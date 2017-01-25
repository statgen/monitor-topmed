#!/usr/bin/perl
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

my $topmedbin = '/usr/cluster/topmed/bin';
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
    topmedgce38push => "$topmedbin/topmed_gcepush.sh",
    topmedgce38pull => "$topmedbin/topmed_gcepull.sh",
    topmedgce38post => "$topmedbin/topmed_gcepost.sh",
    topmedxml    => "$topmedbin/topmed_xml.pl",
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    studies_table => 'studies',
    bamfiles_table => 'bamfiles',
    topdir => '/net/topmed/incoming/topmed',
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
    warn "$Script [options] arrive|verify|bai|qplot|cram|sexpt|sorig|sb37|spush|spull|spost|sb38\n" .
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
#   Push data to Google Cloud for processing
#--------------------------------------------------------------
if ($fcn eq 'spush') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles for Google Cloud
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams known at NCBI and that have not been sent
            my $sql = "SELECT * FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;

                #   Only send data if cram was done
                if ($href->{state_cram} != $COMPLETED) { next; }
                if ($opts{suberr} && $href->{state_gce38push} >= $FAILED) {
                    $href->{state_gce38push} = $REQUESTED;
                }
                if ($href->{state_gce38push} != $NOTSET &&
                    $href->{state_gce38push} != $REQUESTED) { next; }

                #   Send the CRAM to Google
                BatchSubmit("$opts{topmedgce38push} -submit $href->{bamid}");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Pull processed data from Google Cloud
#--------------------------------------------------------------
if ($fcn eq 'spull') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles for Google Cloud
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams known at NCBI and that have not been sent
            my $sql = "SELECT * FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;

                #   Only get data if remap was done
                if ($href->{state_gce38push} != $COMPLETED) { next; }
                if ($opts{suberr} && $href->{state_gce38pull} >= $FAILED) {
                    $href->{state_gce38pull} = $REQUESTED;
                }
                if ($href->{state_gce38pull} != $NOTSET &&
                    $href->{state_gce38pull} != $REQUESTED) { next; }

                #   Send the CRAM to Google
                BatchSubmit("$opts{topmedgce38pull} -submit $href->{bamid}");
            }
        }
    }
    ShowSummary($fcn);
    exit;
}

#--------------------------------------------------------------
#   Post process data we fetched from Google Cloud
#--------------------------------------------------------------
if ($fcn eq 'spost') {
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles for Google Cloud
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams known at NCBI and that have not been sent
            my $sql = "SELECT * FROM $opts{bamfiles_table} " .
                "WHERE runid='$runid'";
            if ($opts{datayear}) { $sql .= " AND datayear=$opts{datayear}"; }
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;

                #   Only post process data that we already fetched
                if ($href->{state_gce38pull} != $COMPLETED) { next; }
                if ($opts{suberr} && $href->{state_gce38post} >= $FAILED) {
                    $href->{state_gce38post} = $REQUESTED;
                }
                if ($href->{state_gce38post} != $NOTSET &&
                    $href->{state_gce38post} != $REQUESTED) { next; }

                #   Send the CRAM to Google
                BatchSubmit("$opts{topmedgce38post} -submit $href->{bamid}");
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
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_monitor.pl - Find runs that need some action

=head1 SYNOPSIS

  topmed_monitor.pl verify
  topmed_monitor.pl -run 20150604 qplot    # Select only samples from one run
  topmed_monitor.pl -center nygc verify    # Select only samples from a center
  topmed_monitor.pl -maxjobs 5 cram        # Only submit a few jobs
  topmed_monitor.pl -datayear 2 bai        # Send only year 2 samples

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

