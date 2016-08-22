#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
###################################################################
#
# Name: topmedcmd.pl
#
# Description:
#   Use this program to update the NHLBI database.
#   This program is tightly coupled with /monitor/topmed/index.php
#
# ChangeLog:
#   $Log: topmedcmd.pl,v $
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
use TopMed_Get;
use My_DB;
use Getopt::Long;
use Cwd qw(realpath abs_path);
use YAML;

my $NOTSET = 0;                     # Not set
my $REQUESTED = 1;                  # Task requested
my $SUBMITTED = 2;                  # Task submitted to be run
my $STARTED   = 3;                  # Task started
my $DELIVERED = 19;                 # Data delivered, but not confirmed
my $COMPLETED = 20;                 # Task completed successfully
my $CANCELLED = 89;                 # Task cancelled
my $FAILED    = 99;                 # Task failed

my %VALIDVERBS = (                  # Valid verbs to database colum
    arrived     => 'state_arrive',
    md5verified => 'state_md5ver',
    backedup    => 'state_backup',  
    baid        => 'state_bai',  
    qploted     => 'state_qplot',
    cramed      => 'state_cram',   
    sentexpt    => 'state_ncbiexpt',
    mapped37    => 'state_b37',   
    mapped38    => 'state_b38',     
    sentorig    => 'state_ncbiorig',
    sentb37     => 'state_ncbib37',
    sentb38     => 'state_ncbib38', 
);
my %VALIDSTATUS = (                 # Valid status for the verbs
   notset    => $NOTSET,
   requested => $REQUESTED,
   submitted => $SUBMITTED,
   started   => $STARTED,
   delivered => $DELIVERED,
   completed => $COMPLETED,
   cancelled => $CANCELLED,
   failed    =>  $FAILED, 
);
my %VALIDOPS = (                    # Used for permit
    all => 1,
    verify => 1,
    backup => 1,
    bai => 1,
    qplot => 1,
    cram => 1,
    nwdid => 1,
    sb37 => 1,
    sb38 => 1,
);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    centers_table => 'centers',
    permissions_table => 'permissions',
    runs_table => 'runs',
    netdir => '/net/topmed',
    qcresultsdir => 'incoming/qc.results',
    incomingdir => 'incoming/topmed',
    backupsdir => 'working/backups/incoming/topmed',
    consoledir => 'working/topmed-output',
    resultsdir => 'working/schelcj/results',
    ascpcmd => "/usr/cluster/bin/ascp -i ".
        "/net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri -l 800M -k 1",
    ascpdest => 'asp-um-sph@gap-submit.ncbi.nlm.nih.gov:protected',
    squeuecmd => "/usr/cluster/bin/squeue -o '%.10i %.18P %.15q %.18j %.8u %.2t %.10M %.6D %R'",
    maxlongjobs => 10,
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose maxlongjobs=i persist with-bamid|bamid
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options] [-persist]";
    my $verbs = join(',', sort keys %VALIDVERBS);
    my $requests = join(',', sort keys %VALIDSTATUS);
    warn "$m mark bamid|nwdid $verbs $requests\n" .
        "  or\n" .
        "$m unmark bamid|nwdid [same list as mark]\n" .
        "  or\n" .
        "$m set bamid|nwdid|dirname colname value\n" .
        "  or\n" .
        "$m show bamid|nwdid colname|run|center|yaml\n" .
        "  or\n" .
        "$m [-with-bamid] showrun runname\n" .
        "  or\n" .
        "$m export\n" .
        "  or\n" .
        "$m send2ncbi files\n" .
        "  or\n" .
        "$m squeue\n" .
        "  or\n" .
        "$m wherepath|where bamid|nwdid bam|backup|qcresults|console|b37|b38\n" .
        "  or\n" .
        "$m wherehost bamid|nwdid bam|backup|qcresults\n" .
        "  or\n" .
        "$m wherefile bamid|nwdid bam|backup|qcresults\n" .
        "  or\n" .
        "$m whatbamid bamname\n" .
        "  or\n" .
        "$m whatnwdid NWDnnnnn\n" .
        "  or\n" .
        "$m permit add operation center run\n" .
        "$m permit remove permitid\n" .
        "$m permit test operation bamid/n" .
        "\nUpdate the topmed database\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

if ($opts{persist}) { DBConnect($opts{realm}); }    # Open database
else  { PersistDBConnect($opts{realm}); }

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'mark')      { Mark(@ARGV); exit; }
if ($fcn eq 'unmark')    { UnMark(@ARGV); exit; }
if ($fcn eq 'set')       { Set(@ARGV); exit; }
if ($fcn eq 'show')      { Show(@ARGV); exit; }
if ($fcn eq 'showrun')   { ShowRun(@ARGV); exit; }
if ($fcn eq 'export')    { Export(@ARGV); exit; }
if ($fcn eq 'send2ncbi') { Send2NCBI(@ARGV); exit; }
if ($fcn eq 'squeue')    { SQueue(@ARGV); exit; }
if ($fcn eq 'where' || $fcn eq 'wherepath') { WherePath(@ARGV); exit; }
if ($fcn eq 'whathost')  { WhatHost(@ARGV); exit; }
if ($fcn eq 'wherefile') { WhereFile(@ARGV); exit; }
if ($fcn eq 'whatbamid') { WhatBAMID(@ARGV); exit; }
if ($fcn eq 'whatnwdid') { WhatNWDID(@ARGV); exit; }
if ($fcn eq 'permit')    { Permit(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   Mark($bamid, $op, $state)
#
#   Set states in the bamfiles database
#==================================================================
sub Mark {
    my ($bamid, $op, $state) = @_;
    $bamid = GetBamid($bamid);
    if ((! exists($VALIDVERBS{$op})) || (! exists($VALIDSTATUS{$state}))) {
        die "$Script - Invalid 'mark' syntax. Try '$Script -help'\n";
    }

    #   Make sure this is a bam we know
    my $sth = ExecSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    #   Set state for the verb
    my $col = $VALIDVERBS{$op};
    my $done = 0;
    if ($state eq 'requested') {
        ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$REQUESTED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'completed') {
        ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$COMPLETED WHERE bamid=$bamid");
        if ($col eq 'state_arrive') {       # hack for Chris until new code in place
            ExecSQL("UPDATE $opts{bamfiles_table} SET datearrived='" . time() . "' WHERE bamid=$bamid");
        }
        if ($col eq 'state_md5ver') {       # hack for Chris until new code in place
            ExecSQL("UPDATE $opts{bamfiles_table} SET datemd5ver='" . time() . "' WHERE bamid=$bamid");
        }
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            ExecSQL("UPDATE $opts{bamfiles_table} SET datemapping='" . time() . "' WHERE bamid=$bamid");
        }
        $done++;
    }
    if ($state eq 'delivered') {
        ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$DELIVERED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'started') {
        ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$STARTED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'failed') {
        ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$FAILED WHERE bamid=$bamid");
        if ($col eq 'state_arrive') {       # hack for Chris until new code in place
            ExecSQL("UPDATE $opts{bamfiles_table} SET datearrived='-1' WHERE bamid=$bamid");
        }
        if ($col eq 'state_md5ver') {       # hack for Chris until new code in place
            ExecSQL("UPDATE $opts{bamfiles_table} SET datemd5ver='-1' WHERE bamid=$bamid");
        }
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            ExecSQL("UPDATE $opts{bamfiles_table} SET datemapping='-1' WHERE bamid=$bamid");
        }
        $done++;
    }
    if ($state eq 'cancelled') {
        ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$CANCELLED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'submitted') {
        ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$SUBMITTED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'notset') {
        ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$NOTSET WHERE bamid=$bamid");
        $done++;
    }
    if ($done) {
        if ($opts{verbose}) { print "$Script  'mark $bamid $op $state'  successful\n"; } 
    }
    else { die "$Script - Invalid state '$state' for '$op'. Try '$Script -help'\n"; }
}

#==================================================================
# Subroutine:
#   UnMark($dirname, $op)
#
#   Reset state in the bamfiles database
#==================================================================
sub UnMark {
    my ($bamid, $op) = @_;
    $bamid = GetBamid($bamid);
    if (! exists($VALIDVERBS{$op})) {
        die "$Script - Invalid 'unmark' syntax. Try '$Script -help'\n";
    }

    #   Make sure this is a bam we know
    my $sth = ExecSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    my $col = $VALIDVERBS{$op};
    ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$NOTSET WHERE bamid=$bamid");
    if ($opts{verbose}) { print "$Script  'unmark $bamid $op'  successful\n"; }
}

#==================================================================
# Subroutine:
#   Export()
#
#   Generate a CSV file of possibly interesting data to STDOUT
#   This is incomplete and will need more attention, but works for Chris now
#==================================================================
sub Export {
    #my ($center, $run) = @_;

    #   Generate header for CSV file
    my @cols = qw(bamname studyname piname expt_sampleid state_b37 state_b38 datayear);
    my $s = 'CENTER,DIRNAME,';
    foreach (@cols) { $s .= uc($_) . ','; }
    print $s . "FULLPATH\n";
    
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles that arrived
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT * FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            my $sth = ExecSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = "$opts{netdir}/$opts{incomingdir}/$centername/$dirname/" .
                    $href->{bamname};
                #   See if this has been verified
                if ($href->{state_md5ver} != $COMPLETED) { next; }
                #   Show data for this BAM
                $s = "$centername,$dirname,";
                foreach (@cols) {
                    if (! defined($href->{$_})) { $s .= ','; }
                    else { $s .= $href->{$_} . ','; }
                 }
                $s .= $f . ' ';
                chop($s);
                print $s . "\n";
            }
        }
    }
}

#==================================================================
# Subroutine:
#   WhatBAMID($bamname)
#
#   Print bamid for a bamfiles entry with this bamname
#==================================================================
sub WhatBAMID {
    my ($bamname) = @_;

    if ($bamname =~ /\/(\S+)$/) { $bamname = $1; }
    my $sth = ExecSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamname='$bamname'");
    if (! $sth) { exit 1; }
    my $href = $sth->fetchrow_hashref;
    if (defined($href->{bamid})) { print $href->{bamid} . "\n"; }
}

#==================================================================
# Subroutine:
#   WhatNWDID($nwdid)
#
#   Print interesting details about an NWDID
#==================================================================
sub WhatNWDID {
    my ($nwdid) = @_;

    if ($nwdid =~ /^\d+/) {             # If bamid, get NWDID
        my $sth = ExecSQL("SELECT expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$nwdid");
        if ($sth) {
            my $href = $sth->fetchrow_hashref;
            $nwdid = $href->{expt_sampleid};
        }
    }
    else {                              # Extrace NWD from whatever was provided
        if ($nwdid =~ /(nwd\d+)/i) { $nwdid = uc($1); }
    }

    #   Reconstruct partial path to BAM
    my $sth = ExecSQL("SELECT runid,bamid FROM $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - NWDID '$nwdid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $bamid = $href->{bamid};
    $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$href->{runid}");
    $href = $sth->fetchrow_hashref;
    my $run = $href->{dirname};
    $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
    $href = $sth->fetchrow_hashref;
    my $center = uc($href->{centername});
    print "$nwdid/$bamid can be found in run '$run' for center $center\n";
}

#==================================================================
# Subroutine:
#   Send2NCBI($files)
#
#   Send files to NCBI
#==================================================================
sub Send2NCBI {
    my ($files) = @_;

    my $cmd = "$opts{ascpcmd} $files $opts{ascpdest}";
    my $errs = "/tmp/Send2NCBI.$$";
    foreach (1 .. 10) {             # Keep trying for auth errors
        my $rc = system("$cmd 2> $errs");
        if (! $rc) {
            unlink($errs);
            exit($rc);
        }
        #   Send failed. If this is an authentication error, wait and retry
        my $sleep = 60;
        if (open(IN, $errs)) {
            while (<IN>) {
                if (/Session Stop/) {
                    print "$Script - $_";
                    last;
                }     
                if (/ascp:/) {
                    print "$Script - $_";
                    if (/authenticate/) {
                        $sleep = 15;
                        last;
                    }
                }
            }
            close(IN);
            if ($sleep) { print "$Script - WARN: ASCP error, wait and retry\n"; }
        }
        sleep($sleep);
    }
    print "$Script - FATAL: Excessive ASCP fatal errors\n";
    unlink($errs);
    exit(3);
}

#==================================================================
# Subroutine:
#   SQueue()
#
#   Print summary of topmed-related jobs
#==================================================================
sub SQueue {
    #my ($bamid, $set) = @_;
    my $cmd = $opts{squeuecmd} . '| grep topmed';
    my %partitions = ();            # Unique partition names
    my %qos = ();                   # Unique QOS names
    my %running = ();
    my %mosixrunning = ();          # Counts of jobs by nomosix users
    my %queued = ();
    my $nottopmed = 0;              # Count of jobs not from my user
    foreach my $l (split('\n', `$cmd`)) {
        my @c = split(' ', $l);
        if ($c[1] eq 'nomosix' && $c[8] =~ /topmed/) {    # nomosix on topmed node
            $mosixrunning{$c[4]}++;    # Count of running jobs for this user
            next;
        }
        if ($c[1] ne 'topmed') { next; }    # Not interested in anything but topmed
        if ($c[4] ne 'topmed' && $c[4] ne 'schelcj') { $nottopmed++; next; }    # Not topmed user
        $partitions{$c[1]} = 1;
        if ($c[5] eq 'PD') {        # Queued
            push @{$queued{$c[1]}{data}},$l;
            $queued{$c[1]}{count}++;
            $qos{$c[2]}{queued}++;
            next;
        }
        if ($c[5] eq 'R') {         # Running
            push @{$running{$c[1]}{data}},$l;
            $running{$c[1]}{count}++;
            $qos{$c[2]}{running}++;
            next;
        }
    }
    
    #   Show summary of partitions
    print "Partition Summary\n";
    foreach my $p (sort keys %partitions) {
        if (! defined($queued{$p}{count}))  { $queued{$p}{count} = 0; }
        if (! defined($running{$p}{count})) { $running{$p}{count} = 0; }
        printf("  %-18s %3d running / %3d queued / %3d not user topmed \n", $p, $running{$p}{count}, $queued{$p}{count}, $nottopmed); 
    }
    print "\n";

    #   Get summary of IO data for key nodes
    my @iostatinfo = ();
    if (opendir(my $dh, "$opts{netdir}/$opts{consoledir}")) {
        @iostatinfo = grep { /iostat.topmed\d*.txt/ } readdir($dh);
        closedir $dh;
        print "I/O Load          Incoming MB/s   Working MB/s\n" .
            "         IOWait   reads  writes   reads  writes    Time\n";
        foreach my $f (sort @iostatinfo) {
            if (open(IN,"$opts{netdir}/$opts{consoledir}/$f")) {
                #   Linux 3.13.0-85-generic (topmed5) 	07/21/2016 	_x86_64_	(120 CPU)
                #
                #   07/21/2016 09:54:01 AM
                #   avg-cpu:  %user   %nice %system %iowait  %steal   %idle
                #          32.50    0.00    2.17    2.95    0.00   62.38
                #
                #   Device:            tps    MB_read/s    MB_wrtn/s    MB_read    MB_wrtn
                #   md20           3288.32        78.27        18.02  498062274  114686675
                #   md21           3196.36        91.78         1.36  584041063    8643316                my ($datestamp, 
                $_ = <IN>;                  # blank
                $_ = <IN>;                  # Device etc
                my ($date, $tod) = split(' ', <IN>);    # Get date and time
                $_ = <IN>;                  # avg-cpu etc
                my ($user, $nice, $system, $iowait) = split(' ', <IN>);     # Get iowait
                $_ = <IN>;                  # blank
                $_ = <IN>;                  # Device etc
                my ($dev1, $tps1, $MB_reads1, $MB_writes1) = split(' ', <IN>);
                my ($dev2, $tps2, $MB_reads2, $MB_writes2) = split(' ', <IN>);
                close(IN);
                if ($f =~ /(topmed\d*)\./) {
                    printf("%-7s %6s%% %7s %7s %7s %7s  %s\n",
                        $1, $iowait, $MB_reads1, $MB_writes1, $MB_reads2, $MB_writes2, $tod);
                }
            }
        }
        print "\n";
    }

    #   Show summary of QOS
    print "QOS Summary\n";
    foreach my $q (sort keys %qos) {
        if (! defined($qos{$q}{running}))  { $qos{$q}{running} = 0; }
        if (! defined($qos{$q}{queued}))   { $qos{$q}{queued} = 0; }
        my $qq = $q;                    # Patch normal to be something meaningful?
        if ($q eq 'normal') { $qq = 'default_qos'; }
        printf("  %-18s %3d running / %3d queued\n", $qq, $qos{$q}{running}, $qos{$q}{queued}); 
    }
    print "\n";

    #   Show summary of running on partitions
    print "Running jobs per host\n";
    my $format = '    %-11s  %s';
    foreach my $p (sort keys %partitions) {
        my %hosts = ();
        foreach my $l (@{$running{$p}{data}}) {
            my @c = split(' ', $l);
            $c[3] =~ s/\d+\-//;             # Remove bamid from jobname
            $c[3] =~ s/^NWD\d+/NWD/;        # This is remapping job
            $hosts{$c[8]}{$c[3]}++;         # Increment $hosts{topmed2}{cram}
        }
        print "   Partition $p:\n";
        my $format = '    %-11s  %s';
        printf ($format . "\n", 'Host', 'Job types and count');
        foreach my $h (sort keys %hosts) {
            my $s = '';
            foreach my $jobname (sort keys %{$hosts{$h}}) {
                $s .= sprintf('  %-12s', $jobname . "  [" . $hosts{$h}{$jobname} . '] ');
            }
            printf ($format . "\n", $h, $s);
        }
    }
    $format = "    %-8s  %s\n";
    print "  Partition nomosix jobs on topmed hosts:\n";
    printf ($format, 'User', 'Job count');
    foreach my $u (sort keys %mosixrunning) {
        printf($format, $u, "  [" . $mosixrunning{$u} . '] ');
    }
    print "\n";

    #   Show longest running jobs
    print "Longest running jobs per partition\n";
    foreach my $p (sort keys %partitions) {
        if (! defined($running{$p}{data})) { next; }
        print "  $p:\n";
        my %longjobs = ();
        foreach my $l (@{$running{$p}{data}}) {
            my @c = split(' ', $l);
            if ($c[3] =~ /^NWD/) { next; }
            my $t = $c[6];                  # Normalize all times to dd-hh:mm:ss
            my $days = '00-';
            if ($t =~ /^(\d+\-)(.+)/) { $days = $1; $t = $2; }
            if ($t =~ /^(\d+):(\d+):(\d+)$/) { $t = sprintf('%02d:%02d:%02d', $1, $2, $3); }
            elsif ($t =~ /^(\d+):(\d+)$/) { $t = sprintf('00:%02d:%02d', $1, $2); }
            elsif ($t =~ /^:(\d+)$/) { $t = sprintf('00:00:%02d', $1); }
            $longjobs{$t} .= '  ' . $l . "\n";      # Descriptions of longest running job
        }
        my $i = $opts{maxlongjobs};
        foreach my $t (reverse sort keys %longjobs) {
            foreach (split("\n", $longjobs{$t})) {
                print ReFormatPartitionData($_);
                $i--;
                if ($i <= 0) { last; }
            }
            if ($i <= 0) { last; }
        }
    }
    print "\n";
    return;
}

#==================================================================
# Subroutine:
#   ReFormatPartitionData($str)
#
#   Return reformatted partition data
#==================================================================
sub ReFormatPartitionData {
    my ($str) = @_;
    my $fmt = '    %-10s %-9s %-15s %-12s %s  %s';
    my @c = split(' ', $str);
    my $nwdid = '?';
    my $dir = '?';
    if ($c[3] =~ /^(\d+)/) {
        my $bamid=$1;
        my $sth = ExecSQL("SELECT runid,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
        if ($sth) {
            my $href = $sth->fetchrow_hashref;
            if (defined($href->{expt_sampleid})) { $nwdid = $href->{expt_sampleid}; }
            $sth = ExecSQL("SELECT dirname FROM $opts{runs_table} WHERE runid=$href->{runid}");
            if ($sth) {
                $href = $sth->fetchrow_hashref;
                $dir = $href->{dirname};
            }
        }
    }
    return sprintf($fmt, $c[0], $nwdid, $c[2], $c[3], $c[6], $dir)  . "\n";
}

#==================================================================
# Subroutine:
#   $path = WherePath($bamid, $set, $extra1)
#
#   Print paths to various things for bamid based on $set
#     bam       Print directory for a bam
#     backup    Print directory for backups file
#     qcresults Print directory where qc.results for a bamfile
#     console   Print directory where SLURM console output lives
#     b37       Print directory for remapped b37
#     b38       Print directory for remapped b38
#==================================================================
sub WherePath {
    my ($bamid, $set, $extra1) = @_;
    if ((! defined($set) || ! $set)) { $set = 'unset'; }    # Default

    $bamid = GetBamid($bamid);

    #   Get values of interest from the database
    my $sth = ExecSQL("SELECT runid,piname,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $piname = $href->{piname};
    my $nwdid = $href->{expt_sampleid};
    my $runid = $href->{runid};
    $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   BAM is in one of those $opts{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamhost = 'none';
        my $bamfdir = abs_path("$opts{netdir}/$opts{incomingdir}/$centername/$rundir");
        print "$bamfdir\n";
        exit;
    }
 
    #   Find where the backup file lives
    if ($set eq 'backup') {
        my $bakbamfdir = abs_path("$opts{netdir}/$opts{backupsdir}/$centername/$rundir");
        print "$bakbamfdir\n";
        exit;
    }

    #   Print where qc.results are 
    if ($set eq 'qcresults') {
        my $qcdir = abs_path("$opts{netdir}/$opts{qcresultsdir}/$centername/$rundir");
        print "$qcdir\n";                       # File might not really exist
        exit;
    }
 
    #   Print where SLURM console output can be found
    if ($set eq 'console') {
        my $dir = abs_path("$opts{netdir}/$opts{consoledir}");
        print "$dir\n";
        exit;
    }
 
    #   Try to guess where the b37 remapped CRAM lives
    if ($set eq 'b37') {
        my $b37fdir = '';
        my $b37file = '';
        foreach ('', '2', '3', '4', '5', '6') {
            $b37fdir = abs_path("$opts{netdir}$_/$opts{resultsdir}/$centername/$piname/$nwdid/bams");
            if ($b37file) {
                print "$b37fdir\n";
                exit;
            }
        }
        exit;
    }

    #   Try to guess where the b38 remapped CRAM lives -- this likely needs to be corrected
    if ($set eq 'b38') {
        die "$Script - b38 paths are not known yet\n";
    }

    die "$Script - Unknown Where option '$set'\n";
}

#==================================================================
# Subroutine:
#   $host = WhatHost($bamid, $set)
#
#   Print host name for the path to various things for bamid based on $set
#     bam       Print host where BAM actually exists
#     backup    Print host for backups file
#     qcresults Print host where qc.results for a bamfile is
#==================================================================
sub WhatHost {
    my ($bamid, $set, $extra1) = @_;
    if ((! defined($set) || ! $set)) { $set = 'unset'; }    # Default

    $bamid = GetBamid($bamid);

    #   Get values of interest from the database
    my $sth = ExecSQL("SELECT runid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};
    $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   BAM is in one of those $opts{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamhost = 'none';
        my $bamfdir = abs_path("$opts{netdir}/$opts{incomingdir}/$centername/$rundir");
        if (! $bamfdir) { exit; }
        if ($bamfdir =~ /\/net\/([^\/]+)\//) { $bamhost = $1; }
        print "$bamhost\n";
        exit;
    }
 
    #   Find where the backup CRAM lives and show host
    if ($set eq 'backup') {
        my $backuphost = 'none';
        my $bakbamfdir = abs_path("$opts{netdir}/$opts{backupsdir}/$centername/$rundir");
        if ($bakbamfdir =~ /\/net\/([^\/]+)\//) { $backuphost = $1; }
        print "$backuphost\n";
        exit;
    }

    #   Print where qc.results are and show host
    if ($set eq 'qcresults') {
        my $qchost = 'none';
        my $qcdir = abs_path("$opts{netdir}/$opts{qcresultsdir}/$centername/$rundir");
        if ($qcdir =~ /\/net\/([^\/]+)\//) { $qcdir = $1; }
        print "$qcdir\n";
        exit;
    }
 
    die "$Script - Unknown Where option '$set'\n";
}

#==================================================================
# Subroutine:
#   $filepath = WhereFile($bamid, $set)
#
#   Print paths to various files for bamid based on $set
#     bam       Print file path where BAM actually exists, no symlink,
#     backup    Print file path for backups
#     qcresults Print file path to selfSM in qc.results for a bamfile
#==================================================================
sub WhereFile {
    my ($bamid, $set) = @_;
    if ((! defined($set) || ! $set)) { $set = 'unset'; }    # Default

    $bamid = GetBamid($bamid);

    #   Get values of interest from the database
    my $sth = ExecSQL("SELECT runid,bamname,cramname,piname,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $bamname = $href->{bamname};
    my $cramname = $href->{cramname} || 'CRAMNAME_not_SET';
    my $runid = $href->{runid};
    $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$centerid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   BAM is in one of those $opts{netdir} trees where without a symlink
    if ($set eq 'bam') {
        my $bamfile = abs_path("$opts{netdir}/$opts{incomingdir}/$centername/$rundir");
        if ($bamfile) { $bamfile .= '/' . $bamname; }
        print "$bamfile\n";
        exit;
    }
 
    #   Find where the backup CRAM lives
    if ($set eq 'backup') {
        my $bakfile = abs_path("$opts{netdir}/$opts{backupsdir}/$centername/$rundir");
        if ($bakfile) { $bakfile .= '/' . $cramname; }
        print "$bakfile\n";
        exit;
    }

    #   Print where qc.results we are interested are
    if ($set eq 'qcresults') {
        my $qcdir = abs_path("$opts{netdir}/$opts{qcresultsdir}/$centername/$rundir");
        $bamname =~ s/\.bam//;                   # Remove the extension
        $bamname =~ s/\.cram//;
        print "$qcdir/$bamname.vb.selfSM\n";
        exit;
    }
 
    die "$Script - Unknown Where option '$set'\n";
}

#==================================================================
# Subroutine:
#   Permit ($fcn, $op, $center, $run)
#
#   Controls a database of enabled or disabled operations for a center/run
#
#   Permit ('test', $op, $bamid)
#
#   Returns a boolean if an operation for a center/run is allowed
#==================================================================
sub Permit {
    my ($fcn) = @_;
    my ($sql, $sth, $href, $rowsofdata);
    shift(@_);

    #   Test is some operation is allowed
    #   topmedcmd.pl permit test bai 4955
    if ($fcn eq 'test') {
        my ($op, $bamid) = @_;

        #   Get runid and centerid for this bam
        my ($centerid, $runid) = GetBamidInfo($bamid);

        #   This is a small table, read it all
        $sth = ExecSQL("SELECT runid,centerid,operation FROM $opts{permissions_table}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { exit(1); }                 # Table empty, permitted
        for (my $i=1; $i<=$rowsofdata; $i++) {
            $href = $sth->fetchrow_hashref;
            if ($centerid eq $href->{centerid} || $href->{centerid} eq '0') {
                if ($runid eq $href->{runid} || $href->{runid} eq '0') {
                    if ($op eq $href->{operation} || $href->{operation} eq 'all') {
                        exit;                           # Operation not permitted
                    }
                }
            }
        }
        exit(1);                                    # TRUE, permitted  
    }

    #   Enable a permission by deleting a database entry
    #   e.g. topmedcmd.pl permit remove id
    if ($fcn eq 'remove') {
        my ($id) = @_;
        $sql = "SELECT * FROM $opts{permissions_table} WHERE id=$id";
        $sth = ExecSQL($sql);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "Permission '$id' does not exist\n"; }
        $href = $sth->fetchrow_hashref;
        my $s = "$href->{centername} / $href->{dirname} / $href->{operation}";
        $sql = "DELETE FROM $opts{permissions_table} WHERE id='$id'";
        $sth = ExecSQL($sql);
        if (! $sth) { die "Failed to remove permission '$fcn' for 'id=$id $s'\nSQL=$sql"; }
        print "Deleted permission control for '$s'\n";
        exit;
    }

    if ($fcn ne 'add') { die "Unknown function '$fcn'\n"; }

    #   Disable a permission by adding to the database entry
    #   e.g. # topmedcmd.pl permit add [bai [broad [2015sep18]]]
    my ($op, $center, $run) = @_;
    if (! $op)       { $op = 'all'; }           # Set defaults
    if (! $center)   { $center = 'all'; }
    if (! $run)      { $run = 'all'; }
    my ($centerid, $runid) = (0, 0);

    if (! exists($VALIDOPS{$op})) { die "Operation '$op' is not known\n"; }

    #   Verify op and run and center setting runid and centerid
    if ($center ne 'all') {
        $sql = "SELECT centerid FROM $opts{centers_table} WHERE centername='$center'";
        $sth = ExecSQL($sql);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "Center '$center' is not known\n"; }
        $href = $sth->fetchrow_hashref;
        $centerid = $href->{centerid};
    }

    if ($run ne 'all') {
        $sql = "SELECT runid FROM $opts{runs_table} WHERE dirname='$run'";
        $sth = ExecSQL($sql);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "Run '$run' is not known\n"; }
        $href = $sth->fetchrow_hashref;
        $runid = $href->{runid};
    }

    $sql = "INSERT INTO $opts{permissions_table} " .
        "(centername,dirname,centerid,runid,operation) " .
        "VALUES ('$center', '$run', $centerid, $runid, '$op')";
    $sth = ExecSQL($sql);
    if (! $sth) { die "Failed to add permission '$fcn' for '$op'\nSQL=$sql"; }
    print "Added permission '$center / $runid / $op'\n";
    exit;
}

#==================================================================
# Subroutine:
#   Set($bamid, $col, $val)
#
#   Set a database column
#==================================================================
sub Set {
    my ($bamid, $col, $val) = @_;
    my ($sth, $rowsofdata, $href);

    #   This could be a run name
    $sth = ExecSQL("SELECT runid FROM $opts{runs_table} WHERE dirname='$bamid'", 0);
    $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        if ($rowsofdata > 1) {
            die "$Script - Eeek, there are $rowsofdata runs named '$bamid'\n";
        }
        $href = $sth->fetchrow_hashref;
        ExecSQL("UPDATE $opts{runs_table} SET $col='$val' WHERE runid='$href->{runid}'");
        return;
    }

    #   This is bamid or nwdid
    $bamid = GetBamid($bamid);

    #   Make sure this is a bam we know
    $sth = ExecSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    if ($col eq 'nwdid') { $col = 'expt_sampleid'; }
    ExecSQL("UPDATE $opts{bamfiles_table} SET $col='$val' WHERE bamid=$bamid");
}

#==================================================================
# Subroutine:
#   ShowRun($runname)
#
#   Generate list of nwdids (and maybe bamids) for a run
#==================================================================
sub ShowRun {
    my ($runname) = @_;
    my $s;

    my $sth = ExecSQL("SELECT runid FROM $opts{runs_table} WHERE dirname='$runname'", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Unknown run '$runname'\n"; }
    my $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};
    $sth = ExecSQL("SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    for (my $i=1; $i<=$rowsofdata; $i++) {
        $href = $sth->fetchrow_hashref;
        $s = $href->{expt_sampleid};
        if ($opts{'with-bamid'}) { $s .= ' ' . $href->{bamid}; }
        print $s . "\n";
    }
}

#==================================================================
# Subroutine:
#   Show($fcn, $bamid, $col)
#
#   Generate list of information from the database
#   or show the value for a column
#==================================================================
sub Show {
    my ($bamid, $col) = @_;
    my ($sth, $rowsofdata, $href);

    #   This could be a run name.g Does not start with NWD, not all digits
    if ($bamid !~ /^NWD/ && $bamid =~ /\D+/) {
        $sth = ExecSQL("SELECT runid,$col FROM $opts{runs_table} WHERE dirname='$bamid'", 0);
        if ($sth) {
            $rowsofdata = $sth->rows();
            if ($rowsofdata) {
                $href = $sth->fetchrow_hashref;
                print $href->{$col} . "\n";
                return;
            }
        }
    }

    #   This is bamid or nwdid
    $bamid = GetBamid($bamid);
    
    $sth = ExecSQL("SELECT bamid,runid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};

    #   Get run if asked for it
    if ($col eq 'run') {
        $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$href->{runid}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - RUNID '$href->runid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        print $href->{dirname} . "\n";
        return;
    }

    #   Get center if asked for it
    if ($col eq 'center') {
        $sth = ExecSQL("SELECT centerid FROM $opts{runs_table} WHERE runid=$href->{runid}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - RUNID '$href->runid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;

        $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - CENTERID '$href->centerid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        print $href->{centername} . "\n";
        return;
    }

    #   Show everything in YAML format
    if ($col eq 'yaml') {
        $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$runid");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - RUNID '$runid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        my $run = $href->{dirname};

        $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - CENTERID '$href->centerid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        my $center = $href->{centername};

        $sth = ExecSQL("SELECT * FROM $opts{bamfiles_table} WHERE bamid=$bamid");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - BAM '$bamid' or column '$col' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        $href->{run} = $run;
        $href->{center} = $center;
        print YAML::Dump($href);
        return;
    }

    #   Get value of column we asked for
    $sth = ExecSQL("SELECT $col FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' or column '$col' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    if (defined($href->{$col})) { print $href->{$col} . "\n"; }
}

#==================================================================
# Subroutine:
#   GetBamid($bamid)
#
#   Return bamid for bamid or expt_sampleid
#==================================================================
sub GetBamid {
    my ($bamid) = @_;
    my ($sth, $rowsofdata, $href);
    if ($bamid =~ /^\d+$/) { return $bamid; }

    if ($bamid =~ /^NWD/){
        $sth = ExecSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE expt_sampleid='$bamid'");
        $rowsofdata = $sth->rows();
        if ($rowsofdata) {
            $href = $sth->fetchrow_hashref;
            return $href->{bamid};
        }
    }
    if ($bamid !~ /^\d+$/){
        die "$Script - Invalid bamid or NWDID ($bamid). Try '$Script -help'\n";
    }
    return $bamid;
}

#==================================================================
# Subroutine:
#   ($centerid, $runid) = GetBamidInfo($bamid)
#
#   Return id for the center and run for a bamid
#==================================================================
sub GetBamidInfo {
    my ($bamid) = @_;

    if ($bamid !~ /^\d+$/) { return (0,0); }     # No bamid, no ids
    my $sth = ExecSQL("SELECT runid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};
    $sth = ExecSQL("SELECT centerid FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' has no center?  How'd that happen?\n"; }
    $href = $sth->fetchrow_hashref;
    return ($href->{centerid}, $runid);
}

#==================================================================
# Subroutine:
#   ExecSQL($sql, $die)
#
#   Execute SQL.  Keep trying if connection lost.
#
#   Returns handle for SQL
#==================================================================
sub ExecSQL {
    my ($sql, $die) = @_;
    if ($opts{persist}) { return PersistDoSQL($opts{realm}, $sql); }
    return DoSQL($sql, $die);
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedcmd.pl - Update the database for NHLBI TopMed

=head1 SYNOPSIS

  topmedcmd.pl mark 33 arrived completed   # BAM has arrived
  topmedcmd.pl mark NWD00234  arrived completed   # Same BAM has arrived
  topmedcmd.pl unmark 33 arrived           # Reset BAM has arrived

  topmedcmd.pl set 33 jobidqplot 123445    # Set jobidqplot in bamfiles
  topmedcmd.pl set NWD123433 jobidqplot 123445    # Set jobidqplot in bamfiles
  topmedcmd.pl set 2016apr20 offsite N     # Set 2016apr20 in runs

  topmedcmd.pl show 2199 state_cram        # Show a column
  topmedcmd.pl show 2199 center            # Show center 
  topmedcmd.pl show NWD00234 run           # Show run
  topmedcmd.pl show 2016apr20 offsite      # Show offsite for run
  topmedcmd.pl show NWD00234 yaml          # Show everything known about a bamid

  topmedcmd.pl showrun 20150604            # Show all NWDID for run
  topmedcmd.pl -bamid showrun 20150604     # Show all NWDID and bamid for run
 
  topmedcmd.pl wherepath 2199 bam          # Returns real path to bam
  topmedcmd.pl wherepath 2199 backup       # Returns path to backups directory
  topmedcmd.pl wherepath 2199 qcresults    # Returns path to directory for qc.results
  topmedcmd.pl wherepath 2199 console      # Returns path to directory for SLURM output
  topmedcmd.pl wherepath NWD00234 b37      # Returns path to remapped b37 directory and to file
  topmedcmd.pl wherepath 2199 b38          # Returns path to remapped b38 directory and to file

  topmedcmd.pl whathost 2199 bam           # Returns host for bam
  topmedcmd.pl whathost 2199 backup        # Returns host for backups directory
  topmedcmd.pl whathost 2199 qcresults     # Returns host for directory for qc.results

  topmedcmd.pl wherefile 2199 bam          # Returns path to bam file (may not exist)
  topmedcmd.pl wherefile 2199 backup       # Returns path for backups file (may not exist)
  topmedcmd.pl wherefile 2199 qcresults    # Returns path for qc.results *.vb.SelfSM file (may not exist)

  topmedcmd.pl permit add bai braod 2015oct18   # Stop bai job submissions for a run
  topmedcmd.pl permit remove 12             # Remove a permit control
  topmedcmd.pl permit test bai 4567         # Test if we should submit a bai job for one bam

  topmedcmd.pl -maxlongjobs 35 squeue       # Show what is running
=head1 DESCRIPTION

This program supports simple commands to set key elements of the NHLBI database.
The queue of tasks is kept in a MySQL database.
See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-bamid  -with-bamid>

Include the bamid for the output from showrun.

=item B<-help>

Generates this output.

=item B<-maxlongjobs N>

Specifies how many of the longest running jobs to show.
This defaults to B<10>.

=item B<-persist>

Specifies that SQL connection errors will not fail, but will be retried many times
before finally failing.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program are grouped into several groups which are used
to deal with specific sets of information in the monitor databases.

B<mark bamid|nwdid dirname  [verb] [state]>
Use this to set the state for a particular BAM file.
You may specify the bamid or the NWDID.
Mark will set a date for the process (e.g. arrived sets state_arrive)
and unmark will set that entry to NULL.
The list of verbs and states can be seen by B<perldoc topmedcmd.pl>.

B<permit enable/disable operation center run>
Use this to control the database which allows one enable or disable topmed operations
(e.g. backup, verify etc) for a center or run.
Use B<all> for all centers or all runs or all operations.

B<permit test operation bamid>
Use this to test if an operation (e.g. backup, verify etc) may be submitted 
for a particular bam.

B<set bamid|nwdid|dirname columnname value>
Use this to set the value for a column for a particular BAM file
or run.

B<send2ncbi filelist>
Use this to copy data to NCBI with ascp.

B<show bamid|nwdid|dirname colname|center|run|yaml>
Use this to show information about a particular bamid (or expt_sampleid)
or run name.
Use 'yaml' to display everything known about the bam of interest.

B<showrun runname>
Use this to show a list of all the NWDIDs for a run.
The option B<-bamid> will cause the bamid to be shown in addition.

B<unmark bamid|nwdid [verb]>
Use this to reset the state for a particular BAM file to the default
database value.

B<whatnwdid bamid|nwdid>
Use this to get some details for a particular bam.

B<wherepath bamid|nwdid bam|backup|qcresults|console|b37|b38>
If B<bam> was specified, display the path to the real bam file.

If B<backup> was specified, display the path to the backup directory.

If B<qcresults> was specified, display the path to the directory where
the qc.results for this bamid will be.

If B<console> was specified, display the path to the directory where
the SLURM console output.

If B<b37> was specified, display the path to the directory of remapped data for build 37 results can be found.

If B<b38> was specified, display the path to the directory of remapped data for build 37 results can be found.

B<whathhost bamid|nwdid bam|backup|qcresults>
returns the host for the bam, backup or qc.results for the bam file.

B<wherefile bamid|nwdid bam|backup|qcresults>
returns the path to the file for the bam, backup or qc.results. This file may not exist


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

