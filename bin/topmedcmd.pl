#!/usr/bin/perl
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
use lib (
  qq($FindBin::Bin),
  qq($FindBin::Bin/../lib),
  qq($FindBin::Bin/../lib/perl5),
  qq($FindBin::Bin/../local/lib/perl5),
  qq(/usr/cluster/topmed/lib/perl5),
  qq(/usr/cluster/topmed/local/lib/perl5),
);

use TopMed_Get;
use My_DB;
use Getopt::Long;
use Cwd qw(realpath abs_path);
use YAML;
use POSIX qw(strftime);

my $NOTSET = 0;                     # Not set
my $REQUESTED = 1;                  # Task requested
my $SUBMITTED = 2;                  # Task submitted to be run
my $STARTED   = 3;                  # Task started
my $DELIVERED = 19;                 # Data delivered, but not confirmed
my $COMPLETED = 20;                 # Task completed successfully
my $CANCELLED = 89;                 # Task cancelled
my $FAILED    = 99;                 # Task failed

#   How to add a column to this program
#       Add verbname => db column name to %VALIDVERBS
#       Add operation verb to %VALIDOPS
my %VALIDVERBS = (                  # Valid verbs to database colum
    arrive     => 'state_arrive',
    verify     => 'state_verify',
    qplot      => 'state_qplot',
    cram       => 'state_cram',   
    gcepush    => 'state_gce38push',
    gcepull    => 'state_gce38pull',
    b37        => 'state_b37',   
    b38        => 'state_b38',     
    bcf        => 'state_bcf',     
    ncbiexpt   => 'state_ncbiexpt',
    ncbiorig   => 'state_ncbiorig',
    ncbib37    => 'state_ncbib37',
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

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    centers_table => 'centers',
    permissions_table => 'permissions',
    runs_table => 'runs',
    netdir => '/net/topmed',
    backupsdir => 'working/backups/incoming/topmed',
    consoledir => 'working/topmed-output',
    ascpcmd => "/usr/cluster/bin/ascp -i ".
        "/net/topmed/incoming/study.reference/send2ncbi/topmed-2-ncbi.pri -l 800M -k 1",
    ascpdest => 'asp-um-sph@gap-submit.ncbi.nlm.nih.gov:protected',
    squeuecmd => "/usr/cluster/bin/squeue -o '%.9i %.15P %.15q %.14j %.8u %.2t %.8M %.6D %R %.9n'",    
    maxlongjobs => 10,
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose maxlongjobs=i persist with-id|bamid emsg=s
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options] [-persist]";
    my $verbs = join(' ', sort keys %VALIDVERBS);
    my $requests = join(' ', sort keys %VALIDSTATUS);
    warn "$m mark bamid|nwdid action newstatus\n" .
        "    Where action: $verbs\n" .
        "    Where newstatus: $requests\n" .
        "  or\n" .
        "$m unmark bamid|nwdid [same list as mark]\n" .
        "  or\n" .
        "$m set bamid|nwdid|dirname colname value\n" .
        "  or\n" .
        "$m setdate bamid|nwdid colname file\n" .
        "  or\n" .
        "$m show bamid|nwdid colname|run|center|yaml\n" .
        "  or\n" .
        "$m list centers\n" .
        "$m list runs centername\n" .
        "$m [-with-id] list samples runname\n" .
        "  or\n" .
        "$m export\n" .
        "  or\n" .
        "$m send2ncbi files\n" .
        "  or\n" .
        "$m squeue\n" .
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
if ($fcn eq 'setdate')   { SetDate(@ARGV); exit; }
if ($fcn eq 'show')      { Show(@ARGV); exit; }
if ($fcn eq 'list')      { List(@ARGV); exit; }
if ($fcn eq 'export')    { Export(@ARGV); exit; }
if ($fcn eq 'send2ncbi') { Send2NCBI(@ARGV); exit; }
if ($fcn eq 'squeue')    { SQueue(@ARGV); exit; }
if ($fcn eq 'whatbamid') { WhatBAMID(@ARGV); exit; }
if ($fcn eq 'whatnwdid') { WhatNWDID(@ARGV); exit; }

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
        if ($col eq 'state_arrive') {       # Used by Kevin for tracking samples
            ExecSQL("UPDATE $opts{bamfiles_table} SET datearrived='" . time() . "' WHERE bamid=$bamid");
        }
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            ExecSQL("UPDATE $opts{bamfiles_table} SET datemapping_b37='" . time() . "' WHERE bamid=$bamid");
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
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            ExecSQL("UPDATE $opts{bamfiles_table} SET datemapping_b37='-1' WHERE bamid=$bamid");
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
        if (exists($opts{emsg})) {
            print $opts{emsg} . "\n";
            ExecSQL("UPDATE $opts{bamfiles_table} SET emsg=\"$op $state -- $opts{emsg}\" WHERE bamid=$bamid");
        }
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
    my @cols = qw(cramname studyname piname expt_sampleid state_b37 state_b38 datayear
        cramflagstat build checksum);
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
                #   See if this has been verified
                if ($href->{state_verify} != $COMPLETED) { next; }
                #   Show data for this CRAM
                my $f = "$opts{netdir}/$opts{backupsdir}/$centername/$dirname/" .
                    $href->{cramname};
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
    my $sth = ExecSQL("SELECT runid,bamid,piname,datayear FROM $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - NWDID '$nwdid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $bamid = $href->{bamid};
    my $datayear = $href->{datayear};
    my $piname = $href->{piname};
    $sth = ExecSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$href->{runid}");
    $href = $sth->fetchrow_hashref;
    my $run = $href->{dirname};
    $sth = ExecSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
    $href = $sth->fetchrow_hashref;
    my $center = uc($href->{centername});
    print "$nwdid $bamid can be found in run $run PI $piname for center $center year $datayear\n";
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
#squeuecmd => "/usr/cluster/bin/squeue -o '%.10i %.18P %.15q %.18j %.8u %.2t %.10M %.9n'",
sub SQueue {
    #my ($bamid, $set) = @_;
    my $cmd = $opts{squeuecmd} . '| grep topmed';
    my %partitions = ();            # Unique partition names
    my %qos = ();                   # Unique QOS names
    my %qosheld = ();               # Count of held jobs in QOS
    my %running = ();
    my %mosixrunning = ();          # Counts of jobs by nomosix users
    my %queued = ();
    my %nottopmed = ();             # Not my user list
    foreach my $l (split('\n', `$cmd`)) {
        my @c = split(' ', $l);
        if ($c[1] eq 'nomosix' && $c[8] =~ /topmed/) {    # nomosix on topmed node
            $mosixrunning{$c[4]}++;    # Count of running jobs for this user
            next;
        }
        if ($c[1] !~ /^topmed/) { next; }    # Only want my partition of interest
        if ($c[4] ne 'topmed' && $c[4] ne 'schelcj') {  # Save partition and not topmed user
            $nottopmed{$c[1]}{$c[4]} = 1;
            next;
        }
        $partitions{$c[1]} = 1;
        if ($c[5] eq 'PD') {        # Queued
            push @{$queued{$c[1]}{data}},$l;
            $queued{$c[1]}{count}++;
            $qos{$c[2]}{queued}++;
            if ($l =~ /held state/) { $qos{$c[2]}{held}++; }
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
        my $s = join(' ',sort keys %{$nottopmed{$p}});
        my $k = scalar(keys %{$nottopmed{$p}});
        printf("  %-18s %3d running / %3d queued / %3d foreign user: $s \n",
            $p, $running{$p}{count}, $queued{$p}{count}, $k);
    }
    print "\n";

    #   Show summary of QOS
    print "QOS Summary\n";
    foreach my $q (sort keys %qos) {
        my $s = $qos{$q}{held} || '';
        if ($s) { $s .= ' jobs held'; }
        if (! defined($qos{$q}{running}))  { $qos{$q}{running} = 0; }
        if (! defined($qos{$q}{queued}))   { $qos{$q}{queued} = 0; }
        my $qq = $q;                    # Patch normal to be something meaningful?
        if ($q eq 'normal') { $qq = 'default_qos'; }
        printf("  %-18s %3d running / %3d queued   %s\n",
            $qq, $qos{$q}{running}, $qos{$q}{queued}, $s); 
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
    print "\n";
    if (%mosixrunning) {
        $format = "    %-8s  %s\n";
        print "  Partition nomosix jobs on topmed hosts:\n";
        printf ($format, 'User', 'Job count');
        foreach my $u (sort keys %mosixrunning) {
            printf($format, $u, "  [" . $mosixrunning{$u} . '] ');
        }
        print "\n";
    }
    
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

    #   Get summary of IO data for key nodes
    my @iostatinfo = ();
    if (opendir(my $dh, "$opts{netdir}/$opts{consoledir}")) {
        @iostatinfo = sort sorthost grep { /iostat.topmed\d*.txt/ } readdir($dh);
        closedir $dh;
        my $in;
        my $prefix = ' ' x 12;
        my $date = $prefix;
        my $hdr = 0;
        my $iohdr = '   %user %iowait  %idle';
        my $iodata = '';
        print "    I/O activity for /incoming and /working filesystems per host\n";
        foreach my $f (@iostatinfo) {
            my $host = $prefix;
            if (open($in,"$opts{netdir}/$opts{consoledir}/$f")) {
                if ($f =~ /(topmed\d*)/) { $host = sprintf('%-12s', $1); }
                #   Linux 4.2.0-42-generic (topmed) 	02/10/2017 	_x86_64_	(120 CPU)
                #
                #   avg-cpu:  %user   %nice %system %iowait  %steal   %idle
                #              0.35    0.00    1.66    0.55    0.00   97.44
                #
                #   Device:            tps    MB_read/s    MB_wrtn/s    MB_read    MB_wrtn
                #   md20              0.40         0.00         0.01          0          0
                #   md21           7782.60       213.60        13.17       1068         65
                while (<$in>) {
                    if (/Linux.+(\d\d\/\d\d\/20\d\d)/) { $date = sprintf('%-12s', $1); next; }
                    if (/Device/ && (! $hdr)) {
                        print $date . ReFormatDeviceData($_) . $iohdr . "\n";
                        $hdr++;
                        next;
                    }
                    if (/avg-cpu/) {
                        $_ = <$in>;             # Get actual numbers
                        $iodata = ReFormatIOData($_);
                        next;
                    }
                    if (/md2\d/) {
                        print $host . ReFormatDeviceData($_) . $iodata ."\n";
                        $host = $prefix;
                        $iodata = '';
                        next;
                    }
                }
                close($in);
            }
        }
        print "\n";
    }
    return;
}

#==================================================================
# Subroutine:
#   sorthost    Used to sort topmed hosts in numnerical order
#
#   Return reformatted device data
#==================================================================
sub sorthost {
    my $cmp1 = 0;
    my $cmp2 = 0;
    if ($a =~ /topmed(\d+)/) { $cmp1 = sprintf('%02d',$1); }
    if ($b =~ /topmed(\d+)/) { $cmp2 = sprintf('%02d',$1); }
    return ($cmp1 <=> $cmp2);
}

#==================================================================
# Subroutine:
#   ReFormatIOData($str)
#
#   Return reformatted device data
#==================================================================
sub ReFormatIOData {
    my ($str) = @_;
    my @c = split(' ', $str);
    return sprintf('   %5s %7s %6s', $c[0], $c[3], $c[5]);
}

#==================================================================
# Subroutine:
#   ReFormatDeviceData($str)
#
#   Return reformatted device data
#==================================================================
sub ReFormatDeviceData {
    my ($str) = @_;
    my @c = split(' ', $str);
    if ($c[0] eq 'md20') { $c[0] = '/incoming'; }
    else {
        if ($c[0] eq 'md21') { $c[0] = '/working'; }
    }
    return sprintf('%-9s %10s %10s', $c[0], $c[2], $c[3]);
}

#==================================================================
# Subroutine:
#   ReFormatPartitionData($str)
#
#   Return reformatted partition data
#==================================================================
sub ReFormatPartitionData {
    my ($str) = @_;
    my $fmt = '    %-9s %-9s %-12s %-10s %-9s %s  %s';
    my @c = split(' ', $str);
    my $nwdid = '?';
    my $dir = '?';
    if ($c[3] =~ /^(\d+)/) {
        my $bamid=$1;
        my $sth = ExecSQL("SELECT runid,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
        if ($sth) {
            my $href = $sth->fetchrow_hashref;
            if (exists($href->{expt_sampleid})) { $nwdid = $href->{expt_sampleid}; }
            if (exists($href->{runid})) {
                $sth = ExecSQL("SELECT dirname FROM $opts{runs_table} WHERE runid=$href->{runid}");
                if ($sth) {
                    $href = $sth->fetchrow_hashref;
                    $dir = $href->{dirname};
                }
            }
        }
    }
    return sprintf($fmt, $c[0], $nwdid, $c[2], $c[3], $c[8], $c[6], $dir)  . "\n";
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

    if (! defined($val)) { $val = ''; }
    if ($val ne 'NULL') { $val = "'$val'"; }    # Use quotes unless this is NULL

    #   This could be a run name
    $sth = ExecSQL("SELECT runid FROM $opts{runs_table} WHERE dirname='$bamid'", 0);
    $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        if ($rowsofdata > 1) {
            die "$Script - Eeek, there are $rowsofdata runs named '$bamid'\n";
        }
        $href = $sth->fetchrow_hashref;
        ExecSQL("UPDATE $opts{runs_table} SET $col=$val WHERE runid='$href->{runid}'");
        return;
    }

    #   This is bamid or nwdid
    $bamid = GetBamid($bamid);

    #   Make sure this is a bam we know
    $sth = ExecSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    if ($col eq 'nwdid') { $col = 'expt_sampleid'; }
    ExecSQL("UPDATE $opts{bamfiles_table} SET $col=$val WHERE bamid=$bamid");
}

#==================================================================
# Subroutine:
#   SetDate($bamid, $col, $val)
#
#   Set a database column to the date for a file
#==================================================================
sub SetDate {
    my ($bamid, $col, $val) = @_;
    my ($sth, $rowsofdata, $href);

    #   This is bamid or nwdid
    $bamid = GetBamid($bamid);

    #   Make sure this is a bam we know
    $sth = ExecSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    my @s = stat($val);
    if (! @s) { die "$Script - '$val' is not a known filename\n"; }
    my $datetime = strftime('%Y-%m-%d %H:%M:%S', localtime($s[10]));
    ExecSQL("UPDATE $opts{bamfiles_table} SET $col='$datetime' WHERE bamid=$bamid");
}

#==================================================================
# Subroutine:
#   List($fcn, $item)
#
#   Generate a list of data from the database
#   Where $fcn may be centers, runs or samples
#==================================================================
sub List {
    my ($fcn, $item) = @_;
    if (! $fcn) { die "$Script - List operator was not provided\n"; }
    my $s;

    if ($fcn eq 'centers') {            # Show all centers
        my $sth = ExecSQL("SELECT centerid,centername FROM $opts{centers_table}");
        my $rowsofdata = $sth->rows();
        for (my $i=1; $i<=$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref;
            $s = $href->{centername};
            if ($opts{'with-id'}) { $s .= ' ' . $href->{centerid}; }
            print $s . "\n";
        }
        return;
    }
    if ($fcn eq 'runs') {               # Show all runs for a center
        if (! $item) { die "$Script - Centername was not provided\n"; }
        my $sth = ExecSQL("SELECT centerid FROM $opts{centers_table} WHERE centername='$item'",0);
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - Unknown centername '$item'\n"; }
        my $href = $sth->fetchrow_hashref;
        my $centerid = $href->{centerid};
        $sth = ExecSQL("SELECT runid,dirname FROM $opts{runs_table} WHERE centerid=$centerid");
        $rowsofdata = $sth->rows();
        for (my $i=1; $i<=$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref;
            $s = $href->{dirname};
            if ($opts{'with-id'}) { $s .= ' ' . $href->{runid}; }
            print $s . "\n";
        }
        return;
    }
    if ($fcn eq 'samples') {            # Show all samples for a run
        if (! $item) { die "$Script - Runname was not provided\n"; }
        my $sth;
        if ($item =~ /^\d+$/) {        # Maybe runid provided
            $sth = ExecSQL("SELECT runid FROM $opts{runs_table} WHERE runid=$item", 0);
        }
        else {
            $sth = ExecSQL("SELECT runid FROM $opts{runs_table} WHERE dirname='$item'", 0);
        }
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - Unknown run '$item'\n"; }
        my $href = $sth->fetchrow_hashref;
        my $runid = $href->{runid};
        $sth = ExecSQL("SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} WHERE runid=$runid");
        $rowsofdata = $sth->rows();
        for (my $i=1; $i<=$rowsofdata; $i++) {
            $href = $sth->fetchrow_hashref;
            $s = $href->{expt_sampleid};
            if ($opts{'with-id'}) { $s .= ' ' . $href->{bamid}; }
            print $s . "\n";
        }
        return;
    }
    die "$Script - Unknown list function '$fcn'\n";
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

  topmedcmd.pl mark 33 arrive completed   # BAM has arrived
  topmedcmd.pl -emsg 'No file found' mark 33 cram failed   # Action failed, set error msg
  topmedcmd.pl mark NWD00234  arrive completed   # Same BAM has arrived
  topmedcmd.pl unmark 33 arrive           # Reset BAM has arrived

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
 
  topmedcmd.pl list centers                # Show all known centers
  topmedcmd.pl list runs broad             # Show all runs for center broad
  topmedcmd.pl list samples 2015dec02      # Show all samples for a run
 
  topmedcmd.pl export                      # Dump the database in a form Chris wanted
 
  topmedcmd.pl send2ncbi files             # Send files to NCBI
 
  topmedcmd.pl whatbamid bamname           # Show bamid for a sample name
 
  topmedcmd.pl whatnwdid bamid|nwdid       # Show details for a sample
 
  topmedcmd.pl -maxlongjobs 35 squeue      # Show what is running

=head1 DESCRIPTION

This program supports simple commands to set key elements of the NHLBI database.
The queue of tasks is kept in a MySQL database.
See B<perldoc DBIx::Connector> for details defining the database.

Functions wherefile, wherepath and whathost were moved into topmedpath.pl.

=head1 OPTIONS

=over 4

=item B<-bamid  -with-id>

Include the bamid or runid for the output from shown.

=item B<-emsg string>

Prints B<string> to STDOUT and saves it in the database.

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

=item B<-with-bamid>

Specifies that B<list samples RUNNAME> should also provide the bamid for the sample.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program are grouped into several groups which are used
to deal with specific sets of information in the monitor databases.

B<export>
Use this to create a CSV file of database columns for Chris.

B<mark bamid|nwdid dirname  [verb] [state]>
Use this to set the state for a particular BAM file.
You may specify the bamid or the NWDID.
Mark will set a date for the process (e.g. arrive sets state_arrive)
and unmark will set that entry to NULL.
The list of verbs and states can be seen by B<perldoc topmedcmd.pl>.

B<set bamid|nwdid|dirname columnname value>
Use this to set the value for a column for a particular BAM file
or run.

B<send2ncbi filelist>
Use this to copy data to NCBI with ascp.

B<show bamid|nwdid|dirname colname|center|run|yaml>
Use this to show information about a particular bamid (or expt_sampleid)
or run name.
Use 'yaml' to display everything known about the bam of interest.

B<list centers|runs|samples  value>
Use this to show a list of all centers, runs for a center or NWDIDs for a run.
The option B<-with-bamid> will cause the bamid to be shown in the last case.

B<unmark bamid|nwdid [verb]>
Use this to reset the state for a particular BAM file to the default
database value.

B<whatbamid bamname>
Use this to show the bamid for a sample using the name of the BAM or CRAM.

B<whatnwdid bamid|nwdid>
Use this to get some details for a particular bam.


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

