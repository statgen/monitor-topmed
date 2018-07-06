#!/usr/bin/perl
###################################################################
#
# Name: topmedcmd.pl
#
# Description:
#   Use this program to update the NHLBI database.
#   This program can work with topmed and inpsyght
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
  qq(/usr/cluster/topmed/lib/perl5),            # So we can run from /tmp
);

use Getopt::Long;
use Cwd qw(realpath abs_path);
use POSIX qw(strftime);
use My_DB;

my $NOTSET = 0;                     # Not set
my $REQUESTED = 1;                  # Task requested
my $SUBMITTED = 2;                  # Task submitted to be run
my $STARTED   = 3;                  # Task started
my $DELIVERED = 19;                 # Data delivered, but not confirmed
my $COMPLETED = 20;                 # Task completed successfully
my $CANCELLED = 89;                 # Task cancelled
my $FAILED    = 99;                 # Task failed

#   How to add a column to this program
#       Define the column in the database
#       Maybe create QOS for this new type of action
#       Add verbname => db column name to %VALIDVERBS
#       Add operation verb to %VALIDOPS in topmedpermit.pl
#       Add operation verb to statecols in topmed_failures.pl
#       Create new shell script for operation (e.g. topmed_cram.sh)
#       Add code for operation to topmed_monitor.pl
#       Add operation,letter to $quickcols,quickletter,validfunctions,TOPMEDJOBNAMES arrays
#           check letter in QuickStatus, update $STATUSLETTERS in index.php
#       Add operation verb to attributes2letter in topmed_status.pl
my %VALIDVERBS = (                  # Valid verbs to database colum
    arrive     => 'state_arrive',
    verify     => 'state_verify',
    gcebackup  => 'state_gcebackup',
    qplot      => 'state_qplot',
    cram       => 'state_cram',   
    gcepush    => 'state_gce38push',
    gcepull    => 'state_gce38pull',
    b37        => 'state_b37',   
    b38        => 'state_b38',     
    bcf        => 'state_gce38bcf',     
    gcecopy    => 'state_gce38copy',     
    gcecpbcf   => 'state_gce38cpbcf',     
    awscopy    => 'state_aws38copy',     
    ncbiexpt   => 'state_ncbiexpt',
    ncbiorig   => 'state_ncbiorig',
    ncbib37    => 'state_ncbib37',
    fix        => 'state_fix',
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
    bamfiles_table => 'bamfiles',
    centers_table => 'centers',
    permissions_table => 'permissions',
    runs_table => 'runs',
    ascpdest => 'asp-um-sph@gap-submit.ncbi.nlm.nih.gov:protected',    
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose with-id|bamid emsg=s project=s
    )) || die "$Script - Failed to parse options\n";
#   Non-typical code for PROJECT so non-project IDs can easily use this
if ($opts{project}) { $ENV{PROJECT} = $opts{project}; }
else {
    if (! exists($ENV{PROJECT})) { $ENV{PROJECT} = 'noproject'; }
}
if (! -d "/usr/cluster/$ENV{PROJECT}") {
   die "$Script - Environment variable PROJECT '$ENV{PROJECT}' incorrect\n";
}
$opts{realm} = "/usr/cluster/$ENV{PROJECT}/etc/.db_connections/$ENV{PROJECT}";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [-p project] [options] ";
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
        "$m show bamid|nwdid colname|run|center\n" .
        "  or\n" .
        "$m list centers\n" .
        "$m list runs centername\n" .
        "$m [-with-id] list samples runname\n" .
        "  or\n" .
        "$m send2ncbi files\n" .
        "  or\n" .
        "$m whatnwdid bamid|nwdid\n" .
        "  or\n" .
        "$m whatrun bamid|nwdid\n" .
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

DBConnect($opts{realm});

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'mark')      { Mark(@ARGV); exit; }
if ($fcn eq 'unmark')    { UnMark(@ARGV); exit; }
if ($fcn eq 'set')       { Set(@ARGV); exit; }
if ($fcn eq 'setdate')   { SetDate(@ARGV); exit; }
if ($fcn eq 'show')      { Show(@ARGV); exit; }
if ($fcn eq 'list')      { List(@ARGV); exit; }
if ($fcn eq 'send2ncbi') { Send2NCBI(@ARGV); exit; }
if ($fcn eq 'whatnwdid') { WhatNWDID(@ARGV); exit; }
if ($fcn eq 'whatrun')   { WhatRun(@ARGV); exit; }

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
    my $sth = DoSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    #   Set state for the verb
    my $col = $VALIDVERBS{$op};
    my $done = 0;
    if ($state eq 'requested') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col=$REQUESTED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'completed') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col=$COMPLETED WHERE bamid=$bamid");
        if ($col eq 'state_arrive') {       # Used by Kevin for tracking samples
            DoSQL("UPDATE $opts{bamfiles_table} SET datearrived='" . time() . "' WHERE bamid=$bamid");
        }
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            DoSQL("UPDATE $opts{bamfiles_table} SET datemapping_b37='" . time() . "' WHERE bamid=$bamid");
        }
        $done++;
    }
    if ($state eq 'delivered') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col=$DELIVERED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'started') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col=$STARTED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'failed') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col=$FAILED WHERE bamid=$bamid");
        if ($col eq 'state_arrive') {       # hack for Chris until new code in place
            DoSQL("UPDATE $opts{bamfiles_table} SET datearrived='-1' WHERE bamid=$bamid");
        }
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            DoSQL("UPDATE $opts{bamfiles_table} SET datemapping_b37='-1' WHERE bamid=$bamid");
        }
        $done++;
    }
    if ($state eq 'cancelled') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col=$CANCELLED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'submitted') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col=$SUBMITTED WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'notset') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col=$NOTSET WHERE bamid=$bamid");
        $done++;
    }
    if ($done) {
        if ($opts{verbose}) { print "$Script  'mark $bamid $op $state'  successful\n"; }
        if (exists($opts{emsg})) {
            print $opts{emsg} . "\n";
            DoSQL("UPDATE $opts{bamfiles_table} SET emsg=\"$op $state -- $opts{emsg}\" WHERE bamid=$bamid");
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
    my $sth = DoSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    my $col = $VALIDVERBS{$op};
    DoSQL("UPDATE $opts{bamfiles_table} SET $col=$NOTSET WHERE bamid=$bamid");
    if ($opts{verbose}) { print "$Script  'unmark $bamid $op'  successful\n"; }
}

#==================================================================
# Subroutine:
#   WhatNWDID($nwdid)
#
#   Print interesting details about an NWDID
#==================================================================
sub WhatNWDID {
    my ($nwdid) = @_;
    my $orignwdid = $nwdid;

    if ($nwdid =~ /^\d+/) {             # If bamid, get NWDID
        my $sth = DoSQL("SELECT expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$nwdid");
        if ($sth) {
            my $href = $sth->fetchrow_hashref;
            $nwdid = $href->{expt_sampleid};
        }
        else { die "$Script - NWDID '$nwdid' is unknown\n"; } 
    }
    else {                              # Extrace NWD from whatever was provided
        if ($nwdid =~ /(nwd\d+)/i) { $nwdid = uc($1); }
    }
    if (! $nwdid) { die "$Script - NWDID/BAMID '$orignwdid' is unknown\n"; }

    #   Reconstruct partial path to BAM
    my $sth = DoSQL("SELECT runid,bamid,piname,datayear FROM $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - NWDID '$nwdid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $bamid = $href->{bamid};
    my $datayear = $href->{datayear};
    my $piname = $href->{piname};
    $sth = DoSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$href->{runid}");
    $href = $sth->fetchrow_hashref;
    my $run = $href->{dirname};
    $sth = DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
    $href = $sth->fetchrow_hashref;
    my $center = uc($href->{centername});
    print "$nwdid $bamid can be found in run $run PI $piname for center $center year $datayear\n";
}

#==================================================================
# Subroutine:
#   WhatRun($bamid)
#
#   Print interesting details about an NWDID
#==================================================================
sub WhatRun {
    my ($bamid) = @_;
    my $sth;

    $bamid = GetBamid($bamid);
    $sth = DoSQL("SELECT runid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    if (! $sth) { die "$Script - Unknown '$bamid'\n"; }
    my $href = $sth->fetchrow_hashref;
    $sth = DoSQL("SELECT centerid,runid,dirname,count,datayear FROM $opts{runs_table} WHERE runid=$href->{runid}");
    $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};
    my $dirname = $href->{dirname};
    my $count = $href->{count};
    $sth = DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
    $href = $sth->fetchrow_hashref;
    my $center = uc($href->{centername});
    print "$bamid in center $center from run $dirname $runid which has $count samples\n";
}

#==================================================================
# Subroutine:
#   Send2NCBI($files)
#
#   Send files to NCBI
#==================================================================
sub Send2NCBI {
    my ($files) = @_;
    my $ascpcmd = "/usr/cluster/bin/ascp -i " .
        "/net/$ENV{PROJECT}/incoming/study.reference/send2ncbi/$ENV{PROJECT}-2-ncbi.pri " .
        "-l 800M -k 1";

    my $cmd = "$ascpcmd $files $opts{ascpdest}";
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
#   Set($bamid, $col, $val)
#
#   Set a database column
#==================================================================
sub Set {
    my ($bamid, $col, $val) = @_;
    my ($sth, $rowsofdata, $href);

    if (! defined($val)) { $val = ''; }
    chomp($val);                                # No trailing newlines please
    if ($val ne 'NULL') { $val = "'$val'"; }    # Use quotes unless this is NULL

    #   This could be a run name
    $sth = DoSQL("SELECT runid FROM $opts{runs_table} WHERE dirname='$bamid'", 0);
    $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        if ($rowsofdata > 1) {
            die "$Script - Eeek, there are $rowsofdata runs named '$bamid'\n";
        }
        $href = $sth->fetchrow_hashref;
        DoSQL("UPDATE $opts{runs_table} SET $col=$val WHERE runid='$href->{runid}'");
        return;
    }

    #   This is bamid or nwdid
    $bamid = GetBamid($bamid);

    #   Make sure this is a bam we know
    $sth = DoSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    if ($col eq 'nwdid') { $col = 'expt_sampleid'; }
    DoSQL("UPDATE $opts{bamfiles_table} SET $col=$val WHERE bamid=$bamid");
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
    $sth = DoSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    my @s = stat($val);
    if (! @s) { die "$Script - '$val' is not a known filename\n"; }
    my $datetime = strftime('%Y-%m-%d %H:%M:%S', localtime($s[10]));
    DoSQL("UPDATE $opts{bamfiles_table} SET $col='$datetime' WHERE bamid=$bamid");
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
        my $sth = DoSQL("SELECT centerid,centername FROM $opts{centers_table}");
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
        my $sth = DoSQL("SELECT centerid FROM $opts{centers_table} WHERE centername='$item'",0);
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - Unknown centername '$item'\n"; }
        my $href = $sth->fetchrow_hashref;
        my $centerid = $href->{centerid};
        $sth = DoSQL("SELECT runid,dirname FROM $opts{runs_table} WHERE centerid=$centerid");
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
            $sth = DoSQL("SELECT runid FROM $opts{runs_table} WHERE runid=$item", 0);
        }
        else {
            $sth = DoSQL("SELECT runid FROM $opts{runs_table} WHERE dirname='$item'", 0);
        }
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - Unknown run '$item'\n"; }
        my $href = $sth->fetchrow_hashref;
        my $runid = $href->{runid};
        $sth = DoSQL("SELECT bamid,expt_sampleid FROM $opts{bamfiles_table} WHERE runid=$runid");
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

    #   This could be a run name. Does not start with NWD, not all digits
    my $b = GetBamid($bamid, 1);    # Do not die, come back with zero if failed
    if (! $b) {                     # Not a sample, must be a run
        $sth = DoSQL("SELECT runid,$col FROM $opts{runs_table} WHERE dirname='$bamid'", 0);
        if ($sth) {
            $rowsofdata = $sth->rows();
            if ($rowsofdata) {
                $href = $sth->fetchrow_hashref;
                print $href->{$col} . "\n";
                return;
            }
        }
        die "$Script - Unknown '$bamid', not a sampleid, bamid or dir of a run\n";
    }
    else { $bamid = $b; }               #   This is bamid or nwdid

    $sth = DoSQL("SELECT bamid,runid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};

    #   Get run if asked for it
    if ($col eq 'run') {
        $sth = DoSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE runid=$href->{runid}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - RUNID '$href->runid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        print $href->{dirname} . "\n";
        return;
    }

    #   Get center if asked for it
    if ($col eq 'center') {
        $sth = DoSQL("SELECT centerid FROM $opts{runs_table} WHERE runid=$href->{runid}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - RUNID '$href->runid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;

        $sth = DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - CENTERID '$href->centerid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        print $href->{centername} . "\n";
        return;
    }

    #   Get value of column we asked for
    $sth = DoSQL("SELECT $col FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' or column '$col' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    if (defined($href->{$col})) { print $href->{$col} . "\n"; }
}

#==================================================================
# Subroutine:
#   GetBamid($bamid, $flag)
#
#   If this is not a sampleid or bamid and FLAG is 1, do not die, return zero
#
#   Return bamid for bamid or expt_sampleid
#==================================================================
sub GetBamid {
    my ($bamid, $flag) = @_;
    my ($sth, $rowsofdata, $href);
    if ($bamid =~ /^\d+$/) { return $bamid; }

    $sth = DoSQL("SELECT bamid FROM $opts{bamfiles_table} WHERE expt_sampleid='$bamid'");
    $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        $href = $sth->fetchrow_hashref;
        return $href->{bamid};
    }
    if ($flag) { return 0; }        # Not bamid or sample and do not die
    die "$Script - Invalid bamid or NWDID ($bamid). Try '$Script -help'\n";
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
    my $sth = DoSQL("SELECT runid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};
    $sth = DoSQL("SELECT centerid FROM $opts{runs_table} WHERE runid=$runid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' has no center?  How'd that happen?\n"; }
    $href = $sth->fetchrow_hashref;
    return ($href->{centerid}, $runid);
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

  topmedcmd.pl -proj MYPROJ show 2199 state_cram        # Show a column
  topmedcmd.pl show 2199 center            # Show center 
  topmedcmd.pl show NWD00234 run           # Show run
  topmedcmd.pl show 2016apr20 offsite      # Show offsite for run

  topmedcmd.pl showrun 20150604            # Show all NWDID for run
  topmedcmd.pl -bamid showrun 20150604     # Show all NWDID and bamid for run
 
  topmedcmd.pl list centers                # Show all known centers
  topmedcmd.pl list runs broad             # Show all runs for center broad
  topmedcmd.pl list samples 2015dec02      # Show all samples for a run
  
  topmedcmd.pl send2ncbi files             # Send files to NCBI
 
  topmedcmd.pl whatnwdid bamid|nwdid       # Show details for a sample
 
  topmedcmd.pl whatrun bamid|nwdid         # Show details of the run for a sample
 
=head1 DESCRIPTION

This program supports simple commands to set key elements of the NHLBI database.
The queue of tasks is kept in a MySQL database.

This program requires you to provide the project name (lower case) as either
environment variable or with the option B<-project>.
Projects must have a directory as /usr/cluster/projectname.

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

=item B<-project projectname>

Allows one to specify the project to work on. You must specify this option
or set the environment variable B<PROJECT}.

=item B<-with-bamid>

Specifies that B<list samples RUNNAME> should also provide the bamid for the sample.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program are grouped into several groups which are used
to deal with specific sets of information in the monitor databases.

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

B<show bamid|nwdid|dirname colname|center|run>
Use this to show information about a particular bamid (or expt_sampleid)
or run name.

B<list centers|runs|samples  value>
Use this to show a list of all centers, runs for a center or NWDIDs for a run.
The option B<-with-bamid> will cause the bamid to be shown in the last case.

B<unmark bamid|nwdid [verb]>
Use this to reset the state for a particular BAM file to the default
database value.

B<whatnwdid bamid|nwdid>
Use this to get some details for a particular bam.

B<whatrun bamid|nwdid>
Use this to get some details of a run for a particular bam.


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

