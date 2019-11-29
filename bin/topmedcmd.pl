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
    backup     => 'state_backup',
    qplot      => 'state_qplot',
    cram       => 'state_cram',   
    gcepush    => 'state_gce38push',
    gcepull    => 'state_gce38pull',
    b37        => 'state_b37',   
    b38        => 'state_b38',     
    bcf        => 'state_gce38bcf',     
    gcecopy    => 'state_gce38copy',     
    gcecpbcf   => 'state_gce38cpbcf',     
    gcecleanup => 'state_gcecleanup',     
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
    runs_table => 'runs',
    runs_pkey => 'runid',
    samples_table => 'bamfiles',
    samples_pkey => 'bamid',
    files_table => 'na',
    files_pkey => 'na',
    datatype => 'genome',
    centers_table => 'centers',
    permissions_table => 'permissions',
    ascpdest => 'asp-um-sph@gap-submit.ncbi.nlm.nih.gov:protected',    
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose with-id|bamid emsg=s project=s datatype=s
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
if ($opts{datatype} eq 'rnaseq') {
	$opts{samples_table} = 'tx_samples';
	$opts{samples_pkey} = 'txseqid';
	$opts{runs_table} = 'tx_projects';
	$opts{runs_pkey} = 'rnaprojectid';
	$opts{files_table} = 'tx_files';
	$opts{files_pkey} = 'fileid';
}

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [-p project] [options] ";
    my $verbs = join(' ', sort keys %VALIDVERBS);
    my $requests = join(' ', sort keys %VALIDSTATUS);
    warn "$m mark sampleid|nwdid action newstatus\n" .
        "    Where action: $verbs\n" .
        "    Where newstatus: $requests\n" .
        "  or\n" .
        "$m unmark sampleid|nwdid [same list as mark]\n" .
        "  or\n" .
        "$m set sampleid|nwdid|dirname colname value\n" .
        "  or\n" .
        "$m setdate sampleid|nwdid colname file\n" .
        "  or\n" .
        "$m show sampleid|nwdid colname|run|center\n" .
        "  or\n" .
        "$m list centers\n" .
        "$m list runs centername\n" .
        "$m [-with-id] list samples runname\n" .
        "  or\n" .
        "$m whatnwdid sampleid|nwdid\n" .
        "  or\n" .
        "$m whatrun sampleid|nwdid\n" .
        "  or\n" .
        "$m permit add operation center run\n" .
        "$m permit remove permitid\n" .
        "$m permit test operation sampleid/n" .
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
if ($fcn eq 'whatnwdid') { WhatNWDID(@ARGV); exit; }
if ($fcn eq 'whatrun')   { WhatRun(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   Mark($sampleid, $op, $state)
#
#   Set states in the database
#==================================================================
sub Mark {
    my ($sampleid, $op, $state) = @_;
    $sampleid = GetSampleid($sampleid);
    if ((! exists($VALIDVERBS{$op})) || (! exists($VALIDSTATUS{$state}))) {
        die "$Script - Invalid 'mark' syntax. Try '$Script -help'\n";
    }

    #   Make sure this is a sample we know
    my $sth = DoSQL("SELECT $opts{samples_pkey} FROM $opts{samples_table} WHERE $opts{samples_pkey}=$sampleid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Sample '$sampleid' is unknown\n"; }

    #   Set state for the verb
    my $col = $VALIDVERBS{$op};
    my $done = 0;
    if ($state eq 'requested') {
        DoSQL("UPDATE $opts{samples_table} SET $col=$REQUESTED WHERE $opts{samples_pkey}=$sampleid");
        $done++;
    }
    if ($state eq 'completed') {
        DoSQL("UPDATE $opts{samples_table} SET $col=$COMPLETED WHERE $opts{samples_pkey}=$sampleid");
        if ($col eq 'state_arrive') {       # Used by Kevin for tracking samples
            DoSQL("UPDATE $opts{samples_table} SET datearrived='" . time() . "' WHERE $opts{samples_pkey}=$sampleid");
        }
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            DoSQL("UPDATE $opts{samples_table} SET datemapping_b37='" . time() . "' WHERE $opts{samples_pkey}=$sampleid");
        }
        $done++;
    }
    if ($state eq 'delivered') {
        DoSQL("UPDATE $opts{samples_table} SET $col=$DELIVERED WHERE $opts{samples_pkey}=$sampleid");
        $done++;
    }
    if ($state eq 'started') {
        DoSQL("UPDATE $opts{samples_table} SET $col=$STARTED WHERE $opts{samples_pkey}=$sampleid");
        $done++;
    }
    if ($state eq 'failed') {
        DoSQL("UPDATE $opts{samples_table} SET $col=$FAILED WHERE $opts{samples_pkey}=$sampleid");
        if ($col eq 'state_arrive') {       # hack for Chris until new code in place
            DoSQL("UPDATE $opts{samples_table} SET datearrived='-1' WHERE $opts{samples_pkey}=$sampleid");
        }
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            DoSQL("UPDATE $opts{samples_table} SET datemapping_b37='-1' WHERE $opts{samples_pkey}=$sampleid");
        }
        $done++;
    }
    if ($state eq 'cancelled') {
        DoSQL("UPDATE $opts{samples_table} SET $col=$CANCELLED WHERE $opts{samples_pkey}=$sampleid");
        $done++;
    }
    if ($state eq 'submitted') {
        DoSQL("UPDATE $opts{samples_table} SET $col=$SUBMITTED WHERE $opts{samples_pkey}=$sampleid");
        $done++;
    }
    if ($state eq 'notset') {
        DoSQL("UPDATE $opts{samples_table} SET $col=$NOTSET WHERE $opts{samples_pkey}=$sampleid");
        $done++;
    }
    if ($done) {
        if ($opts{verbose}) { print "$Script  'mark $sampleid $op $state'  successful\n"; }
        if (exists($opts{emsg})) {
            print $opts{emsg} . "\n";
            DoSQL("UPDATE $opts{samples_table} SET emsg=\"$op $state -- $opts{emsg}\" WHERE $opts{samples_pkey}=$sampleid");
        }
    }
    else { die "$Script - Invalid state '$state' for '$op'. Try '$Script -help'\n"; }
}

#==================================================================
# Subroutine:
#   UnMark($dirname, $op)
#
#   Reset state in the database
#==================================================================
sub UnMark {
    my ($sampleid, $op) = @_;
    $sampleid = GetSampleid($sampleid);
    if (! exists($VALIDVERBS{$op})) {
        die "$Script - Invalid 'unmark' syntax. Try '$Script -help'\n";
    }

    #   Make sure this is a sample we know
    my $sth = DoSQL("SELECT $opts{samples_pkey} FROM $opts{samples_table} WHERE $opts{samples_pkey}=$sampleid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Sample '$sampleid' is unknown\n"; }

    my $col = $VALIDVERBS{$op};
    DoSQL("UPDATE $opts{samples_table} SET $col=$NOTSET WHERE $opts{samples_pkey}=$sampleid");
    if ($opts{verbose}) { print "$Script  'unmark $sampleid $op'  successful\n"; }
}

#==================================================================
# Subroutine:
#   WhatNWDID($nwdid)
#
#   Print interesting details about an NWDID
#==================================================================
sub WhatNWDID {
    my ($nwdid) = @_;
    my $sampleid = GetSampleid($nwdid);

    #   Reconstruct partial path to sample
    my $sth = DoSQL("SELECT $opts{runs_pkey},piname,expt_sampleid,datayear FROM $opts{samples_table} WHERE $opts{samples_pkey}='$sampleid'");
    my $rowsofdata = $sth->rows();
    my $href = $sth->fetchrow_hashref;
    my $datayear = $href->{datayear};
    my $piname = $href->{piname};
    $nwdid = $href->{expt_sampleid};
    $sth = DoSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE $opts{runs_pkey}=$href->{$opts{runs_pkey}}");
    $href = $sth->fetchrow_hashref;
    my $run = $href->{dirname};
    $sth = DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
    $href = $sth->fetchrow_hashref;
    my $center = uc($href->{centername});
    print "$nwdid $sampleid can be found in run $run PI $piname for center $center year $datayear\n";
}

#==================================================================
# Subroutine:
#   WhatRun($sampleid)
#
#   Print interesting details about an NWDID
#==================================================================
sub WhatRun {
    my ($sampleid) = @_;
    my $sth;

    $sampleid = GetSampleid($sampleid);
    $sth = DoSQL("SELECT $opts{runs_pkey} FROM $opts{samples_table} WHERE $opts{samples_pkey}=$sampleid");
    if (! $sth) { die "$Script - Unknown '$sampleid'\n"; }
    my $href = $sth->fetchrow_hashref;
    $sth = DoSQL("SELECT centerid,$opts{runs_pkey},dirname,count,datayear FROM $opts{runs_table} WHERE $opts{runs_pkey}=$href->{$opts{runs_pkey}}");
    $href = $sth->fetchrow_hashref;
    my $id = $href->{$opts{runs_pkey}};
    my $dirname = $href->{dirname};
    my $count = $href->{count};
    $sth = DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
    $href = $sth->fetchrow_hashref;
    my $center = uc($href->{centername});
    print "$sampleid in center $center from run $dirname $id which has $count samples\n";
}

#==================================================================
# Subroutine:
#   Set($sampleid, $col, $val)
#
#   Set a database column
#==================================================================
sub Set {
    my ($sampleid, $col, $val) = @_;
    my ($sth, $rowsofdata, $href);

    if (! defined($val)) { $val = ''; }
    chomp($val);                                # No trailing newlines please
    if ($val ne 'NULL') { $val = "'$val'"; }    # Use quotes unless this is NULL

    #   This could be a run name
    $sth = DoSQL("SELECT $opts{runs_pkey} FROM $opts{runs_table} WHERE dirname='$sampleid'", 0);
    $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        if ($rowsofdata > 1) {
            die "$Script - Eeek, there are $rowsofdata runs named '$sampleid'\n";
        }
        $href = $sth->fetchrow_hashref;
        DoSQL("UPDATE $opts{runs_table} SET $col=$val WHERE $opts{runs_pkey}='$href->{$opts{runs_pkey}}'");
        return;
    }

    #   This is sampleid or expt_sampleid
    $sampleid = GetSampleid($sampleid);

    #   Make sure this is a sample we know
    $sth = DoSQL("SELECT $opts{samples_pkey} FROM $opts{samples_table} WHERE $opts{samples_pkey}=$sampleid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Sample '$sampleid' is unknown\n"; }

    if ($col eq 'nwdid') { $col = 'expt_sampleid'; }
    DoSQL("UPDATE $opts{samples_table} SET $col=$val WHERE $opts{samples_pkey}=$sampleid");
}

#==================================================================
# Subroutine:
#   SetDate($sampleid, $col, $val)
#
#   Set a database column to the date for a file
#==================================================================
sub SetDate {
    my ($sampleid, $col, $val) = @_;
    my ($sth, $rowsofdata, $href);

    #   This is sampleid or expt_sampleid
    $sampleid = GetSampleid($sampleid);

    #   Make sure this is a sample we know
    $sth = DoSQL("SELECT $opts{samples_pkey} FROM $opts{samples_table} WHERE $opts{samples_pkey}=$sampleid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Sample '$sampleid' is unknown\n"; }

    my @s = stat($val);
    if (! @s) { die "$Script - '$val' is not a known filename\n"; }
    my $datetime = strftime('%Y-%m-%d %H:%M:%S', localtime($s[10]));
    DoSQL("UPDATE $opts{samples_table} SET $col='$datetime' WHERE $opts{samples_pkey}=$sampleid");
}

#==================================================================
# Subroutine:
#   List($fcn, $item)
#
#   Generate a list of data from the database
#   Where $fcn may be centers, runs , samples or files for a sample
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
        $sth = DoSQL("SELECT $opts{runs_pkey},dirname FROM $opts{runs_table} WHERE centerid=$centerid");
        $rowsofdata = $sth->rows();
        for (my $i=1; $i<=$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref;
            $s = $href->{dirname};
            if ($opts{'with-id'}) { $s .= ' ' . $href->{$opts{runs_pkey}}; }
            print $s . "\n";
        }
        return;
    }
    if ($fcn eq 'samples') {            # Show all samples for a run
        if (! $item) { die "$Script - Runname was not provided\n"; }
        my $sth;
        if ($item =~ /^\d+$/) {        # Maybe runid provided
            $sth = DoSQL("SELECT $opts{runs_pkey} FROM $opts{runs_table} WHERE $opts{runs_pkey}=$item", 0);
        }
        else {
            $sth = DoSQL("SELECT $opts{runs_pkey} FROM $opts{runs_table} WHERE dirname='$item'", 0);
        }
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - Unknown run '$item'\n"; }
        my $href = $sth->fetchrow_hashref;
        my $id = $href->{$opts{runs_pkey}};
        $sth = DoSQL("SELECT $opts{samples_pkey},expt_sampleid FROM $opts{samples_table} WHERE $opts{runs_pkey}=$id");
        $rowsofdata = $sth->rows();
        for (my $i=1; $i<=$rowsofdata; $i++) {
            $href = $sth->fetchrow_hashref;
            $s = $href->{expt_sampleid};
            if ($opts{'with-id'}) { $s .= ' ' . $href->{$opts{samples_pkey}}; }
            print $s . "\n";
        }
        return;
    }
    if ($fcn eq 'files') {            	# Show all files for a sample
        if (! $item) { die "$Script - Sample was not provided\n"; }
       	my $sampleid = GetSampleid($item);
        my $sth = DoSQL("SELECT $opts{samples_pkey},filename,checksum FROM $opts{files_table} WHERE $opts{samples_pkey}=$sampleid");
        my $rowsofdata = $sth->rows();
        for (my $i=1; $i<=$rowsofdata; $i++) {
            my $href = $sth->fetchrow_hashref;
            my $s = '';
            if ($opts{'with-id'}) { $s .= $href->{$opts{samples_pkey}} . ' '; }
            $s = $href->{filename} . ' ' . $href->{checksum};
            print $s . "\n";
        }
        return;
    }
    die "$Script - Unknown list function '$fcn'\n";
}

#==================================================================
# Subroutine:
#   Show($fcn, $sampleid, $col)
#
#   Generate list of information from the database
#   or show the value for a column
#==================================================================
sub Show {
    my ($sampleid, $col) = @_;
    my ($sth, $rowsofdata, $href);

    #   This could be a run name. Does not start with NWD, not all digits
    my $b = GetSampleid($sampleid, 1);    # Do not die, come back with zero if failed
    if (! $b) {                     # Not a sample, must be a run
        $sth = DoSQL("SELECT $opts{runs_pkey},$col FROM $opts{runs_table} WHERE dirname='$sampleid'", 0);
        if ($sth) {
            $rowsofdata = $sth->rows();
            if ($rowsofdata) {
                $href = $sth->fetchrow_hashref;
                print $href->{$col} . "\n";
                return;
            }
        }
        die "$Script - Unknown '$sampleid', not a sampleid, id or dir of a run\n";
    }
    else { $sampleid = $b; }               #   This is sampleid or nwd|tor

    $sth = DoSQL("SELECT $opts{samples_pkey},$opts{runs_pkey} FROM $opts{samples_table} WHERE $opts{samples_pkey}=$sampleid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Sample '$sampleid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $id = $href->{$opts{runs_pkey}};

    #   Get run if asked for it
    if ($col eq 'run') {
        $sth = DoSQL("SELECT centerid,dirname FROM $opts{runs_table} WHERE $opts{runs_pkey}=$href->{$opts{runs_pkey}}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - RUNID '$href->$opts{runs_pkey}' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        print $href->{dirname} . "\n";
        return;
    }

    #   Get center if asked for it
    if ($col eq 'center') {
        $sth = DoSQL("SELECT centerid FROM $opts{runs_table} WHERE $opts{runs_pkey}=$href->{$opts{runs_pkey}}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - RUNID '$href->$opts{runs_pkey}' is unknown\n"; }
        $href = $sth->fetchrow_hashref;

        $sth = DoSQL("SELECT centername FROM $opts{centers_table} WHERE centerid=$href->{centerid}");
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - CENTERID '$href->centerid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        print $href->{centername} . "\n";
        return;
    }

    #   Get value of column we asked for
    $sth = DoSQL("SELECT $col FROM $opts{samples_table} WHERE $opts{samples_pkey}=$sampleid");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Sample '$sampleid' or column '$col' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    if (defined($href->{$col})) { print $href->{$col} . "\n"; }
}

#==================================================================
# Subroutine:
#   GetSampleid($sampleid, $flag)
#
#   If this is not a sampleid or id and FLAG is 1, do not die, return zero
#
#   Return pkeyid for expt_sampleid
#==================================================================
sub GetSampleid {
    my ($sampleid, $flag) = @_;
    my ($sth, $rowsofdata, $href);
    if ($sampleid =~ /^\d+$/) { return $sampleid; }

    $sth = DoSQL("SELECT $opts{samples_pkey} FROM $opts{samples_table} WHERE expt_sampleid='$sampleid'");
    $rowsofdata = $sth->rows();
    if ($rowsofdata) {
        $href = $sth->fetchrow_hashref;
        return $href->{$opts{samples_pkey}};
    }
    if ($flag) { return 0; }        # Not sampleid or sample and do not die
    die "$Script - Invalid sampleid or NWD/TOR-ID ($sampleid). Try '$Script -help'\n";
}

#==================================================================
# Subroutine:
#   ($centerid, $id) = GetSampleidInfo($sampleid)
#
#   Return id for the center and run_id for a sample
#==================================================================
sub GetSampleidInfo {
    my ($sampleid) = @_;

    if ($sampleid !~ /^\d+$/) { return (0,0); }     # No sampleid, no ids
    my $sth = DoSQL("SELECT $opts{runs_pkey} FROM $opts{samples_table} WHERE $opts{samples_pkey}=$sampleid");
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Sample '$sampleid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $id = $href->{$opts{runs_pkey}};
    $sth = DoSQL("SELECT centerid FROM $opts{runs_table} WHERE $opts{runs_pkey}=$id");
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - Sample '$sampleid' has no center?  How'd that happen?\n"; }
    $href = $sth->fetchrow_hashref;
    return ($href->{centerid}, $id);
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedcmd.pl - Update the database for NHLBI TopMed

=head1 SYNOPSIS

  topmedcmd.pl mark 33 arrive completed   # sample has arrived
  topmedcmd.pl -emsg 'No file found' mark 33 cram failed   # Action failed, set error msg
  topmedcmd.pl mark NWD00234  arrive completed   # Same sample has arrived
  topmedcmd.pl unmark 33 arrive           # Reset sample has arrived

  topmedcmd.pl set 33 jobidqplot 123445    # Set jobidqplot in database
  topmedcmd.pl set NWD123433 jobidqplot 123445    # Set jobidqplot in database
  topmedcmd.pl set 2016apr20 offsite N     # Set 2016apr20 in database

  topmedcmd.pl -proj MYPROJ show 2199 state_cram        # Show a column
  topmedcmd.pl show 2199 center            # Show center 
  topmedcmd.pl show NWD00234 run           # Show run
  topmedcmd.pl show 2016apr20 offsite      # Show offsite for run

  topmedcmd.pl showrun 20150604            # Show all NWDID for run
  topmedcmd.pl -bamid showrun 20150604     # Show all NWDID and bamid for run
 
  topmedcmd.pl list centers                # Show all known centers
  topmedcmd.pl list runs broad             # Show all runs for center broad
  topmedcmd.pl list samples 2015dec02      # Show all samples for a run
  
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

=item B<-datatype genome|rnaseq>

Specifies this is a particular datatype of data. The default is 'genome'.

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

