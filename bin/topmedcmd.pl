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
use Cwd qw(realpath);

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
    incomingdir => 'incoming/topmed',
    backupsdir => 'working/backups/incoming/topmed',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options]";
    my $verbs = join(',', sort keys %VALIDVERBS);
    my $requests = join(',', sort values %VALIDSTATUS);
    warn "$m mark bamid|NWDnnnnn $verbs $requests\n" .
        "  or\n" .
        "$m unmark bamid [same list as mark]\n" .
        "  or\n" .
        "$m set bamid colname value\n" .
        "  or\n" .
        "$m show arrived\n" .
        "  or\n" .
        "$m show bamid colname|run|center\n" .
        "  or\n" .
        "$m export\n" .
        "  or\n" .
        "$m where bamid\n" .
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

#   Our network can get flakey and the DBConnect can fail
#   Catch when this happens and wait a bit and try again
my $dbh;
my $sleeptime = 10;
for (1 .. 10) {
    eval { $dbh = DBConnect($opts{realm}); };
    if ($@) {                           # Failed, wait a bit and try again
        print "Datbase connection failed, wait and retry\n";
        sleep($sleeptime);
        $sleeptime += 10;
    }
    else { last; }
}
if ($@) { die $@ . "\n"; }
#my $dbh = DBConnect($opts{realm});

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'mark')     { Mark(@ARGV); exit; }
if ($fcn eq 'unmark')   { UnMark(@ARGV); exit; }
if ($fcn eq 'set')      { Set(@ARGV); exit; }
if ($fcn eq 'show')     { Show(@ARGV); exit; }
if ($fcn eq 'export')   { Export(@ARGV); exit; }
if ($fcn eq 'where')    { Where(@ARGV); exit; }
if ($fcn eq 'whatnwdid')  { WhatNWDID(@ARGV); exit; }
if ($fcn eq 'permit')   { Permit(@ARGV); exit; }

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
    if ($bamid =~ /NWD/) {              # Special Hack for Chris cause he does not same BAMID
        my $sth = DoSQL("SELECT bamid from $opts{bamfiles_table} WHERE expt_sampleid='$bamid'", 0);
        if ($sth) {
            my $href = $sth->fetchrow_hashref;
            $bamid = $href->{bamid};
        }
    }
    if (($bamid !~ /^\d+$/) || (! exists($VALIDVERBS{$op})) || (! exists($VALIDSTATUS{$state}))) {
        die "$Script Invalid 'mark' syntax. Try '$Script -help'\n";
    }

    #   Make sure this is a bam we know
    my $sth = DoSQL("SELECT bamid from $opts{bamfiles_table} WHERE bamid=$bamid", 0);
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
        if ($col eq 'state_arrive') {       # hack for Chris until new code in place
            DoSQL("UPDATE $opts{bamfiles_table} SET datearrived='" . time() . "' WHERE bamid=$bamid");
        }
        if ($col eq 'state_md5ver') {       # hack for Chris until new code in place
            DoSQL("UPDATE $opts{bamfiles_table} SET datemd5ver='" . time() . "' WHERE bamid=$bamid");
        }
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            DoSQL("UPDATE $opts{bamfiles_table} SET datemapping='" . time() . "' WHERE bamid=$bamid");
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
        if ($col eq 'state_md5ver') {       # hack for Chris until new code in place
            DoSQL("UPDATE $opts{bamfiles_table} SET datemd5ver='-1' WHERE bamid=$bamid");
        }
        if ($col eq 'state_b37') {          # hack for Chris until new code in place
            DoSQL("UPDATE $opts{bamfiles_table} SET datemapping='-1' WHERE bamid=$bamid");
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
    }
    else { die "$Script Invalid state '$state' for '$op'. Try '$Script -help'\n"; }
}

#==================================================================
# Subroutine:
#   UnMark($dirname, $op)
#
#   Reset state in the bamfiles database
#==================================================================
sub UnMark {
    my ($bamid, $op) = @_;
    if (($bamid !~ /^\d+$/) || (! exists($VALIDVERBS{$op}))) {
        die "$Script Invalid 'unmark' syntax. Try '$Script -help'\n";
    }

    #   Make sure this is a bam we know
    my $sth = DoSQL("SELECT bamid from $opts{bamfiles_table} WHERE bamid=$bamid", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    my $col = $VALIDVERBS{$op};
    DoSQL("UPDATE $opts{bamfiles_table} SET $col=$NOTSET WHERE bamid=$bamid");
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
    my @cols = qw(bamname studyname piname expt_sampleid);
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
            my $sth = DoSQL($sql);
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
#   WhatNWDID($nwdid)
#
#   Print interesting details about an NWDID
#==================================================================
sub WhatNWDID {
    my ($nwdid) = @_;

    if ($nwdid =~ /^\d+/) {             # If bamid, get NWDID
        my $sth = DoSQL("SELECT expt_sampleid from $opts{bamfiles_table} WHERE bamid=$nwdid", 0);
        if ($sth) {
            my $href = $sth->fetchrow_hashref;
            $nwdid = $href->{expt_sampleid};
        }
    }
  
    #   Reconstruct partial path to BAM
    my $sth = DoSQL("SELECT runid from $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - NWDID '$nwdid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    $sth = DoSQL("SELECT dirname from $opts{runs_table} WHERE runid=$href->{runid}");
    $href = $sth->fetchrow_hashref;
    print "'$nwdid' can be found in run '$href->{dirname}'\n";
}

#==================================================================
# Subroutine:
#   Where($bamid)
#
#   Print path to bamid, path to the backup directory, bamname, real host of bamname, numeric part of realhost
#==================================================================
sub Where {
    my ($bamid) = @_;

    #   Reconstruct partial path to BAM
    my $sth = DoSQL("SELECT runid,bamname from $opts{bamfiles_table} WHERE bamid=$bamid", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $bamname = $href->{bamname};
    my $runid = $href->{runid};
    $sth = DoSQL("SELECT centerid,dirname from $opts{runs_table} WHERE runid=$runid", 0);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' run '$runid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $rundir = $href->{dirname};
    my $centerid = $href->{centerid};
    $sth = DoSQL("SELECT centername from $opts{centers_table} WHERE centerid=$centerid", 0);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' center '$centerid' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    my $centername = $href->{centername};

    #   The BAM is in one of those $opts{netdir} trees where center is not a symlink
    my $bamfdir = '';
    foreach ('', '2', '3', '4') {
        $bamfdir = "$opts{netdir}$_/$opts{incomingdir}/$centername";
        if (! -l $bamfdir) { last; }        # Found non-symlink to center directory
    }
    if (! $bamfdir) { die "BAMID=$bamid Unable to find real directory for '$centername'\n"; }
    $bamfdir .= '/' . $rundir;

    #   Because we can't really plan very well, some files in a directory might actually
    #   live on some other host. All center directories and files should 'look' the
    #   same on all hosts, but some directories can have files scattered all over the place.
    my $realbamname = realpath("$bamfdir/$bamname");    # e.g. /net/topmed2/incoming/topmed/broad/2015jun22/xxxx
    my $realhost = 'unknown';
    if ($realbamname =~ /^($opts{netdir}\d*)\//)  { $realhost = substr($1,5); } # e.g. topmedN
    my $realhostindex = '';
    if ($realhost =~ /(\d)/)  { $realhostindex = $1; }     # e.g. 2

    #   The backup for the BAM is in another tree, also not a symlink
    my $bakbamfdir = $opts{netdir};
    if ($realhostindex eq '') { $bakbamfdir = $opts{netdir} . '2'; }
    if ($realhostindex eq '2') { $bakbamfdir = $opts{netdir} . ''; }
    if ($realhostindex eq '3') { $bakbamfdir = $opts{netdir} . '3'; }   # No idea what to do here
    if ($realhostindex eq '4') { $bakbamfdir = $opts{netdir} . 'x4'; }  # No idea what to do here
    $bakbamfdir = "$bakbamfdir/$opts{backupsdir}/$centername";
    if (-l $bakbamfdir) { die "BAMID=$bamid bamdir=$bamfdir Backup directory may not be a symlink for '$bakbamfdir'\n"; }
    $bakbamfdir .= '/' . $rundir;

    print "$bamfdir $bakbamfdir $bamname $realhost $realhostindex\n";
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
        $sth = DoSQL("SELECT runid,centerid,operation FROM $opts{permissions_table}", 0);
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
        $sth = DoSQL($sql, 0);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "Permission '$id' does not exist\n"; }
        $href = $sth->fetchrow_hashref;
        my $s = "$href->{centername} / $href->{dirname} / $href->{operation}";
        $sql = "DELETE FROM $opts{permissions_table} WHERE id='$id'";
        $sth = DoSQL($sql, 0);
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
        $sth = DoSQL($sql, 0);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "Center '$center' is not known\n"; }
        $href = $sth->fetchrow_hashref;
        $centerid = $href->{centerid};
    }

    if ($run ne 'all') {
        $sql = "SELECT runid FROM $opts{runs_table} WHERE dirname='$run'";
        $sth = DoSQL($sql, 0);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "Run '$run' is not known\n"; }
        $href = $sth->fetchrow_hashref;
        $runid = $href->{runid};
    }

    $sql = "INSERT INTO $opts{permissions_table} " .
        "(centername,dirname,centerid,runid,operation) " .
        "VALUES ('$center', '$run', $centerid, $runid, '$op')";
    $sth = DoSQL($sql, 0);
    if (! $sth) { die "Failed to add permission '$fcn' for '$op'\nSQL=$sql"; }
    print "Added permission '$center / $runid / $op'\n";
    exit;
}

#==================================================================
# Subroutine:
#   ($centerid, $runid) = GetBamidInfo($bamid)
#
#==================================================================
sub GetBamidInfo {
    my ($bamid) = @_;

    my $sth = DoSQL("SELECT runid FROM $opts{bamfiles_table} WHERE bamid=$bamid", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;
    my $runid = $href->{runid};
    $sth = DoSQL("SELECT centerid FROM $opts{runs_table} WHERE runid=$runid", 0);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' has no center?  How'd that happen?\n"; }
    $href = $sth->fetchrow_hashref;
    return ($href->{centerid}, $runid);
}

#==================================================================
# Subroutine:
#   Set($bamid, $col, $val)
#
#   Set a database column
#==================================================================
sub Set {
    my ($bamid, $col, $val) = @_;
    if ($bamid !~ /^\d+$/){
        die "$Script Invalid 'set' syntax. Try '$Script -help'\n";
    }

    #   Make sure this is a bam we know
    my $sth = DoSQL("SELECT bamid from $opts{bamfiles_table} WHERE bamid=$bamid", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }

    if ($col eq 'nwdid') { $col = 'expt_sampleid'; }
    DoSQL("UPDATE $opts{bamfiles_table} SET $col='$val' WHERE bamid=$bamid");
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
    if ($bamid eq 'arrived') { return ShowArrived($bamid); }

    if ($bamid !~ /^\d+$/){
        die "$Script Invalid 'show' syntax. Try '$Script -help'\n";
    }

    my $sth = DoSQL("SELECT bamid,runid from $opts{bamfiles_table} WHERE bamid=$bamid", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' is unknown\n"; }
    my $href = $sth->fetchrow_hashref;

    #   Get run if asked for it
    if ($col eq 'run') {
        $sth = DoSQL("SELECT centerid,dirname from $opts{runs_table} WHERE runid=$href->{runid}", 0);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - RUNID '$href->runid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        print $href->{dirname} . "\n";
        return;
    }

    #   Get center if asked for it
    if ($col eq 'center') {
        $sth = DoSQL("SELECT centerid from $opts{runs_table} WHERE runid=$href->{runid}", 0);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - RUNID '$href->runid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;

        $sth = DoSQL("SELECT centername from $opts{centers_table} WHERE centerid=$href->{centerid}", 0);
        $rowsofdata = $sth->rows();
        if (! $rowsofdata) { die "$Script - CENTERID '$href->centerid' is unknown\n"; }
        $href = $sth->fetchrow_hashref;
        print $href->{centername} . "\n";
        return;
    }

    #   Get value of column we asked for
    $sth = DoSQL("SELECT $col from $opts{bamfiles_table} WHERE bamid=$bamid", 0);
    $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$Script - BAM '$bamid' or column '$col' is unknown\n"; }
    $href = $sth->fetchrow_hashref;
    print $href->{$col} . "\n";
}

#==================================================================
# Subroutine:
#   ShowArrived($fcn)
#
#   Generate list of information from the database
#==================================================================
sub ShowArrived {
    my ($fcn) = @_;
    #   Get all the known centers in the database
    my $centersref = GetCenters();
    foreach my $cid (keys %{$centersref}) {
        my $centername = $centersref->{$cid};
        my $runsref = GetRuns($cid) || next;
        #   For each run, see if there are bamfiles that arrived
        foreach my $runid (keys %{$runsref}) {
            my $dirname = $runsref->{$runid};
            #   Get list of all bams that have not yet arrived properly
            my $sql = "SELECT bamid,bamname,state_arrived FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" .
                    $href->{bamname};
                #   See if this has arrived. Few states possible
                if ($href->{state_arrive} != $COMPLETED) { next; }
                #   Run the command
                print "$href->{bamid} $centername $dirname $f\n";
            }
        }
    }
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
 
  topmedcmd.pl where 2199                  # Returns path to bam, path to backup, bamname

  topmedcmd.pl permit add bai braod 2015oct18   # Stop bai job submissions for a run
  topmedcmd.pl permit remove 12             # Remove a permit control
  topmedcmd.pl permit test bai 4567         # Test if we should submit a bai job for one bam

=head1 DESCRIPTION

This program supports simple commands to set key elements of the NHLBI database.
The queue of tasks is kept in a MySQL database.
See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-runs NAME[,NAME,...]>

Specifies a specific set of runs on which to run the action,
e.g. B<2015jun05.weiss.02,2015jun05.weiss.03>.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

Parameters to this program are grouped into several groups which are used
to deal with specific sets of information in the monitor databases.

B<mark bamid dirname  [verb] [state]>
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

B<set bamid columnname value>
Use this to set the value for a column for a particular BAM file.

B<show arrived>
Use this to show the bamids for all BAMs that are marked arrived.

B<unmark bamid [verb]>
Use this to reset the state for a particular BAM file to the default
database value.

B<whatnwdid bamid|NWDnnnnnn>
Use this to get some details for a particular bam.

<where bamid>
Use this to display the directory of the bam file, 
the path to the backup direcotry,
the name of the bam without any path,
the real host where the file leaves (e.g. B<topmed3>),
and the index of the hostname (e.g. B<2> for topmed2).


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

