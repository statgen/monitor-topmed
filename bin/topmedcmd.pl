#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
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
use File::Basename;
use Cwd;
use Getopt::Long;
use DBIx::Connector;

my($me, $mepath, $mesuffix) = fileparse($0, '\.pl');
(my $version = '$Revision: 1.6 $ ') =~ tr/[0-9].//cd;
if (substr($mepath,0,1) ne '/') {
    my $d = getcwd();
    chdir($mepath) ||
        die "$me$mesuffix Unable to CD to '$d': $!\n";
    $mepath = getcwd();
    chdir($d);
}

my %VALIDVERBS = (                  # Valid verbs to database colum
    arrived => 'datearrived',
    md5verified => 'datemd5ver',
    baid    => 'datebai',
    qploted => 'dateqplot',
    backedup => 'datebackup',
    cramed => 'datecram',
    cp2ncbi => 'datecp2ncbi',
);
my %VALIDSTATUS = (                 # Valid status for the verbs
    requested => 1,
    completed => 1,
    started => 1,
    submitted => 1,
    failed => 1,
    cancelled => 1,
);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    centers_table => 'centers',
    runs_table => 'runs',
    topdir => '/net/topmed/incoming/topmed',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s runs=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$me$mesuffix [options]";
    warn "$m mark bamid arrived|md5verified|baid|qploted|backedup|cramed|cp2ncbi requested|submitted|completed|started|failed|cancelled\n" .
        "  or\n" .
        "$m unmark bamid [same list as mark]\n" .
        "  or\n" .
        "$m set bamid colname value\n" .
        "  or\n" .
        "$m show arrived\n" .
        "  or\n" .
        "$m export\n" .
        "\nVersion $version\n" .
        "Update the topmed database\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

my ($dbc);                      # For access to DB
my $dbh = DBConnect($opts{realm});
if ($opts{verbose}) { print "$me$mesuffix Version $version, realm=$opts{realm}\n"; }

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'mark')     { Mark(@ARGV); exit; }
if ($fcn eq 'unmark')   { UnMark(@ARGV); exit; }
if ($fcn eq 'set')      { Set(@ARGV); exit; }
if ($fcn eq 'show')     { Show(@ARGV); exit; }
if ($fcn eq 'export')   { Export(@ARGV); exit; }

die "$me$mesuffix  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   Mark($bamid, $op, $state)
#
#   Set dates in the status database
#==================================================================
sub Mark {
    my ($bamid, $op, $state) = @_;
    if (($bamid !~ /^\d+$/) || (! exists($VALIDVERBS{$op})) || (! exists($VALIDSTATUS{$state}))) {
        die "$me$mesuffix Invalid 'mark' syntax. Try '$me$mesuffix -help'\n";
    }

    #   Make sure this is a bam we know
    my $sth = DoSQL("SELECT bamid from $opts{bamfiles_table} WHERE bamid=$bamid", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$me$mesuffix - BAM '$bamid' is unknown\n"; }

    #   Set time for the verb based on this algorithm
    #   $t > 0     Task completed
    #   $t < 0     Task started
    #   $t = -1    Task failed
    #   $t = 0     Task requested
    #   $t = 1     Task cancelled
    #   $t = 2     Task submitted
    my $col = $VALIDVERBS{$op};
    my $val;
    my $done = 0;
    if ($state eq 'requested') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col='0' WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'completed') {
        $val = time();
        DoSQL("UPDATE $opts{bamfiles_table} SET $col='$val' WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'started') {
        $val = -time();
        DoSQL("UPDATE $opts{bamfiles_table} SET $col='$val' WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'failed') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col='-1' WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'cancelled') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col='1' WHERE bamid=$bamid");
        $done++;
    }
    if ($state eq 'submitted') {
        DoSQL("UPDATE $opts{bamfiles_table} SET $col='2' WHERE bamid=$bamid");
        $done++;
    }
    if ($done) {
        if ($opts{verbose}) { print "$me$mesuffix  'mark $bamid $op $state'  successful\n"; } 
    }
    else { die "$me$mesuffix Invalid state '$state' for '$op'. Try '$me$mesuffix -help'\n"; }
}

#==================================================================
# Subroutine:
#   UnMark($dirname, $op)
#
#   Reset dates in the statusdatabase
#==================================================================
sub UnMark {
    my ($bamid, $op) = @_;
    if (($bamid !~ /^\d+$/) || (! exists($VALIDVERBS{$op}))) {
        die "$me$mesuffix Invalid 'unmark' syntax. Try '$me$mesuffix -help'\n";
    }

    #   Make sure this is a bam we know
    my $sth = DoSQL("SELECT bamid from $opts{bamfiles_table} WHERE bamid=$bamid", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$me$mesuffix - BAM '$bamid' is unknown\n"; }

    my $col = $VALIDVERBS{$op};
    DoSQL("UPDATE $opts{bamfiles_table} SET $col=NULL WHERE bamid=$bamid");
    if ($opts{verbose}) { print "$me$mesuffix  'unmark $bamid $op'  successful\n"; }
}

#==================================================================
# Subroutine:
#   Export()
#
#   Generate a CSV file of possibly interesting data to STDOUT
#==================================================================
sub Export {
    #my ($center, $run) = @_;

    #   Generate header for CSV file
    my @cols = qw(bamid bamname studyname piname bamsize datemapping);
    my $s = 'CENTER,DIRNAME,FULLPATH';
    foreach (@cols) { $s .= uc($_) . ','; }
    chop($s);
    print $s . "\n";
    
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
                my $f = $opts{topdir} . "/$centername/$dirname/" .
                    $href->{bamname};
                #   See if this has arrived. Few states possible
                if (defined($href->{datearrived}) && ($href->{datearrived} ne '')) {
                    if ($href->{datearrived} =~ /\D/) { next; }  # Not numeric
                    if ($href->{datearrived} < 10) { next; }     # Already arrived
                }
                #   Convert mapping state into a string
                my $state = '';
                if (defined($href->{datemapping}) && $href->{datemapping} ne '') {
                    if ($href->{datemapping} eq '0')  { $state = 'requested'; }
                    if ($href->{datemapping} eq '1')  { $state = 'cancelled'; }
                    if ($href->{datemapping} eq '2')  { $state = 'submitted'; }
                    if ($href->{datemapping} < 0)     { $state = 'started'; }
                    if ($href->{datemapping} eq '-1') { $state = 'failed'; }
                    if ($href->{datemapping} > 10)    { $state = 'completed'; }
                }
                $href->{datemapping} = $state;               

                #   Show data for this BAM
                $s = "$centername,$dirname,$f,";
                foreach (@cols) {
                    if (! defined($href->{$_})) { $s .= ','; }
                    else { $s .= $href->{$_} . ','; }
                 }
                chop($s);
                print $s . "\n";
            }
        }
    }


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
        die "$me$mesuffix Invalid 'set' syntax. Try '$me$mesuffix -help'\n";
    }

    #   Make sure this is a bam we know
    my $sth = DoSQL("SELECT bamid from $opts{bamfiles_table} WHERE bamid=$bamid", 0);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { die "$me$mesuffix - BAM '$bamid' is unknown\n"; }

    DoSQL("UPDATE $opts{bamfiles_table} SET $col='$val' WHERE bamid=$bamid");
}

#==================================================================
# Subroutine:
#   Show($fcn)
#
#   Generate list of information from the database
#==================================================================
sub Show {
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
            my $sql = "SELECT bamid,bamname,datearrived FROM " .
                $opts{bamfiles_table} . " WHERE runid='$runid'";
            my $sth = DoSQL($sql);
            my $rowsofdata = $sth->rows();
            if (! $rowsofdata) { next; }
            for (my $i=1; $i<=$rowsofdata; $i++) {
                my $href = $sth->fetchrow_hashref;
                my $f = $opts{topdir} . "/$centername/$dirname/" .
                    $href->{bamname};
                #   See if this has arrived. Few states possible
                if (defined($href->{datearrived}) && ($href->{datearrived} ne '')) {
                    if ($href->{datearrived} =~ /\D/) { next; }  # Not numeric
                    if ($href->{datearrived} < 10) { next; }     # Not arrived
                }
                #   Run the command
                print "$href->{bamid} $centername $dirname $f\n";
            }
        }
    }
}

#==================================================================
# Subroutine:
#   DBConnect($realm)
#
#   Connect to our database using realm '$realm'. Return a DB handle.
#   Get the connection information from DBIx::Connector
#   Fully qualified realm file may be provided
#==================================================================
sub DBConnect {
    my ($realm) = @_;
    if (! $realm) { return 0; }
    #   Get the connection information FROM DBIx::Connector
    #   Fully qualified realm file may be provided
    if ($realm =~ /^(\/.+)\/([^\/]+)$/) {
        my $d = $1;
        $realm = $2;
        $dbc = new DBIx::Connector(-realm => $realm, -connection_dir => $d,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    else {
        $dbc = new DBIx::Connector(-realm => $realm,
            -dbi_options => {RaiseError => 1, PrintError => 1});
    }
    return $dbc->connect();
}

#==================================================================
# Subroutine:
#   DoSQL - Execute an SQL command
#
# Arguments:
#   sql - String of SQL to run
#   die - boolean if we should die on error
#
# Returns:
#   SQL handle for subsequent MySQL actions
#   Does not return if error detected
#==================================================================
sub DoSQL {
    my ($sql, $die) = @_;
    if (! defined($die)) { $die = 1; }
    if ($opts{verbose}) { warn "DEBUG: SQL=$sql\n"; }
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    if ($DBI::err) {
        if (! $die) { return 0; }
        die "$me$mesuffix SQL failure: $DBI::errstr\n  SQL=$sql\n";
    }
    return $sth;
}

#==================================================================
# Subroutine:
#   GetCenters - Get list of centers
#
# Arguments:
#   none
#
# Returns:
#   Reference to hash of center ids  to center names
#==================================================================
sub GetCenters {
#    my ($d) = @_;
    my %center2name = ();

    #   Get all the known centers in the database
    my $sql = "SELECT centerid,centername FROM $opts{centers_table}";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { warn "$me$mesuffix No centers found\n"; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{center} && $href->{centername} ne $opts{center}) { next; }
        $center2name{$href->{centerid}} = $href->{centername};
    }
    if ($opts{center}) {            # Only do one center
        if (! %center2name) { die "$me$mesuffix Center '$opts{center}' unknown\n"; }
        print "$me$mesuffix - Running on center '$opts{center}' only\n";
    }
    return \%center2name;
}

#==================================================================
# Subroutine:
#   GetRuns - Get list of runs for a center
#
# Arguments:
#   cid = center id
#
# Returns:
#   Reference to hash of run ids  to run dirnames
#==================================================================
sub GetRuns {
    my ($cid) = @_;
    my %run2dir = ();

    my $sql = "SELECT runid,dirname FROM $opts{runs_table} WHERE centerid=$cid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if ((! $rowsofdata) && $opts{verbose}) {
        warn "$me$mesuffix Found no runs for center '$cid'\n";
        return \%run2dir;
    }
    my %theseruns = ();
    if ($opts{runs}) {              # We only want some runs
        $opts{runs} =~ s/,/ /g;
        foreach my $r (split(' ', $opts{runs})) {
            $theseruns{$r} = 1;
            if (! $opts{onlyrun}) { 
                print "$me$mesuffix - Using only run '$r'\n";
                $opts{onlyrun}++;
            }
        }
    }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        if ($opts{runs} && (! exists($theseruns{$href->{dirname}}))) { next; }
        $run2dir{$href->{runid}} = $href->{dirname};
    }
    return \%run2dir;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedcmd.pl - Update the database for NHLBI TopMed

=head1 SYNOPSIS

  topmedcmd.pl mark 33 arrived completed   # BAM has arrived
  topmedcmd.pl unmark 33 arrived           # Reset BAM has arrived

  topmedcmd.pl set 33 jobidqplot 123445    # set jobidqplot in bamfiles

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
This defaults to B<.db_connections/monitor> in the same directory as
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

=over 4

=item B<mark bamid dirname  [verb] [state]>
Use this to set the state for a particular BAM file.
Mark will set a date for the process (e.g. arrived sets datearrived)
and unmark will set that entry to NULL.

Here is the list of supported verbs:
 arrived
 md5verified
 backedup
 cramed
 baid
 qploted
 cp2ncbi

Here is the list of supported states:
 requested
 submitted
 completed
 started
 failed
 cancelled

Note that this database field can have some surprising values. 
 If time is -1, the task failed.
 If time is zero, the task has been requested.
 If time is 1, the task was cancelled.
 If time is 2, the task was submitted to the batch system.
 If time is positive, the task was completed.
 q
 If time is negative, the task was started.

=item B<set bamid columnname value>
Use this to set the value for a column for a particular BAM file.

=item B<unmark bamid [verb]>
Use this to reset the state for a particular BAM file to NULL, the default
database value.

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

