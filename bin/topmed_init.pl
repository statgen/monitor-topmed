#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name:	topmed_init.pl
#
# Description:
#   Use this program to check for files that are arriving
#   and initialize the NHLBI TOPMED database
#
#   Requires XML::Simple    (apt-get install libxml-simple-perl)
#
# ChangeLog:
#   $Log: nhlbi_init.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use DBIx::Connector;
use XML::Simple;
use POSIX qw(strftime);


my($me, $mepath, $mesuffix) = fileparse($0, '\.pl');
(my $version = '$Revision: 1.6 $ ') =~ tr/[0-9].//cd;
my $dbc;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my %opts = (
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    centers_table => 'centers',
    runs_table => 'runs',
    bamfiles_table => 'bamfiles',
    topdir => '/net/topmed/incoming/topmed',
    fakedcomment => 'Created from MD5 data',
    runcount => 0,
    bamcount => 0,
    bamcountruns => '',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose=i topdir=s center=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$me$mesuffix [options] updatedb\n" .
        "\nVersion $version\n" .
        "Monitor NHLBI data arriving in a directory (default=$opts{topdir}').\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

my $dbh = DBConnect($opts{realm});
if ($opts{verbose}) { print "$me$mesuffix Version $version, realm=$opts{realm}\n"; }

my $nowdate = time();

if ($fcn ne 'updatedb') { die "$me$mesuffix  - Invalid function '$fcn'\n"; }
chdir($opts{topdir}) ||
    die "$me$mesuffix Unable to CD to '$opts{topdir}': $!\n";

#--------------------------------------------------------------
#   Get all the known centers in the database
#--------------------------------------------------------------
my $centersref = GetCenters();
if ($opts{center}) {            # Debugging hook, only do one center
    foreach my $cid (keys %{$centersref}) { # Remove all but center we want
        if ($centersref->{$cid} ne $opts{center}) {
            delete($centersref->{$cid});
        }
    }
    if (! %$centersref) { die "$me$mesuffix Center '$opts{center}' unknown\n"; }
    print "$me$mesuffix - Running on center '$opts{center}' only\n";
}

#--------------------------------------------------------------
#   For each center watch for a new run to arrive
#--------------------------------------------------------------
foreach my $centerid (keys %{$centersref}) {
    my $c = $centersref->{$centerid};
    my $d = $opts{topdir} . '/' . $c;
    if (! chdir($d)) {
        warn "$me$mesuffix Unable to CD to '$d': $!\n";
        next;
    }
    #   Get all the known batch runs for this center
    my $sql = "SELECT runid,dirname,xmlfound FROM $opts{runs_table} WHERE centerid=$centerid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %knownruns = ();
    my $dir;
    for (my $i=1; $i<=$rowsofdata; $i++) {
        foreach my $href ($sth->fetchrow_hashref) {
            $dir = $href->{dirname};            
            $knownruns{$dir}{runid} = $href->{runid};
            $knownruns{$dir}{xmlfound} = $href->{xmlfound};
        }
    }
    #   Get list of all runs for this center
    #   Find new ones and add details to database
    my $dirsref = GetDirs('.');
    my $runid;
    foreach my $d (@{$dirsref}) {
        $runid = $knownruns{$d}{runid} || CreateRun($centerid, $d);
        if (! defined($runid)) {
            warn "How can runid not be defined?  dir=$d\n";
            next;
        }

        #   If XML exists and has not been found before, process it
        #   For now, always try to add new bams
        #   Someday this isn't supposed to be necessary because the XML
        #   will assure us that we know immediately what BAMs are expected
        if (! $knownruns{$d}{xmlfound}) { AddXML($runid, $d); }
        AddBams($runid, $d);
    }
}

$nowdate = strftime('%Y/%m/%d %H:%M', localtime);

if ($opts{runcount}) { print "$nowdate  Added $opts{runcount} runs\n"; }
if ($opts{bamcount}) { print "$nowdate  Added $opts{bamcount} bams from: $opts{bamcountruns}\n"; }
exit;

#==================================================================
# Subroutine:
#   CreateRun - Add details on this run to the database
#   It's not complete, but it can get us going
#
# Arguments:
#   d - directory (e.g. run name)
#   cid - center id
#
# Returns:
#   runid
#==================================================================
sub CreateRun {
    my ($cid, $d) = @_;

    my $sql = "INSERT INTO $opts{runs_table} " .
        "(centerid,dirname,comments,bamcount,dateinit) " .
        "VALUES($cid,'$d','',0,'$nowdate')";
    my $sth = DoSQL($sql);
    my $runid = $sth->{mysql_insertid};
    warn "Added run '$d'\n";
    $opts{runcount}++;
    #   Try to force permissions so things can work later. Can't trust the users
    chmod(0775, $d) || print "Unable to force permissions for '$d'. Too bad for you.\n"; 
    return $runid;
}

#==================================================================
# Subroutine:
#   AddBams - Add details for bams in this directory to the database
#
# Arguments:
#   runid - run id
#   d - directory
#
# Returns:
#   Boolean if any were added or not
#==================================================================
sub AddBams {
    my ($runid, $d) = @_;

    if (! -d $d) {
        warn "Unable to read directory '$d': $!";
        return 0;
    }
    
    #   Get all the known bams for this run
    my $sql = "SELECT bamname FROM $opts{bamfiles_table} WHERE runid=$runid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %knownbams = ();
    for (my $i=1; $i<=$rowsofdata; $i++) {
        foreach my $href ($sth->fetchrow_hashref) {
            $knownbams{$href->{bamname}} = 1;
        }
    }

    #   Merge all the MD5 files together and get list of bams and checksums
    #   There is noconsistency what people do here. Some files are bam checksum,
    #   some are the other way around. Other files have three field and again
    #   the order is anything you can imagine. Grrrr
    my %bams = ();
    my $bamcount = 0;
    if (! open(IN, "cat $d/*.md5 $d/*.MD5 *$d/Manifest.txt 2>/dev/null |")) {
         warn "No MD5 files found in directory '$d'. Maybe later, eh?\n";
         return 0;
    }
    my ($fn, $checksum, $f3);
    while (<IN>) {                  # Sometimes it's file checksum, sometimes not
        ($checksum, $fn, $f3) = split(' ',$_);
        if (! $fn) { warn "Surprising md5 syntax line: $_"; next; }
        if ($f3) { ($fn, $checksum) = ($checksum, $f3); }
        if ($fn !~ /bam$/) { ($checksum, $fn) = ($fn, $checksum); }
        if ($fn =~ /\//) { $fn = basename($fn); }
        if ($fn !~ /bam$/) {
             if ($opts{verbose}) { warn "Ignoring md5 syntax line: $_"; }
             next;
        }
        $bams{$fn} = $checksum;
        $bamcount++;
     }
    close(IN);
    if (! %bams) {
        if ($opts{verbose}) { warn "No bams found in '$d'\n"; }
        return 0;
    }

    #   Generate the bam database entry for each bam
    my $newbams = 0;
    foreach my $f (keys %bams) {
        if (exists($knownbams{$f})) { next; }   # Skip known bams
        $sql = "INSERT INTO $opts{bamfiles_table} " .
            "(runid,bamname,checksum,refname,expt_refname,expt_sampleid,dateinit) " .
            "VALUES($runid,'$f','$bams{$f}','UNKNOWN','UNKNOWN', 'UNKNOWN', $nowdate)";
        $sth = DoSQL($sql);
        $newbams++;
    }

    #   If we added bams, change the bamcount
    if ($newbams) {
        my $n = 
        $sql = "UPDATE $opts{runs_table}  SET bamcount=$bamcount WHERE runid=$runid";
        $sth = DoSQL($sql);
        $opts{bamcount} += $newbams;
        $opts{bamcountruns} .= $d . ' ';
    }
    return 1;
}


#==================================================================
# Subroutine:
#   AddXML - Update details on this run to the database from the XML
#   It's not complete since no one provides an XML file yet
#
# Arguments:
#   runid - run id
#   d - directory (e.g. run name)
#
# Returns:
#   boolean if the database was updated
#==================================================================
sub AddXML {
    my ($runid, $d) = @_;
    return 0;
    
#==================================================================
#   This has not been implemented yet 
#==================================================================
my $cid = 0;
    if (! opendir(DIR, $d)) {
        warn "Unable to read directory '$d': $!";
        return 0;
    }
    my @xmllist = grep { (/^\w.+\.xml/) } readdir(DIR);
    closedir DIR;

    #   There are two types of XML files, run.xml and experiment.xml
    #   Do the experiment data first, then run data
    #   Catch XML parse errors this first time through. We will
    #   re-parse these files again later, but only need catch errors
    #   the first time through.
    my %expt = ();
    my $emsg;
    foreach my $f (@xmllist) {
        $emsg = "XML parse error in '$d/$f': ";
        my $xml = new XML::Simple();
        my $myxml = eval { $xml->XMLin("$d/$f") };
        if ($@) {
            warn $emsg . "Parse XML file error:$@\n";
            exit(3);
        }

        if ($myxml->{RUN}[0]) { next; }
        if (! $myxml->{EXPERIMENT}[0]) {
            warn $emsg . "Unexpected XML file format\n";
            next;
        }
        if ($opts{verbose}) { warn "Parse experiment file '$d/$f'\n"; }
        #   Collect all the EXPERIMENT data for each BAM
        for (my $i=0; $i<=$#{$myxml->{EXPERIMENT}}; $i++) {
            my $alias   = $myxml->{EXPERIMENT}[$i]->{alias} ||
                warn $emsg . "alias not found\n";
            my $xref   = ${$myxml->{EXPERIMENT}[$i]->{DESIGN}->{SAMPLE_DESCRIPTOR}}{refname} ||
                warn $emsg . "refname not found\n";
            my $sampleid = ${$myxml->{EXPERIMENT}[$i]->{STUDY_REF}}{accession} ||
                warn $emsg . "accession not found\n";
            $expt{$alias}{xref} = $xref;
            $expt{$alias}{sampleid} = $sampleid;
        }
    }
    #   Now process the run files
    my %bams = ();
    foreach my $f (@xmllist) {
        my $xml = new XML::Simple();
        my $myxml = $xml->XMLin("$d/$f");
        if ($myxml->{EXPERIMENT}[0]) { next; }
        if (! $myxml->{RUN}[0]) {
            warn $emsg . "Unexpected XML file format\n";
            next;
        }
        #   Collect all the RUN data for each BAM
        $emsg = "XML parse error in '$d/$f': ";
        if ($opts{verbose}) { warn "Parse run file '$d/$f'\n"; }
        for (my $i=0; $i<=$#{$myxml->{RUN}}; $i++) {
            my $fn   = ${$myxml->{RUN}[$i]->{DATA_BLOCK}->{FILES}->{FILE}}{filename} ||
                warn $emsg . "filename not found\n";
            my $meth = ${$myxml->{RUN}[$i]->{DATA_BLOCK}->{FILES}->{FILE}}{checksum_method} ||
                warn $emsg . "checksum_method not found\n";
            my $cs   = ${$myxml->{RUN}[$i]->{DATA_BLOCK}->{FILES}->{FILE}}{checksum} ||
                warn $emsg . "checksum not found\n";
            my $rref = ${$myxml->{RUN}[$i]->{EXPERIMENT_REF}}{refname} ||
                warn $emsg . "refname not found\n";
            if ($meth ne 'MD5') {
                warn $emsg . "unsupport checksum method found '$meth' in '$d/$f'\n";
                return 0;
            }
            if (exists($bams{$fn})) {
                warn $emsg . "Duplicate BAM '$fn' defined in '$d/$f'\n";
                return 0;
            }
            if (! exists($expt{$rref})) {
                warn $emsg . "refname '$rref' not found in any EXPERIMENT FILE\n";
                return 0;
            }
            $bams{$fn}{checksum} = $cs;
            $bams{$fn}{refname} = $rref;
            $bams{$fn}{expt_refname} = $expt{$rref}{xref};
            $bams{$fn}{expt_sampleid} = $expt{$rref}{sampleid};
            $bams{$fn}{dateinit} = $nowdate;
        }
    }

    #   We now have all the information for the bams in this run
    #   Generate the run database entry. This entry might already
    #   exist, because we faked the entry with MD5 data.
    my @bamlist = keys %bams;
    if (! @bamlist) {
        warn "No BAMs found in XML files in '$d'\n";
        return;
    }
    my $sql = "SELECT runid FROM $opts{runs_table} WHERE centerid=$cid AND dirname='$d'";
    my $sth = DoSQL($sql,0);
    my $rowsofdata = $sth->rows();
    if ($rowsofdata) {                  # Update an existing row
        my $href = $sth->fetchrow_hashref;
        my $runid = $href->{runid};
        foreach my $bamfn (keys %bams) {
            my @cols = keys %{$bams{$bamfn}};
            $sql = "UPDATE $opts{runs_table} SET "; 
            foreach my $col (@cols) {
                $sql .= "$col='" . $bams{$bamfn}{$col} . "',";
            }
            $sql .= "dateinit='$nowdate' WHERE runid=$runid";
        }
    }
    else {
        #  Add a new row
        $sql = "INSERT INTO $opts{runs_table} (centerid,dirname,bamcount,comments,dateinit) " .
            "VALUES($cid,'$d'," . scalar(@bamlist) . ",'Created FROM XML data','$nowdate')";
        $sth = DoSQL($sql);
        my $runid = $sth->{mysql_insertid};
        #   Generate the bam database entry for each bam
        foreach my $bamfn (keys %bams) {
            my @cols = keys %{$bams{$bamfn}};
            my $sql = "INSERT INTO $opts{bamfiles_table} (runid,bamname," .
                join(',',@cols) . ") VALUES($runid,'$bamfn',";
            foreach my $col (@cols) {
                $sql .= "'" . $bams{$bamfn}{$col} . "',";
            }
            chop($sql);                         # Remove trailing ,
            $sql .= ')';
        }
    }
    $sth = DoSQL($sql);
    return;
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
    my %centers = ();

    #   Get all the known centers in the database
    my $sql = "SELECT centerid,centername FROM $opts{centers_table}";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    if (! $rowsofdata) { warn "$me$mesuffix No centers in '$opts{topdir}'\n"; }
    for (my $i=1; $i<=$rowsofdata; $i++) {
        my $href = $sth->fetchrow_hashref;
        $centers{$href->{centerid}} = $href->{centername};
    }
    return \%centers;
}

#==================================================================
# Subroutine:
#   GetDirs - Get list of non-dotted directories
#
# Arguments:
#   dirname
#
# Returns:
#   Reference to array of dir names
#==================================================================
sub GetDirs {
    my ($d) = @_;

    opendir(DIR, $d) ||
        die "Unable to read directory '$d': $!";
    my @dirlist = grep { (/^\w/ && -d "$d/$_") } readdir(DIR);
    closedir DIR;
    return \@dirlist;
}

#==================================================================
# Subroutine:
#   GetFiles - Get list of files (non-dotted)
#
# Arguments:
#   dirname
#
# Returns:
#   Reference to array of file names
#==================================================================
sub GetFiles {
    my ($d) = @_;

    opendir(DIR, $d) ||
        die "Unable to read directory '$d': $!";
    my @filelist = grep { (/^\w/ && -f "$d/$_") } readdir(DIR);
    closedir DIR;
    return \@filelist;
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
    if ($opts{verbose} > 1) { warn "DEBUG: SQL=$sql\n"; }
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    if ($DBI::err) {
        if (! $die) { return 0; }
        die "$me$mesuffix SQL failure: $DBI::errstr\n  SQL=$sql\n";
    }
    return $sth;
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

monitor.pl - Monitor

=head1 SYNOPSIS

  monitor.pl monitor            # Check directories, runs once
  monitor.pl -daemon monitor    # Check directories, runs repeatedly

=head1 DESCRIPTION

This program monitors directories for incoming data and then
updates a database with various status values.
The program B<moncmd.pl> supports simple commands to set key elements of the database.

The typical mode is to run this in 'daemon mode' by launching
it with a crontab entry, perhaps like this:

  01 06 * * * /usr/local/bin/monitor.pl -daemon monitor 2>&1 /dev/null

This will run for some time, checking the status of the directory specified
in by the B<-dir> option and updating the status in a database every
once in a while.

The queue of tasks is kept in a MySQL database.
See B<perldoc DBIx::Connector> for details defining the database.

=head1 OPTIONS

=over 4

=item B<-daemon>

Specifies the script should not exit after checking the directory.
It will wait for B<waitsecs> seconds and then relook at the directory again.
The program exits only after running for B<timeup>.

=item B<-dir>

Specifies the directory to be monitored.
This defaults to /data/CoreDump

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<webcd>.

=item B<-timeup nnn[MH]>

Specifies how long this program will remain running before exitting.
The number may have 'M' (minutes) or 'H' (hours) appended.
This defaults to B<6H>.
This is only used if B<daemon> is specified.

=item B<-verbose N>

Provided for developers to see additional information.

=item B<-waitsecs nnn>

Specifies how many seconds this program will wait before checking the database queue.
This defaults to B<300>.
This is only used if B<daemon> is specified.

=back

=head1 PARAMETERS

=over 4

=item B<monitor>

This directs the program to monitor the B<-dir> directory for changes
and update the database table specified by B<-realm>.

=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2010 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

