#!/usr/bin/perl
###################################################################
#
# Name:	topmed_init.pl
#
# Description:
#   Use this program to check for files that are arriving
#   and initialize the NHLBI TOPMED database
#   This program can work with topmed and inpsyght genome data
#	as well as RNS sequence and methylation data
#
# ChangeLog:
#   $Log: topmed_init.pl,v $
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
);
use My_DB;
use Topmed::Get;
use Getopt::Long;

use POSIX qw(strftime);
use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
#   Pre-check options for project
for (my $i=0; $i<=$#ARGV; $i++) {
    if (($ARGV[$i] eq '-p' || $ARGV[$i] eq '-project') && defined($ARGV[$i+1])) {
        $ENV{PROJECT} = $ARGV[$i+1];
        last;
    }
}
if (! -d "/usr/cluster/$ENV{PROJECT}") { die "$Script - Environment variable PROJECT '$ENV{PROJECT}' incorrect\n"; }

my $NOTSET = 0;

our %opts = (
	datatype => 'genome',		# What kind of data we are handling
    realm => "/usr/cluster/$ENV{PROJECT}/etc/.db_connections/$ENV{PROJECT}",
    centers_table => 'centers',
    runs_table => 'runs',
    runs_pkey => 'runid',
    samples_table => 'bamfiles',
    samples_pkey => 'bamid',
    topdir => "/net/$ENV{PROJECT}/incoming/$ENV{PROJECT}",
    runcount => 0,
    count => 0,
    countruns => '',
    arrivedsecs => 86400*21,	# If no new files in a bit, stop looking at run
    ignorearrived => 0,     	# If set, ignored arrived database field
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose center=s run=s ignorearrived project=s datatype=s
    )) || die "Failed to parse options\n";

#	Set %opts based on datatype
if ($opts{datatype} eq 'rnaseq') {
	$opts{runs_table} = 'tx_projects',
	$opts{runs_pkey} = 'rnaprojectid',
	$opts{samples_table} = 'tx_samples',
	$opts{samples_pkey} = 'txseqid',
	$opts{files_table} = 'tx_files',
	$opts{files_pkey} = 'fileid',
	$opts{topdir} = "/net/$ENV{PROJECT}/incoming/$ENV{PROJECT}/rnaseq";
}
if ($opts{datatype} =~ /^meth/) {
	$opts{datatype} = 'methyl',
	$opts{runs_table} = 'methyl_projects',
	$opts{runs_pkey} = 'methylprojectid',
	$opts{samples_table} = 'methyl_samples',
	$opts{samples_pkey} = 'methylid',
	$opts{files_table} = 'methyl_batch',
	$opts{files_pkey} = 'methylbatchid',
	$opts{topdir} = "/net/$ENV{PROJECT}/incoming/$ENV{PROJECT}/methylation";
}
if ($opts{verbose}) { warn "$Script - processing datatype=$opts{datatype}\n"; }

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] [-datatype genome|rnaseq|methyl] updatedb\n" .
        "Monitor TOPMed data arriving in a directory (default=$opts{topdir}').\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

my $dbh = DBConnect($opts{realm});
my $nowdate = time();
if ($fcn ne 'updatedb') { die "$Script  - Invalid function '$fcn'\n"; }
chdir($opts{topdir}) ||
    die "$Script Unable to CD to '$opts{topdir}': $!\n";

#--------------------------------------------------------------
#   For each center watch for a new data to arrive
#--------------------------------------------------------------
my $centersref = GetCenters();
foreach my $centerid (keys %{$centersref}) {
    my $c = $centersref->{$centerid};
    my $d = $opts{topdir} . '/' . $c;
    if (! chdir($d)) {
        #warn "$Script Unable to CD to '$d': $!\n";
        next;
    }
    #   Get all the known batch data for this center that are not fully arrived
    my $sql = "SELECT $opts{runs_pkey},dirname,arrived FROM $opts{runs_table} WHERE centerid=$centerid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %knownruns = ();
    my $dir;
    for (my $i=1; $i<=$rowsofdata; $i++) {
        foreach my $href ($sth->fetchrow_hashref) {
            $dir = $href->{dirname};            
            $knownruns{$dir}{$opts{runs_pkey}} = $href->{$opts{runs_pkey}};
            $knownruns{$dir}{arrived} = $href->{arrived};
            if ($opts{ignorearrived}) { $knownruns{$dir}{arrived} = 'N'; }
        }
    }
    #   Get list of all runs for this center
    #   Find new ones and add details to database
    my $dirsref = GetDirs('.');
    my $runid;
    foreach my $d (@{$dirsref}) {
    	if ($opts{verbose}) { warn "$Script - checking '$d'\n"; }
        $runid = $knownruns{$d}{$opts{runs_pkey}} || CreateRun($centerid, $d);
        if (! defined($runid)) {
            if ($opts{verbose}) { warn "$Script ID not defined for dir=$d\n"; }
            next;
        }
        if ($opts{run} && $opts{run} ne $d) { next; }
        #   Check if this run has arrived, no need to look at it further
        if (defined($knownruns{$d}{arrived}) && $knownruns{$d}{arrived} eq 'Y') { next; }

		#warn "$Script - d=$d $opts{datatype} would have been called. arrived='$knownruns{$d}{arrived}'<br>\n"; next;
        #	Add database entries for each sample in this set
        if ($opts{datatype} eq 'genome') { AddBams($runid, $d); next;}
        if ($opts{datatype} eq 'rnaseq') { AddRNAProjects($runid, $d); next;}
        if ($opts{datatype} eq 'methyl') { AddMethylBatches($runid, $d); next;}
        warn "$Script Unable to process data '$opts{datatype}' for dir=$d\n";
    }
}

$nowdate = strftime('%Y/%m/%d %H:%M', localtime);
if ($opts{runcount}) { print "$nowdate  Added $opts{runcount} runs\n"; }
if ($opts{count})    { print "$nowdate  Added $opts{count} samlpes from: $opts{countruns}\n"; }
exit;

#==================================================================
# Subroutine:
#   AddMethylBatches - Add database entries for each set of files/samples in this directory
#
# Arguments:
#   methylprojectid - parent id
#   d - directory
#
# Returns:
#   Boolean if any samples were added or not
#==================================================================
sub AddMethylBatches {
    my ($methylprojectid, $d) = @_;
    if (! -d $d) {
        print "$Script - Unable to read directory '$d': $!\n";
        return 0;
    }
   	if ($opts{verbose}) { warn "$Script - AddProjects($methylprojectid, $d)\n"; }

	#------------------------------------------------------------------
    #	We assume if there is a Manifest, all the files are in place.
    #	Read the MD5 and collect into useful data structures
	#------------------------------------------------------------------
    my $manifestfile = "$d/Manifest.txt";
    my %batchnames = ();
	if (! open(IN, $manifestfile)) {
		warn "$Script - Unable to read file '$manifestfile': $!\n";
		return 0;
	}	
	my $totalbatch = 0;
	while (my $l = <IN>) {          # Read md5 checksum file
		chomp($l);
		my ($name, $checksum) = NormalizeMD5Line($l, $manifestfile, 1);
		if (! $checksum) {
			warn "$Script - Unable to parse checksum from '$manifestfile'\n";
			return 0;
		}
		my ($bname, $n) = ('', 99);
		if ($name =~ /^(\S+)-LEVEL(\d)/) { ($bname,$n) = ($1, $2); }
		if ($n == 99) {
			warn "$Script - zip name could not be parsed '$manifestfile'\n";
			return 0;
		}
		my $size = (stat "$d/$name")[7];
		if (! $size) {
			warn "$Script - zip name '$d/$name' does not exist\n";
			return 0;
		}

		#	Save the important stuff
		my $k = 'file' . $n; 
		$batchnames{$bname}{$k}{filename} = $name;
		$batchnames{$bname}{$k}{checksum} = $checksum;
		$batchnames{$bname}{$k}{size} = $size;
		if ($n == 0) { $totalbatch++; }
	}
	close(IN);

	#------------------------------------------------------------------
	#	Manifest parsed. Create methyl_batch entries
	#------------------------------------------------------------------
	foreach my $bname (keys %batchnames) {
		foreach my $n (qw(0 1 2 3)) {
			my $k = 'file' . $n; 
			if (! exists($batchnames{$bname}{$k})) {
				warn "$Script - Batch '$bname is missing LEVEL $n entry\n";
				return 0;
			}
		}
	}
	my $methylbatchid;
	foreach my $bname (keys %batchnames) {
		#   This needs to be restartable, so I must first check if this file already exists
		my $sql = "SELECT $opts{files_pkey} from $opts{files_table} WHERE " .
		    "methylprojectid='$methylprojectid' AND batchname='$bname'";
		my $sth = DoSQL($sql);
		my $href = $sth->fetchrow_hashref;
		my $x = $href->{$opts{files_pkey}};
		if ($x) {
		    if ($opts{verbose}) { warn "$Script - Ignoring existing file '$bname' [$x]\n"; }
		    next;
        }
		$sql = "INSERT INTO $opts{files_table} " .
			"(methylprojectid,batchname,count,dateinit, " .
			"file0, file0checksum, file0size," .
			"file1, file1checksum, file1size," .
			"file2, file2checksum, file2size," .
			"file3, file3checksum, file3size)";
		my $values = " VALUES($methylprojectid,'$bname',0,$nowdate,";
		for (my $i=0; $i<=3; $i++) {
			my $k = 'file' . $i; 
			$values .= "'" . $batchnames{$bname}{$k}{filename} . 
				"','" . $batchnames{$bname}{$k}{checksum} .
				"'," . $batchnames{$bname}{$k}{size} . ',';
		}
		chop($values);
		$values .= ")";
		$sth = DoSQLv($sql . $values);
		if ($sth) { $methylbatchid  = SQL_Last_Insert($sth); }	# Index to this sample
		if (! $methylbatchid) {
			warn "$Script - INSERT for '$bname failed\n";
			return 0;
		}
		#	Save index to entry so we can create samples
		$batchnames{$bname}{id} = $methylbatchid
	}

	#------------------------------------------------------------------
	#	Batch entries created, Create samples by tearing apart LEVEL0 zip
	#------------------------------------------------------------------
	my $tmpdir = "/run/shm/$$";
	mkdir $tmpdir;
	my $totalsamples = 0;
	foreach my $bname (keys %batchnames) {
		my $cmd = "unzip -d $tmpdir $d/$batchnames{$bname}{file0}{filename} '*/Manifest*.txt'";
		my $manfile;
		foreach my $l (split("\n", `$cmd`)) {
			if ($l =~ /inflating:\s+(\S+)/) { $manfile = $1; last; }
		}
		#	Be sure to remove DOS \r characters. Some files have them, some not
		#	Some files have NO \n, only \r !  Piece of crap data !
		if (! open(IN, $manfile)) {
			warn "$Script - Unable to read LEVEL0 Manifest file '$manfile'\n";
			return 0;
		}
		open(OUT, '>' . $manfile . '.out');
		while (my $l = <IN>) { $l =~ s/\r/\n/g; print OUT $l; }
		close(IN);
		close(OUT);
		if (! open(IN, $manfile . '.out')) {
			warn "$Script - Unable to read LEVEL0 Manifest file '$manfile'\n";
			return 0;
		}
		# 	Parse the header line to figure out what fields are in what columns
		my $hdr = <IN>;
		chomp($hdr);
		my %cols = (
			'grn.idat'        => 99,		# 99 means the column is not known
			'red.idat'        => 99,
			'.sdf'            => 99,
			'meth.raw.csv'    => 99,
			'unmeth.raw.csv'  => 99,
			'beta.raw.csv'    => 99,
			'pvalue.csv'      => 99,
			'rgset.snp.csv'   => 99,
			'meth.noob.csv'   => 99,
			'unmeth.noob.csv' => 99,
			'beta.noob.csv'   => 99,
			'noob.correct.csv'=> 99,
			'sampleid'        => 1,			# I hope these are always here
			'array'           => 2,
			'p05'             => 3,
			'p01'             => 4,
		);
		my %cols2dbcol = (					# Map hdr column to database column
			'grn.idat'        => 'grn_idat',
			'red.idat'        => 'red_idat',
			'.sdf'            => 'sdf',
			'meth.raw.csv'    => 'meth_raw',
			'unmeth.raw.csv'  => 'unmeth_raw',
			'beta.raw.csv'    => 'beta_raw',
			'pvalue.csv'      => 'pvalue',
			'rgset.snp.csv'   => 'rgset_snp',
			'meth.noob.csv'   => 'meth_noob',
			'unmeth.noob.csv' => 'unmeth_noob',
			'beta.noob.csv'   => 'beta_noob',
			'noob.correct.csv'=> 'noob_correct',
			'sampleid'        => 'expt_sampleid',
			'array'           => 'array',
			'p05'             => 'p05',
			'p01'             => 'p01',
		);
		my @c = split("\t", $hdr);
		for (my $i=1; $i<=$#c; $i++) {
			foreach my $col (keys %cols) {
				if ($c[$i] =~ /\($col/i) { $cols{$col} = $i; last; }
			}
		}
		#	See if we found a column index for every %cols
		foreach my $c (keys %cols) {
			if ($cols{$c} eq '99') {
				warn "$Script - Unable to find index for column '$c' in LEVEL0 Manifest file '$manfile'\nHDR=$hdr\n";
				return 0;
			}
		}
		#------------------------------------------------------------------
		#	Manifest hdr parsed. Create samples entries
		#------------------------------------------------------------------
		my $methylid;
		my $samplecount = 0;
		while (my $l = <IN>) {          # Read line from Manifest file
			if ($l =~ /^</) { next; }	# Ignore start/end like <B0001> and </B0001>
			if ($l !~ /\S/) { next; }	# and blank lines
			chomp($l);
			my @c = split("\t", $l);
			if ($c[$cols{'grn.idat'}] !~ /^\S+\-(\S+)\-(\S{3})_(\S{6})/) {
				warn "$Script - Unable to parse grn.idat column '" . $c[$cols{'grn.idat'}] .
					"' from:\nLINE=$l\n";
				if ($c[$cols{'grn.idat'}] !~ /Failure/) { return 0; }
				next;
			}
		    #   This needs to be restartable, so I must first check if this sample already exists
		    my $sql = "SELECT $opts{files_pkey} from $opts{files_table} WHERE " .
		        "methylprojectid=$batchnames{$bname}{id} AND expt_sampleid='" . $c[$cols{'sampleid'}] . "'";
		    my $sth = DoSQL($sql);
		    if ($sth) {
		        if ($opts{verbose}) { warn "$Script - Ignoring existing sample '" . $c[$cols{'sampleid'}] . "'\n"; }
		        next;
            }
			my ($type, $version, $rowcol) = ($1, $2, $3);
			$sql = "INSERT INTO $opts{samples_table} " .
				"(methylbatchid,expt_sampleid,array,p05,p01,type,version,rowcol,";
			my $values .= "VALUES(" . $batchnames{$bname}{id} . "," .
				"'$c[$cols{'sampleid'}]'," .
				"'$c[$cols{'array'}]'," .
				"$c[$cols{'p05'}]," .
				"$c[$cols{'p01'}]," .
				"'$type', '$version', '$rowcol',";
			foreach my $k (keys %cols){
				if ($cols{$k} <= 4) { next; }
				#	Mismatch between %cols index and SQL column name :-(
				$sql .= $cols2dbcol{$k} . ',';	# Get DB column
				if ($c[$cols{$k}]) { $values .= "'Y',"; }
				else { $values .= "'N',"; }
			}
			chop($sql);					# Drop trailing comma
			$sql .= ') ';
			chop($values);
			$values .= ')';
			$sth = DoSQL($sql . $values);
			$methylid  = SQL_Last_Insert($sth);
			if (! $methylid) {	#	Sample exists I hope. Need to know index for it
				warn "$Script - No database record created for '$c[$cols{'sampleid'}]'\n";
				return 0;
			}
			$totalsamples++;
		}
		close(IN);

		#	Update count of samples for this batch
		my $sql = "UPDATE $opts{files_table} SET count=$samplecount WHERE $opts{files_pkey}=$batchnames{$bname}{id}";
		if ($opts{verbose}) { print "SQL=$sql\n"; }
		my $sth = DoSQLv($sql);	
		$samplecount++;               # fix samplecount and/or totalsamples
	}

	#	Update count of samples for this project
	my $sql = "UPDATE $opts{runs_table} SET arrived='Y' WHERE methylprojectid=$methylprojectid";
	if ($opts{verbose}) { print "SQL=$sql\n"; }
	my $sth = DoSQLv($sql);
	$sql = "UPDATE $opts{runs_table} SET count=$totalsamples WHERE methylprojectid=$methylprojectid";
	if ($opts{verbose}) { print "SQL=$sql\n"; }
	$sth = DoSQL($sql);

	#	Give summary of what happened
	system("rm -rf $tmpdir");
	print "Created $totalsamples samples in $totalbatch batch files\n";
    return $totalsamples;
}

#==================================================================
# Subroutine:
#   AddRNAProjects - Add database entries for each set of samples in this directory
#
# Arguments:
#   rnaprojectid - parent id
#   d - directory
#
# Returns:
#   Boolean if any samples were added or not
#==================================================================
sub AddRNAProjects {
    my ($rnaprojectid, $d) = @_;
    if (! -d $d) {
        print "$Script - Unable to read directory '$d': $!\n";
        return 0;
    }
   	if ($opts{verbose}) { warn "$Script - AddProjects($rnaprojectid, $d)\n"; }
	if (! $rnaprojectid) { die "$Script - AddRNAProjects called without a valid rnaprojectid\n"; }
	my $nwgcref = GetNWGC($d, 'lookup.*.tab', 'Manifest.txt');
	if (! $nwgcref) {
		warn "$Script - Unable to parse NWGC from data in '$d'\n";
		return 0;
	}

	#------------------------------------------------------------------
    #	We assume if there is a Manifest, all the files are in place.
    #	Read the MD5 and collect into useful data structures
    #	Assumes each file in this data is uniquely named
	#------------------------------------------------------------------
    my $manifestfile = "$d/Manifest.txt";
    my %expt_sampleids = ();
    my %fullfiletodata = ();
    my $inputstyle = 'olduw';				# Which format is this Manifest
	if (! open(IN, $manifestfile)) {
		warn "$Script - Unable to read file '$manifestfile': $!\n";
		return 0;
	}
	while (my $l = <IN>) {          # Read md5 checksum file
		chomp($l);
		my ($fullfn, $check) = NormalizeMD5Line($l, $manifestfile, 1);
		if (! $check) {
			warn "$Script - Unable to parse checksum from '$manifestfile'\n";
			return 0;
		}
		if ($fullfn =~ /^TOR/) { $inputstyle = 'broad'; }
		if ($inputstyle eq 'broad' && $fullfn !~ /\.bam|\.bai$/) { next; }	# Ignore files we are not interested in

		#	Until we know otherwise, assume first part of $fn is the sampleid
		#	Broad uses TOR + sampleid
		my $nwgc = '';
		if ($fullfn =~ /^(\w+)\./) { $nwgc = $1; }
		if (! $nwgc) {
			warn "$Script - Unable to parse filename from '$manifestfile'  Line=$l\n";
			return 0;
		}
		my $sampleid = $nwgcref->{$nwgc}->{expt_sampleid};
		if (! $sampleid) {
			warn "$Script - Unable to find expt_sampleid for '$nwgc' in '$manifestfile'\nLINE=$l\n";
			return 0;
		}

		#	Save what we will need, list of all sampleids and all files for samples
		$fullfiletodata{$fullfn} = $check . ' ' . $nwgc . ' ' . $sampleid;
		$expt_sampleids{$sampleid} = $nwgcref->{$nwgc};	# Cheap way to get unique list
	}
	close(IN);

	#------------------------------------------------------------------
	#	Create all sample entries from manifest, capture txseqid
	#------------------------------------------------------------------
	my $totalsamples = 0;
	my $txseqid;
	foreach my $samp (sort keys %expt_sampleids) {
		$totalsamples++;
		#   This needs to be restartable, so I must first check if this sample already exists
		my $sql = "SELECT $opts{samples_pkey} from $opts{samples_table} WHERE " .
		    "rnaprojectid=$rnaprojectid AND expt_sampleid='$samp'";
		my $sth = DoSQL($sql);
		my $href = $sth->fetchrow_hashref;
		$txseqid = $href->{txseqid};
		if ($txseqid) {
		    if ($opts{verbose}) { warn "$Script - Ignoring existing sample '$samp' [$txseqid]\n"; }
		    $expt_sampleids{$samp}{txseqid} = $txseqid;
		    next;
        }
		#   Create database record for this sample. It might already exist
		$sql = "INSERT INTO $opts{samples_table} " .
			"(rnaprojectid,expt_sampleid,rnasubject,fileprefix,notes,dateinit) " .
			"VALUES($rnaprojectid,'$samp','$expt_sampleids{$samp}->{rnasubject}'," .
			"'$expt_sampleids{$samp}->{fileprefix}'," .
			"'$expt_sampleids{$samp}->{notes}',$nowdate)";
		$sth = DoSQLv($sql);
		if ($sth) { $txseqid = SQL_Last_Insert($sth); }		# Index to this sample
		if (! $txseqid) {	#	Sample exists I hope. Need to know index for it
			$sql = "SELECT txseqid FROM $opts{samples_table} " .
				"WHERE expt_sampleid='$samp'";
			$sth = DoSQL($sql);
			$href = $sth->fetchrow_hashref;
			$txseqid = $href->{txseqid};
		}
		if (! $txseqid) {
			warn "$Script - No database record created for '$samp'\n";
			return 0;
		}
		$expt_sampleids{$samp}{txseqid} = $txseqid;
	}
	#	Update count of samples for this project
	my $sql = "UPDATE $opts{runs_table} SET count=$totalsamples WHERE rnaprojectid=$rnaprojectid";
	if ($opts{verbose}) { print "SQL=$sql\n"; }
	my $sth = DoSQLv($sql);

	#------------------------------------------------------------------
	#	Create all file entries
	#------------------------------------------------------------------
	my $totalfilecount = 0;
	my $samplefilecount = 0;
	my $lasttxseqid = '';
	foreach my $fullfn (sort keys %fullfiletodata) {
		my ($check, $fileprefix, $sampleid) = split(' ', $fullfiletodata{$fullfn});
		my $txseqid = $expt_sampleids{$sampleid}{txseqid};
		#	Is this a new sample? If so, update file count for previous sample
		if ($txseqid ne $lasttxseqid) {
			if ($lasttxseqid) {				# Correct count for last sample
				my $sql = "UPDATE $opts{samples_table} SET count=$samplefilecount WHERE txseqid=$lasttxseqid";
				if ($opts{verbose}) { print "SQL=$sql\n"; }
				my $sth = DoSQLv($sql);
				$samplefilecount = 0;
			}
		}
		$lasttxseqid = $txseqid;

		#	Now create database record in $opts{files_table} for this file
		$samplefilecount++;
		$totalfilecount++;
		if (AddRNAFile($txseqid, $fullfn, 'N', $check)) { $samplefilecount--; $totalfilecount--; };
		if ($fullfn !~ /\.tar/) { next; }

		# 	This is a tar file, create records for each file in tar
		my $path = "$d/releaseFiles/$fullfn";
		if (! -f $path) {
			warn "$Script - TAR file '$path' was not found\n";
			return 0;
		}
		warn "  Processing tar $path for sample index $txseqid\n";
		my $z = '';
		if ($fullfn =~ /.gz$/) { $z = '-z'; }
		foreach my $l (split("\n", `tar $z -t -f $path`)) {	# Get filenames in tar
			#if ($l =~ /.+\/($fileprefix\/.+)/) { $l = $1; } 	# Strip useless prefix
			if ($l =~ /\/$/) { next; }		# Not interested in directory names
			if (AddRNAFile($txseqid, $l, 'Y', ' ')) { $samplefilecount++; $totalfilecount++; };
		}
	}
	#	Update count of files for last sample
	$sql = "UPDATE $opts{runs_table} SET arrived='Y' WHERE rnaprojectid=$rnaprojectid";
	if ($opts{verbose}) { print "SQL=$sql\n"; }
	$sth = DoSQLv($sql);
	$sql = "UPDATE $opts{samples_table} SET count=$samplefilecount WHERE txseqid=$txseqid";
	if ($opts{verbose}) { print "SQL=$sql\n"; }
	$sth = DoSQLv($sql);

	#	Give summary of what happened
	print "Created $totalsamples samples with $totalfilecount files in project '$d'\n";
    return $totalsamples;
}

#==================================================================
# Subroutine:
#   GetNWGC - Parse tab file which maps NWGC ids to TOR ids
#	Column names must map to our database columns
#
# Arguments:
#   dir - directory to RNA data
#   oldtabfile - tab file to parse  (probably early architecture)
#	manifest - Manifest for data, infer NWGC from filename
#
# Returns:
#   Reference to a hash of fileprefix to a hash of other column data
#==================================================================
sub GetNWGC {
    my ($dir, $oldtabfile, $manifest) = @_;

	my %arrayindex2dbcol = ();
	#	Should only be one tab file if this is early architecture
	my @file = split("\n", `/bin/ls $dir/$oldtabfile 2>&1`);
	if ($#file == 0) {
		if ($file[0] !~ /No such/) { return ParseTab($file[0]); }
	}

	#	Broad architecture is later and probably the one we'll use
	return ParseBroad("$dir/$manifest");
}

#==================================================================
# Subroutine:
#   ParseBroad - Parse manifest extracting NWGC ids to TOR ids
#
# Arguments:
#   file - file to parse (broad architecture)
#
# Returns:
#   Reference to a hash of nwgc to a hash of other column data
#==================================================================
sub ParseBroad {
    my ($file) = @_;
	my %parsedfile = ();

	#	This is Broad style input where the TOR name has the NWGC id in it
	if (! open(IN, $file)) {
		warn "$Script - Unable to read MANIFEST file '$file': $!\n";
		return \%parsedfile;;
	}
	# 	Read tab delimited manifest e.g.   md5 TOR123456.stuff.bam
	while (my $l = <IN>) {
		chomp($l);
		my ($md5, $f) = split(' ', $l);
		my $nwgc;
		if ($f =~ /\.bai|gz|tsv\s*$/) { next; }	# Ignore files we are not interested in
		if ($f =~ /^(TOR\d+)\..+\.bam\s*$/) { $nwgc = $1; }
		if (! $nwgc) {
			warn "$Script - Ignoring parse of NWGC from '$file' line: $l\n";
			next;
		}
		my %tabhash = (					# Hard code these, since not CSV input
			notes => '',
			fileprefix => 'bam_files',
			rnasubject => '',
			expt_sampleid => $nwgc,
		);
		$parsedfile{$nwgc} = \%tabhash;	# Save hash of this data
	}
	close(IN);
	return \%parsedfile;				# Return reference of hash of hashes
}

#==================================================================
# Subroutine:
#   ParseTab - Parse tab file which maps NWGC ids to TOR ids
#	Column names must map to our database columns
#
# Arguments:
#   tabfile - tab file to parse  (early architecture)
#
# Returns:
#   Reference to a hash of nwgc to a hash of other column data
#==================================================================
sub ParseTab {
    my ($tabfile) = @_;
	my %parsedtabfile = ();

	#------------------------------------------------------------------
    #	Figure out which columns map to which database columns
	#------------------------------------------------------------------
	my %arrayindex2dbcol = ();
	if (! open(IN, $tabfile)) {
		warn "$Script - Unable to read TAB file '$tabfile': $!\n";
		return \%parsedtabfile;;
	}
	# First line is header which will vary by the whim of each project
	my %tabcol2dbcol = (
		notes => 'notes',
		'nwgc sample id' => 'fileprefix',
		'investigator id' => 'rnasubject',
		torid => 'expt_sampleid',
	);
	my $l = lc(<IN>);
	chomp($l);
	my @hdrcols = split("\t", $l);
	for (my $i=0; $i<=$#hdrcols; $i++) {
		if (! exists($tabcol2dbcol{$hdrcols[$i]})) {
			warn "$Script - Ignoring column '$hdrcols[$i]' in TAB file '$tabfile' HDR='$l'\n";
			next;
		}
		$arrayindex2dbcol{$i} = $tabcol2dbcol{$hdrcols[$i]};	# Index to db column name
	}
	if (scalar(keys %tabcol2dbcol) != scalar(keys %arrayindex2dbcol)) {
		warn "$Script - Unable to which TAB file columns to map to which database " .
				"columns. TAB file='$tabfile'. HDR=$l\n";
	}

	# 	Now we know which tab column maps to which database column
	#	Save as hash of hashes.  fileprefix to hash of other fields
	while ($l = <IN>) {
		chomp($l);
		my @c = split("\t", $l);
		my %tabhash = ();
		my $fp = '';
		foreach my $cindex (keys %arrayindex2dbcol) {
			my $dbcol = $arrayindex2dbcol{$cindex};	# DB column name
			if ($dbcol eq 'fileprefix') { $fp = $c[$cindex]; }
			$tabhash{$dbcol} = $c[$cindex];
		}
		$parsedtabfile{$fp} = \%tabhash;		# Save hash of this line
	}
	close(IN);
	return \%parsedtabfile;				# Return reference of hash of hashes
}

#==================================================================
# Subroutine:
#   AddRNAFile - Add new or update a filename in $opts{files_table}
#
# Arguments:
#   txseqid - index to a sampleid
#	filename - name of file for sample
#	intar - flag if this file is in a tar
#	checksum - checksum for this file
#
# Returns:
#   Boolean INSERT was successful
#==================================================================
sub AddRNAFile {
    my ($txseqid, $filename, $intar, $checksum) = @_;

	#	This filename may not even be in the database yet
    my $sql = "SELECT fileid,txseqid FROM $opts{files_table} " .
    	"WHERE filename='$filename'";
    my $sth = DoSQL($sql,0);
    my $href = $sth->fetchrow_hashref;
    if ($href) {					# filename already exists
    	warn "  File '$filename' already in database (fileid=$href->{fileid} " .
    		"txsequid=$href->{txseqid})\n";
    	return 0;
    }

   	$sql = "INSERT INTO $opts{files_table} " .
    	"(txseqid,filename,intar,checksum,dateinit) " .
       	"VALUES($txseqid,'$filename','$intar', '$checksum', $nowdate)";
   	if ($opts{verbose}) { print "SQL=$sql\n"; }
   	$sth = DoSQLv($sql);
    return 1;
}

#==================================================================
# Subroutine:
#   CreateRun - Add details on this run to the database
#   Make sure this directory is runnable. If not, no new runid
#
# Arguments:
#   d - directory (e.g. run name)
#   cid - center id
#
# Returns:
#   runid or undef
#==================================================================
sub CreateRun {
    my ($cid, $d) = @_;
    #   Runs with a magic name are ignored
    if ($d eq 'upload') { return undef(); }

    #   Do nothing until this is owned by the right user
    my @s = stat($d);
    if ($s[5] != 2307982 && $s[5] != 1106) {
        if ($opts{verbose}) { print "$Script - Ignoring directory '$d' until it is owned by topmed\n"; }
        return undef();
    }

    #   Try to write in this directory. Can't trust the users
    my $touchfile = '.test';
    if ( ! open(TOUCH, '>' . $touchfile)) {
        if ($opts{verbose}) { print "$Script - Ignoring non-writable directory '$d'\n"; }
        return undef();
    }
    if ( ! print TOUCH 'test if writable') {
        if ($opts{verbose}) { print "$Script - Ignoring non-writable directory '$d'\n"; }
        return undef();
    }
    close(TOUCH);
    unlink($touchfile);

    #   Directory is writable, create SQL record
    my $sql = "INSERT INTO $opts{runs_table} " .
        "(centerid,dirname,comments,count,dateinit) " .
        "VALUES($cid,'$d','',0,'$nowdate')";
    my $sth = DoSQLv($sql);
    my $runid = $sth->{mysql_insertid};
    print "$Script - Added run/project '$d' [ $runid ]\n";
    $opts{runcount}++;
    return $runid;
}

#==================================================================
# Subroutine:
#   AddBams - Add details for bams in this directory to the database
#       In order to know if we've found all the bam files, we look
#       for any BAMs older than the MD5 files. If so, we parse the
#       md5 files again.
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
        print "$Script - Unable to read directory '$d': $!\n";
        return 0;
    }
   	if ($opts{verbose}) { warn "$Script - AddBams($runid, $d)\n"; }

    #   Get all the known bams for this directory/run   Get NWD name and original name
    my $sql = "SELECT bamname,bamname_orig FROM $opts{samples_table} WHERE $opts{samples_pkey}=$runid";
    my $sth = DoSQL($sql);
    my $rowsofdata = $sth->rows();
    my %knownbams = ();
    my %knownorigbam = ();
    for (my $i=1; $i<=$rowsofdata; $i++) {
        foreach my $href ($sth->fetchrow_hashref) {
            $knownbams{$href->{bamname}} = 1;
            if (exists($href->{bamname_orig}) && defined($href->{bamname_orig}) &&
                $href->{bamname_orig} ne '' && $href->{bamname_orig} ne ' ') {
                $knownorigbam{$href->{bamname_orig}} = 1;
            }
        }
    }

    #   Get list of all MD5 files and reprocess them. Don't get excited when
    #   the BAM has already been added to bamfiles. This will not go on forever.
    my @md5files = ();
    if ( -r "$d/Manifest.txt") { push @md5files,'Manifest.txt'; }
    else {
        opendir(my $dh, $d) ||
            die "$Script - Unable to read directory '$d'\n";
        while (readdir $dh) {
            if (! /\.md5$/) { next; }
            push @md5files,$_;
        }
        closedir $dh;
    }
    if (! @md5files) {
        if ($opts{verbose}) { print "$Script - No new MD5 files found in $d\n"; }
        return 0;
    }

    #   There is no consistency what people do here.
    #   Foreach md5 file, get the BAM name and create the bamfiles record
    #   and rename the md5 file so we do not process it again
    my $newbams = 0;
    my ($fn, $checksum);
    foreach my $origf (@md5files) {
        my $f = "$d/$origf";
        if (! open(IN, $f)) {
            warn "$Script - Unable to read MD5 file '$f': $!\n";
            next;
        }
        my $badmd5 = 0;
        while (my $l = <IN>) {          # Read md5 checksum file
            ($fn, $checksum) = NormalizeMD5Line($l, $f);
            if (! $checksum) { $badmd5++; next; }
            #   Don't make duplicate records
            if (exists($knownbams{$fn})) { next; }       # Skip known bams
            if (exists($knownorigbam{$fn})) { next; }    # Skip known original bams

            #   New BAM, create database record for it. Actual BAM may not exist yet
            $sql = "INSERT INTO $opts{samples_table} " .
                "(runid,bamname,checksum,dateinit) " .
                "VALUES($runid,'$fn','$checksum', $nowdate)";
            if ($opts{verbose}) { print "SQL=$sql\n"; }
            $sth = DoSQLv($sql);
            $newbams++;
        }
        close(IN);
        if ($badmd5) { next; }              # This MD5 was in error
    }

    #   If we added bams, change the count
    if ($newbams) {
        print "$Script - $newbams new bams found in '$d'\n";
        $opts{count} += $newbams;            # Stats for ending msg
        $opts{countruns} .= $d . ' ';
    }

    #   Get number of database records
   	#	We no longer have data arriving piecemeal - set arrived
    $sql = "SELECT $opts{samples_pkey} FROM $opts{samples_table} WHERE $opts{runs_pkey}=$runid";
    $sth = DoSQL($sql);
    my $numbamrecords = $sth->rows();
    $sql = "UPDATE $opts{runs_table}  SET arrived='Y',count=$numbamrecords WHERE runid=$runid";
    $sth = DoSQLv($sql);

    #   Last sanity check, see if number of BAM files matches records
    #   This might not always be accurate, but ...
    if ($newbams) {
        my $n = `ls $d/N*.bam $d/N*.cram 2>/dev/null | wc -l`;
        chomp($n);
        if ($n eq $numbamrecords) { print "$Script - Congratulations, # bams = # database records\n"; }
        else {
            print "$Script - Warning, # bams [$n] != # database records [$numbamrecords].  " .
            "If data is incoming, this might be OK\n";
        }
    }
    return 1;								# No need to check/set arrived

    #   If ALL the samples has been processed as arrived, maybe we do not
    #   need to look at this run any more.
    $sql = "SELECT $opts{samples_pkey} FROM $opts{samples_table} WHERE state_arrive!=$NOTSET";
    $sth = DoSQL($sql);
    my $numberarrived = $sth->rows();
    if (! $numberarrived) { return 1; }         # No sample processed, keep looking

    #   At least one sample was marked as arrived, see if there has been little changed
    #   in this directory in a long time, then make this run as 'arrived' so we'll
    #   look at this run again
    my $oldestbamdate = OldestBAM($d);
    if ((time() - $oldestbamdate) > $opts{arrivedsecs}) {
        $sql = "UPDATE $opts{runs_table}  SET arrived='Y' WHERE runid=$runid";
        $sth = DoSQLv($sql);
        print "$Script - Run '$d' has finally arrived. Look at it no more\n";
    }
    return 1;
}

#==================================================================
# Subroutine:
#   Get date of oldest BAM in this directory
#   filedate = OldestBAM($d)
#
# Arguments:
#   d - directory to search
#
# Returns:
#   Date of oldest BAM or zero
#==================================================================
sub OldestBAM {
    my ($d) = @_;
    my $oldestbamdate = 0;
    opendir(my $dh, $d) ||
        die "$Script - Unable to read directory '$d'\n";
    while (readdir $dh) {
        if ((! /\.bam$/) && (! /\.cram$/)) { next; }
        my @stats = stat("$d/$_");
        if ($oldestbamdate < $stats[9]) { $oldestbamdate = $stats[9]; }
    }
    closedir $dh;
    return $oldestbamdate;
}

#==================================================================
# Subroutine:
#   NormalizeMD5Line - Return filename and checksum
#
# Arguments:
#   l - line from an md5 file (well, what they pretend is an md5 file)
#   f - name of file (for error msgs only)
#	nocheck - boolean to avoid checking filetype
#
# Returns:
#   array of filename and checksum or NULL
#==================================================================
sub NormalizeMD5Line {
    my ($l, $f, $nocheck) = @_;
    if (! defined($nocheck)) { $nocheck = 0; }
    my @retvals = ();
    if ($l =~ /^#/) { return @retvals; }

    my ($checksum, $fn) = split(' ',$l);
    #   Do sanity check trying to guess the format of their md5 file
    if ($checksum =~ /\./) { ($fn, $checksum) = ($checksum, $fn); }
    if (length($checksum) != 32) {
        print "$Script - Invalid checksum '$checksum' in '$f'. Line: $l";
        return @retvals;
    }
    if ($fn =~ /\//) { $fn = basename($fn); }   # Path should not be qualified, but ...
    if ((! $nocheck) && $fn !~ /bam$/ && $fn !~ /cram$/) {	# File better be a bam or cram
        #   Only generate error message for true errors, not dumb ones :-(
        if ($fn !~ /bai$/) { print "$Script - Invalid BAM name '$fn' in '$f'. Line: $l"; }
        return @retvals;
    }
    @retvals = ($fn, $checksum);
    return @retvals;
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
#   DoSQLv - Run SQL with verbose check
#
# Arguments:
#   sql
#
# Returns:
#   boolean if it worked
#==================================================================
sub DoSQLv {
    my ($sql) = @_;
    if ($opts{verbose}) {
    	print "SQL: $sql\n";
    	return 0;
    }
    return DoSQL($sql);
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_init.pl - check for files that are arriving and initialize the database

=head1 SYNOPSIS

  topmed_init.pl updatedb

=head1 DESCRIPTION

This program monitors directories for incoming data and then
creates records for new run directories and/or new BAMs
(based on new MD5 files being discovered).

=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name on which to run the action, e.g. B<uw>.
This is useful for testing.
The default is to run against all centers.

=item B<-datatype TYPE>

Specifies the type of data to check.
The default is B<genome>.

=item B<-help>

Generates this output.

=item B<-ignorearrived>

Specifies we should re-look at the run because something new has arrived.

=item B<-project PROJECT>

Specifies these commands are to be used for a specific project.
Warning, this can only be abbreviated as B<-p> or <-project>.
The default is to use the environment variable PROJECT.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<topmed>.

=item B<-runs NAME>

Specifies a specific run on which to run the action.
This is useful for testing.
The default is to run against all runs for the center.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<updatedb>

This directs the program to monitor the B<-dir> directory for changes
and update the database table specified by B<-realm>.

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

