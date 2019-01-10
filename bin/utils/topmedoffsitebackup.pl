#!/usr/bin/perl
###################################################################
#
# Name: topmedoffsitebackup.pl
#
# Description:
#   Use this program to figure out if the backup directory link
#   for a run that is backed up offsite is correct.
#   If the incoming cram files have been backed up offsite, then
#   the backup directory link should point to the the incoming files.
#   This is needed so that when we send data to NCBI, we can find
#   the crams to be sent.
#
# ChangeLog:
#   $Log: topmedoffsitebackup.pl,v $
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

use Getopt::Long;
use File::Basename;
use Topmed::Get;
use My_DB;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    topmedpath => '/usr/cluster/topmed/bin/topmedpath.pl -p topmed',
    centers_table => 'centers',
    bamfiles_table => 'bamfiles',
    runs_table => 'runs',
	tmpfile => "/tmp/$Script.temp.sh",
    topdir => "/net/topmed/working/backups/incoming/topmed",
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose
)) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    warn "$Script [options] check\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;
if ($fcn ne 'check') { die "$Script - Unknown directive '$fcn'. Try $Script --help\n"; }

DBConnect($opts{realm});

#--------------------------------------------------------------
#   Look for runs without a proper backup directory
#--------------------------------------------------------------
#	Create bash for loop in shell script
my $out;
open($out, '>' . $opts{tmpfile}) ||
	die "$Script - Unable to create '$opts{tmpfile})': $!\n";
print $out "#!/bin/bash\n" .
	"#  cd $opts{topdir} || exit 5\n" .
	"for c in \$1; do\n" .
	"  cd $opts{topdir}\n" .
	"  cd \$c\n " .
	"    for r in *; do\n" .
	"       n=`ls \$r/*.src.cram | wc -l`\n" .
	"       cd \$r\n" .
	"       echo \$c/\$r  \$n  `pwd -P`\n" .
	"       cd ..\n" .
	"     done\n" .
	"done\n";
close($out);

#   Get all the known centers in the database
my $centersref = GetCenters();
foreach my $cid (keys %{$centersref}) {
	my $runs = 0;
	my $errs = 0;
    my $center = $centersref->{$cid};
	print "Checking $center\n";
	#   Get complete listing of all local backup directories 
	#	and counts of src.cram files for this center.
	my %backupruns = ();
	my $cmd = "bash $opts{tmpfile} $center";
	my $in;
	open ($in, "$cmd |") ||
		die "$Script - Unable to run command: $cmd\n";
	while (<$in>) {
		chomp;
		my @cols = split(' ', $_);
		if ($cols[1] eq '0') {
			warn "$Script - run $cols[0] found zero backup files\n";
			$errs++;
			next;
		}
		$backupruns{$cols[0]} = $cols[2];
	}
	close($in);

	#	For each run in this center, see if the local backup directory
	#	has a src.cram for one sample and see if topmedpath returns
	#	the correct place for this sample
    my $runsref = GetRuns($cid);
    if (! $runsref) {
    	warn "$Script - found no runs in database for center 'center'\n";
    	$errs++;
    	next;
    }
    foreach my $runid (keys %{$runsref}) {
        my $dirname = $runsref->{$runid};
        my $centerrun = "$center/$dirname";
        my $sql = "SELECT offsite FROM $opts{runs_table} WHERE runid=$runid";
        my $sth = DoSQL($sql);
        my $href = $sth->fetchrow_hashref;
        my $offsite = $href->{offsite};
        $runs++;
        #	Get one sample for sanity checking
        $sql = "SELECT bamid,cramname FROM " . $opts{bamfiles_table} . 
            " WHERE runid=$runid LIMIT 1";
        $sth = DoSQL($sql);
        if (! $sth) {
    		warn "$Script - no sample for $centerrun\n";
    		$errs++;
    		next;
    	}
        $href = $sth->fetchrow_hashref;
        if (! $href) {
    		warn "$Script - no sample for $centerrun\n";
    		$errs++;
    		next;
    	}
    	my $bamid = $href->{bamid};
    	my $cramname = $href->{cramname};

		# 	Verify this run has a backup directory
		if (! exists($backupruns{$centerrun})) {
    		warn "$Script - no backup directory found for for $centerrun\n";
    		$errs++;
    		next;
    	}
    	#	Where is this backup directory supposed to be
    	#	Does it match what we have in $opts{topdir} tree?
    	my $formalpath = `$opts{topmedpath} wherepath $bamid localbackup`;
    	chomp($formalpath);
    	if ($formalpath ne $backupruns{$centerrun}) {
    		warn "$Script - topmedpath wherepath does not match tree in $opts{topdir}\n" .
    			"    $formalpath ne $backupruns{$centerrun}\n";
    		$errs++;
    		next;
    	}
		#	Now see if one of our samples is there
		if (! -f "$formalpath/$cramname") {
    		warn "$Script - Missing backup file $formalpath/$cramname\n";
    		$errs++;
    		next;
    	}
    }
	print "  Found $runs runs and detected $errs errors\n";
}
unlink($opts{tmpfile});
exit;



#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedoffsitebackup.pl - See if the backup directory link for a run is correct

=head1 SYNOPSIS

  topmedoffsitebackup.pl check

=head1 DESCRIPTION

Use this program to figure out if the backup directory link
for a run that is backed up offsite is correct.
If the incoming cram files have been backed up offsite, then
the backup directory link should point to the the incoming files.
This is needed so that when we send data to NCBI, we can find
the crams to be sent.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-realm NAME>

Specifies the realm name to be used.
This defaults to B<$opts{realm}> in the same directory as
where this program is to be found.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS


B<check>
Use this to look for crams that are not in the correct place.

=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2016, 2019 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

