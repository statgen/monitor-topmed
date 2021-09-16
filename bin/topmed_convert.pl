#!/usr/bin/perl
###################################################################
#
# Name:	topmed_convert.pl
#
# Description:
#   Use this program to rewrite RNA seq summary data files.
#	Given a map file it can convert NWGC numeric IDs into 
#	TOR values as well as drop 'unknown' TOR columns.
#
# ChangeLog:
#   $Log: topmed_convert.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use FindBin qw($Script);
use File::Spec;
use Getopt::Long;
use File::Basename;
#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    prefix => 'id2tor',
    verbose => 0,
);
Getopt::Long::GetOptions( \%opts,qw(
    help verbose prefix=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
	warn "$Script [options] mapfile file1 ... fileN\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $mapfile = shift @ARGV;
my ($IN, $OUT);
my %id2tor = ();				# Global Map of id to TOR

#	First parameter is a two column tab delimited map file providing a list of TOR IDs 
#	of interest.  This can also be used to map decimal IDs to their proper TOR.
LoadMapFile($mapfile);
Convert($opts{prefix}, \@ARGV);
exit;

#==================================================================
# Subroutine:
#   Convert
#
# Arguments:
#	$prefix - Prefix rewrite of file with this + '.'
#   reffiles - reference to array of paths to data files to process
#
# Returns:
#   Boolean if anything failed
#==================================================================
sub Convert {
    my ($prefix, $reffiles) = @_;
    if (! $reffiles) {
        die "$Script - No input files provided\n";
    }
	#------------------------------------------------------------------
    #	Parse each input file. These may have 1 or 3 header lines
    #	If necessary convert ID to TOR
    #	Drop any columns that do not have our TOR data
    #	Write new file to current working directory
	#------------------------------------------------------------------
	my @keepcol;	 	# List of columns to keep
	my @newids;
	my $l;
	foreach my $f (@$reffiles) {
		if ($f =~ /.gz$/) {
			open($IN, "gunzip -c $f |") ||
    			die "$Script - Unable to read file '$f: $!\n";
    	}
    	else {
    		open($IN, $f) ||
				die "$Script - Unable to read file '$f': $!\n";
    	}
		if ($opts{verbose}) { print "  Reading file '$f'\n"; }

		# 	Figure out what kind of file this is, headers can differ
		#	*.gct files look like
		#	#1.2
		#	NROWs  NIDs     (e.g. number of data rows and number of IDs)
		#	Name	Description	ID1 etc
		#
		#	Non GCT files just look like
		#	gene_id  transcript_id(s)  ID1 etc
	 	#	  (with a possibly empty header for the first column)
		$l = <$IN>;
		$l =~ s/\r//g;
		chomp($l);
		my $hdr = '';			# If GCT (3 line header), this is non-blank
		if ($l =~ /^#/) {		# Is this a GCT header ?
			$hdr = $l . "\n";
			$l = <$IN>;			# Second line of GCT
			if ($l !~ /^(\d+)/) { die "  Unexpected header in '$f':\n$l\n"; }
			$hdr .= $1 . "\t";	# Must append number of ID columns later
			$l = <$IN>;			# Get third line of file 
			$l =~ s/\r//g;
			chomp($l);
		}

		#	Now have first line of data with column headers and list of IDs/TORs
		#	If this is a TOR we want (e.g. from map file) keep this column
		#	Looks like:  '	someword	218260	218261	218262 ...'
		#	  or
		#	Looks like:  'Name Description	218260	218261	218262 ...'
		my @c = split("\t", $l);
		#	Always keep column zero and maybe column 1
		@keepcol = (0);
		@newids = ($c[0]);			# Build new line by column here
		my $starti = 1;
		if ($c[$starti] !~ /^\d{6}/ && $c[$starti] !~ /^TOR/) {  # Not ID/TOR
			push @keepcol,$starti;
			push @newids, $c[$starti];
			$starti++;
		}

		# 	Look at each ID/TOR column in this header line
		#	Ignore those we do not want (e.g. not in mapfile)
		my $countdropped = 0;
		my $countkept = 0;	
		for (my $i=$starti; $i<=$#c; $i++) {
			if ($id2tor{$c[$i]}) {		# If this is ID/TOR we want
				$c[$i] = $id2tor{$c[$i]};	# Convert possible numeric to TOR
				$countkept++;
				push @newids, $c[$i];		# Keep header entry
				push @keepcol, $i;		# Remember to keep data for this column
				if ($opts{verbose}) { print "    Keep column '$c[$i]' (col=$i)\n"; }
			}
			else {
				$countdropped++;
				if ($opts{verbose}) { print "    Drop ID $c[$i] (col=$i)\n"; }
			}
		}
		print "  Dropped $countdropped columns, keeping $countkept columns\n";

		if (! $countkept) { die "$Script - No columns from '$mapfile' found in '$f'\n"; }

		#	Prepare to generate output file in CWD
		#	If original input file was compressed, compress this one
		my ($v,$d,$filebit) = File::Spec->splitpath($f);
		my $ofile =  $prefix . '.' . $filebit;
		if ($f =~ /.gz$/) {
			open ($OUT, '| gzip >' . $ofile) ||
				die "$Script - Unable to create '$ofile': $!\n";
		}
		else {
			open ($OUT, '>' . $ofile) ||
				die "$Script - Unable to create '$ofile': $!\n";
		}
		if ($opts{verbose}) { print "  Creating '$ofile' ...\n"; }

		# 	If GCT file, we have a header to complete
		if ($hdr) { print $OUT $hdr . $countkept . "\n"; }
		print $OUT join("\t", @newids) . "\n";

		#	Read body of input data, drop columns not in our map
		while ($l = <$IN>) {
			$l =~ s/\r//g;
			chomp($l);
			my @c = split("\t", $l);
			@newids = ();						# Build new line
			for (my $i=0; $i<=$#keepcol; $i++) {
				if (defined($c[$keepcol[$i]])) { push @newids,$c[$keepcol[$i]]; }
			}
			print $OUT join("\t", @newids) . "\n";
		}
		# 	Output file created
		close($OUT);
		print "  Created '$ofile\n";
	}

	return;
}

#==================================================================
# Subroutine:
#   LoadMapFile - Initialize %id2tor from a tab delimited map file
#
#	Format of file:
#	   header which is ignored
#	   ID   TOR
#
#	If ID is a digital string, then TOR must be a TOR id for this ID
#   To make testing easier, if there is only one column, assume ID = TOR
#
# Arguments:
#   mapfile - path to file containing IDs and TORs
#
#==================================================================
sub LoadMapFile {
    my ($mapfile) = @_;

	open($IN, $mapfile) ||
		die "$Script - Unable to read map file '$mapfile': $!\n";
	<$IN>;						# Skip header
	my $k = 0;
	while (my $l = <$IN>) {
		my @c = split(' ', $l);
		if ($c[1]) { $id2tor{$c[0]} = $c[1]; }
		else { $id2tor{$c[0]} = $c[0]; }
		$k++;
	}
	close($IN);
	if (! $k) { die "  Unable to find any ID to TOR map data\n"; }
	print "  Loaded $k ID/TOR map entries from $mapfile\n";
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_convert.pl - Rewrite RNA seq files

=head1 SYNOPSIS

  topmed_convert.pl  mapfile file1 file2 ...

=head1 DESCRIPTION

This program reads a list of ID/TOR files of interest and converts
a new file of RNA seq data containing only those data in the list.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-prefix STRING>

The new RNA data file is created in the current working directory 
and is prefixed with 'string'.  This defaults to B<id2tor.>.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

=over 4

=item B<mapfile file1 file2 ...> 

Reads the list of ID/TOR names in B<mapfile> and then for each file provided
writes a new RNA seq file with only data for those ID/TORs.
The new files are prefixed with B<-prefix> and written to the current working directory.

=back

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0.  Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2021 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

