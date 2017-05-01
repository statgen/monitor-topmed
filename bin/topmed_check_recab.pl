#!/usr/bin/perl
###################################################################
#
# Name: topmed_check_recab.pl
#
# Description:
#   Apparently samtools or something else in the remapping pipeline
#   can screw up the resulting cram.
#
#   The `SM:` field in the `@RG` headers should only ever reference
#   a single NWD. Anything more than one ID is an error state and
#   the remapping has failed. 
#
#   This pgm should be invoked in a bash shell like this:
#       set -o pipefail
#       samtools view -H <sample.cram> | grep '^@RG' | topmed_check_recab.pl
#       if [ "$?" != "0" ]; then
#           echo "Recab header was messed up"
#           exit 4
#       fi
#
#   Unfortunately, it's not a simple awk on the columns because it is
#   not always in the same position and I even found weird line breaks
#   within some sample headers (which means a double chomp() might be needed.)
#
# ChangeLog:
#   $Log: topmed_check_recab.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;

#--------------------------------------------------------------
#   Parse samtools header information for multiple sample ids
#--------------------------------------------------------------
my %id = ();
while (<>) {
  chomp;
  my $h = grep {/^SM:/} split(/\t/);
  $id{$h} = 1;
}       

my @samples = keys %id;
my $n = scalar(@samples);
if ($n == 1) { exit 0; }
print "Incorrect number of samples found: " . join(' ', @samples) . "\n";
exit(1);
