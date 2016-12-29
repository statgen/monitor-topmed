#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
###################################################################
#
# Name: topmedncbiemial.pl
#
# Description:
#   Use this program to parse mail from NCBI and save it
#
# ChangeLog:
#   $Log: topmedncbiemial.pl,v $
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
use Getopt::Long;
use Cwd qw(abs_path);

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    logfile => '/net/topmed/working/topmed-output/XMLfiles/ncbiemail.log.txt',
    maildir => '/net/topmed/working/topmed-output/ncbimail',
    mailext => '.txt',                  # Mail files end with this
    keep => 0,                          # Do not delete input file
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help logfile=s maildir=s keep verbose 
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options] parsemail";
    warn "$Script [options] parsemail\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;
$opts{logfile} = abs_path($opts{logfile});  # Don't get fooled by relative paths
$opts{maildir} = abs_path($opts{maildir});

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'parsemail') { Parse($opts{maildir}, $opts{logfile}); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   Parse($dir, $log)
#
#   Parse Email from NCBI
#==================================================================
sub Parse {
    my ($dir, $log) = @_;
    my %mon2num = (
        Jan => '01',
        Feb => '02',
        Mar => '03',
        Apr => '04',
        May => '05',
        Jun => '06',
        Jul => '07',
        Aug => '08',
        Sep => '09',
        Oct => '10',
        Nov => '11',
        Dec => '12'
    );
    #   Make sure the log file does not get deleted by find as an old file
    system("touch $log");

    #   Get list of files to parse
    opendir(my $dh, $dir) ||
        die "$Script - Unable to CD to directory of mail '$dir': $!\n";
    my @files = grep { -f "$dir/$_" } readdir($dh);
    closedir $dh;
    chdir($dir) ||
        die "$Script - Unable to CD to directory of mail '$dir': $!\n";

    #   Read all mail, extract lines of interest and save them
    #   Anything we cannot handle, we rename
    my $good = 0;
    my $notsogood = 0;
    my ($msgid, $logtext, $in, $date);
    foreach my $f (@files) {
        if ($f =~ /\.processed$/) { next; }     # We already did these
        if ($f =~ /\.unparsable$/) { next; }
        if ($f !~ /$opts{mailext}$/) {
            warn "$Script - Skipping file '$dir/$f', not Email\n";
            next;
        }
        if ($opts{verbose}) { print $f . "\n"; }
        $msgid = $logtext = $date ='';
        open($in, $f) ||
            die "$Script - Unable to read file '$dir/$f': $!\n";
        while(<$in>) {
            if (/^Date:\s+\w+, (\d+) (\w+) (\d+) (\S+)/) {      # Get date of email
                $date = $3 . '-' . $mon2num{$2} . '-' . $1 . ' ' . $4;
                next;
            }
            if (/^To: unattended/) {        # Start of interesting stuff
                while(<$in>) {
                    if (/^Message-Id: <(.+)>/) {
                        if ($msgid) { die "$Script - '$dir/$f' had two Message-Id lines\n"; }
                        $msgid = $1;
                        next;
                    }
                    if (/^This email has been sent/) {
                        s/\r//g;                        # Remove DOS crap
                        if ($logtext) {
                            die "$Script - '$dir/$f' has multiple inform lines:\n" .
                                "==>${logtext}==>$_";
                        }
                        #   Fix up the logtext making it shorter
                        $logtext = $_;
                        if ($logtext =~ /.+(NWD.+)/) {
                            $logtext = $1;
                            $logtext =~ s/\)//g;
                            $logtext =~ s/\(//g;
                            $logtext .= "\n";
                        }
                        next;
                    }
                }
            }
        }
        close($in);

        #   Maybe this is not a file we are interested in, save it
        if (! $msgid) {
            my $newf = $f . '.' . time() . '.unparsable';
            print "$Script - File '$dir/$newf' is not parsable\n";
            rename($f, $newf) ||
                die "$Script - Unable to rename '$dir/$newf' as unparsable file\n";
            $notsogood++;
            next;
        }
        #   Maybe this is not a file with an inform message, save it
        if (! $logtext) {
            my $newf = $f . '.' . time() . '.nomsg';
            print "$Script - File '$dir/$newf' had no msg for us\n";
            rename($f, $newf) ||
                die "$Script - Unable to rename '$dir/$newf' as nomsg file\n";
            $notsogood++;
            next;
        }
        #   Yes, we found something of interest
        #   Fix up the logtext making it shorter
        if ($logtext =~ /.+(NWD.+)/) {
            $logtext = $1;
            $logtext =~ s/\)//g;
            $logtext =~ s/\(//g;
            $logtext .= "\n";
        }
        if ($date) { chomp($logtext); $logtext .= ' ' . $date . "\n"; }
        Append($log, $logtext);
        $good++;
        if (! $opts{keep}) {                # Do not keep, rename Email file
            my $newf = $f . '.' . time() . '.processed';
            rename($f, $newf) ||
                die "$Script - Unable to rename '$dir/$newf' as processed file\n";
        }
    }
    print "Parsed $good files, found $notsogood files with issues\n";
    if ($good) { print "Appended to log $log\n"; }
}

#==================================================================
# Subroutine:
#   Append($file, $txt)
#
#   Append $txt to $file
#==================================================================
sub Append {
    my ($file, $txt) = @_;
    open(OUT, '>>' . $file) ||
        die "$Script - Unable to append to '$file': $!\n";
    print OUT $txt;
    close(OUT);
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedncbiemail.pl - Extract information from NCBI Emails

=head1 SYNOPSIS
 
  topmedncbiemail.pl parsemail

=head1 DESCRIPTION

This program reads an Email from NCBI and extracts the key information
from the Email into a log file for easy checking of the status of an NWDID.

=head1 OPTIONS

=over 4

=item B<-help>

Generates this output.

=item B<-keep>

Do not rename the Email file, but leave it alone.
The default is to rename the Email file with an extension '.processed'
which can be deleted by a cron job.

=item B<-logfile file>

Append the important Email information to a log file.
This defaults to B</net/topmed/working/topmed-output/XMLfiles/ncbiemail.log.txt>

=item B<-maildir file>

Search for B<*.txt> files in this directory.
This defaults to B</net/topmed/working/topmed-output/ncbimail>

=item B<-verbose>

Provided for developers to see additional information.

=back


=head1 PARAMETERS

=over 4

=item B<parsemail>

Will parse all Email files extracting the information wanted.

=back


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

