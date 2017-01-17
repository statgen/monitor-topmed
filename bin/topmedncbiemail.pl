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
use My_DB;
use POSIX qw(strftime tmpnam);

my $NOTSET = 0;                     # Not set
my $REQUESTED = 1;                  # Task requested
my $SUBMITTED = 2;                  # Task submitted to be run
my $STARTED   = 3;                  # Task started
my $DELIVERED = 19;                 # Data delivered, but not confirmed
my $COMPLETED = 20;                 # Task completed successfully
my $CANCELLED = 89;                 # Task cancelled
my $FAILED    = 99;                 # Task failed

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    logfile => '/net/topmed/working/topmed-output/XMLfiles/ncbiemail.log.txt',
    maildir2 => '/net/topmed/working/topmed-output/ncbimail',
    maildir => '/net/dumbo/home/tpg/ncbimail',
    mailext => '.mail',                 # Mail files end with this
    realm => '/usr/cluster/monitor/etc/.db_connections/topmed',
    bamfiles_table => 'bamfiles',
    keep => 0,                          # Do not delete input file
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help logfile=s maildir=s clean keep verbose 
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options] parsemail";
    warn "$Script [options] parsemail\n" .
        "or\n" .
        "$Script [options] summarizemail\n" .
        "or\n" .
        "$Script [options] checkmail\n" .
        "\n" .
        "$Script [options] clean\n" .
        "\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;
$opts{logfile} = abs_path($opts{logfile});  # Don't get fooled by relative paths
$opts{maildir} = abs_path($opts{maildir});

my $nowdate = strftime('%Y/%m/%d %H:%M', localtime);

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn =~ /^par/)     { Parse($opts{maildir}, $opts{logfile}); exit; }
if ($fcn =~ /^sum/)     { Summary($opts{logfile}); exit; }
if ($fcn =~ /^check/)   { Check($opts{logfile}); exit; }

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
    print "$nowdate Parsing mail in $dir\n";
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
#   Summary($log)
#
#   Summarize the log file of Email from NCBI
#==================================================================
sub Summary {
    my ($log) = @_;
    my %count = ();
    my $in;

    print "$nowdate Summary of $log\n";

    open($in, $log) ||
        die "$Script - Unable to open file '$log': $!\n";
    while (my $l = <$in>) {
        if ($l !~ /^(NWD\d+)\.(\S+)\s+(.+)\s+20\d\d/) { next; }
        my $nwdid = $1;
        my $ext = $2;
        my $msg = $3;
        $count{$ext}++;
        my $s = $ext . ' ' . $msg;
        $count{$s}++;
    }
    close($in);

    #   Show summary
    foreach my $k (keys %count) {
        print "$k   $count{$k}\n";
    }
}
    
#==================================================================
# Subroutine:
#   Check($log)
#
#   Parse log file and check consistency of NCBI states indicated
#   by the log file with the monitor database
#==================================================================
sub Check {
    my ($log) = @_;
    my $in;
    my $needscorrection = 0;
    my $correct = 0;
    my $duplicates = 0;
    DBConnect($opts{realm});        # Open database

    print "$nowdate Checking states of samples in database against $log\n";

    my %alreadyseen = ();
    my $cleanlog = '';
    open($in, $log) ||
        die "$Script - Unable to open file '$log': $!\n";
    while (my $l = <$in>) {
        if ($l !~ /^(NWD\d+)\.(\S+)\s+(.+)\s+20\d\d/) { next; }
        my $nwdid = $1;
        my $ext = $2;
        my $msg = $3;
        #   Figure out what database column is to be checked
        my $dbcol = '';
        if ($ext =~ /expt.submit/)      { $dbcol = 'state_ncbiexpt'; }
        if ($ext =~ /secondary.submit/) { $dbcol = 'state_ncbiorig'; }
        if ($ext =~ /remap.submit/)     { $dbcol = 'state_ncbib37'; }
        if (! $dbcol) {
            print "Unable to determine a column for '$ext'. Line=$l";
            next;
        }
        my $sth = ExecSQL("SELECT $dbcol FROM $opts{bamfiles_table} WHERE expt_sampleid='$nwdid'");
        my $rowsofdata = $sth->rows();
        if (! $rowsofdata) {
            print "Unable to find '$nwdid' in the monitor database\n";
            next;
        }
        my $k = $nwdid . '.' . $ext . ' ' . $msg;
        if (exists($alreadyseen{$k})) {
            $duplicates++;
            if ($opts{verbose}) { print "   Ignored duplicate: $k\n"; }
            next;
        }
        $alreadyseen{$k}++;
        if ($opts{clean}) { $cleanlog .= $l; }      # Save if rewriting log file
        
        my $href = $sth->fetchrow_hashref;
        #   Now based on the NCBI Email text, see if our database agrees
        if ($msg =~ /has been released/) {  # Informational only, nothing done or to do
            $correct++;
            next;
        }

        if ($msg =~ /has errors/) {
            if ($href->{$dbcol} ne $FAILED) {
                print "NWDID=$nwdid has errors, but $dbcol not correct ($href->{$dbcol})\n";
                $needscorrection++;
            }
            else { $correct++; }
            next;
        }

        if ($msg =~ /waiting for runs/) {
            if ($href->{$dbcol} ne $COMPLETED) {
                print "NWDID=$nwdid is waiting for runs, but $dbcol not correct ($href->{$dbcol})\n";
                $needscorrection++;
            }
            else { $correct++; }
            next;
        }
        print "Unable to check msg=$msg  Line=$l";
    }
    close($in);

    #   Maybe rewrite the log file
    if ($cleanlog) {
        my $f = $log . '.clean';
        my $out;
        open($out, '>' . $f) ||
            die "$Script - Unable to rewrite file '$f': $!\n";
        print $out $cleanlog;
        close($out);
        rename($log, "$log.prev") ||
            die "$Script - Unable to rename log file: $log\n";
        rename($f, $log) ||
            die "$Script - Unable to rename $f to $log\n";
        print "Log file '$log' rewritten to remove duplicates\n";
    }

    print "Monitor database fields that are correct=$correct\n";
    print "Monitor database fields that are incorrect=$needscorrection\n";
    print "Ignored $duplicates duplicate Email entries\n";
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
# Subroutine:
#   ExecSQL($sql, $die)
#
#   Execute SQL.  Keep trying if connection lost.
#
#   Returns handle for SQL
#==================================================================
sub ExecSQL {
    my ($sql, $die) = @_;
    if ($opts{persist}) { return PersistDoSQL($opts{realm}, $sql); }
    return DoSQL($sql, $die);
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

