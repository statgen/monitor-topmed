#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl -I/usr/cluster/monitor/lib/perl5 -I /usr/cluster/monitor/bin
###################################################################
#
# Name: topmedaspera.pl
#
# Description:
#   Use this program to fetch files from a remote center.
#   The approach is to fetch the Manifest.txt file and
#   then fetch each file individually from the remote site.
#   This has the advantage that we can fetch data concurrently
#   and once the ascp command finishes, we know the file
#   has arrived and will not be changed.
#
# ChangeLog:
#   $Log: topmedaspera.pl,v $
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
use Cwd 'abs_path';

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    url => 'https://shares.broadinstitute.org/',
    center => 'broad',
    database => 'dbfile',                       # As found in $basedir 
    logdir => "/net/topmed/working/topmed-output/aspera",
    keyfile => 'keyfile',                       # Name of keyfile in $basedir
    hooksfile => '/tmp/ascp-shares.sh.hooks',   # File created by our ascpshares
    ascpshares => '/usr/cluster/monitor/bin/ascp-shares.sh',
    asperabroad => '/usr/cluster/monitor/bin/aspera.broad.sh',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help realm=s verbose url=s logdir=s
    user=s password=s decryptpassword=s
    target=s rmtdir=s dest=s
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 1 || $opts{help}) {
    warn "$Script [options] -user SN0096478 -password FA3XCLEI7JZAMX8 -decryptpassword yk2iCI7zefH3nmu -rmtdir SN0096478 manifest targetdir\n" .
        "  or\n" .
       "$Script [options] -dest 2016jun23 fetch file\n" .
       "\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

#--------------------------------------------------------------
#   Fetch data for a run
#--------------------------------------------------------------
#   First time, you must get the Manifest, then do fetch
if ($fcn eq 'manifest') {      # Start getting the manifest. Must be like SN0096478/Manifest
    if ((! exists($opts{user})) ||
        (! exists($opts{password})) ||
        (! exists($opts{rmtdir})) ||
        (! exists($opts{decryptpassword})) ) {
        die "$Script important options were not specified.  Try $Script --help\n";
    }

    if ($ARGV[0] !~ /^\//) {
        die "$Script - file to fetch ($ARGV[0]) may not begin with a '/'\n";
    }
    GetManifest($ARGV[0]);
    exit;
}

#   Fetch a file.  If Manifest was done, the dbfile was created with all the important details
if ($fcn eq 'fetch') {      # Start getting the manifest. Must be like SN0096478/Manifest
    if (! exists($opts{dest})) {
        die "$Script important options were not specified.  Try $Script --help\n";
    }
    FetchFile($ARGV[0]);
    exit;
}

die "$Script - Incorrect directive '$fcn', try $Script --help\n";
exit;

#==================================================================
# Subroutine:
#   GetManifest($tdir)
#
#   Fetch the Manifest file. This has the side effect of
#   saving key information for subsequent data being copied to the target.
#   We invoke ascp with an exclude of a single character file/directory
#   and this gets one or more files, possibly called Manifest.txt.
#   If that works, we can get a list of all files to be fetched.
#==================================================================
sub GetManifest {
    my ($tdir) = @_;

    # Figure out where these files should go
    if (! -d $tdir) {
        die "$Script - Directory '$tdir' does not exist\n";
    }
    my $fqtdir = abs_path($tdir);           # All our data goes here
    system("mkdir -p $fqtdir") && die "$Script - unable to create directory '$fqtdir'\n";
    system("mkdir -p $opts{logdir}") && die "$Script - unable to create directory '$opts{logdir}'\n";

    #   Squash the target directory so we have a unique identifier
    my $base = `basename $fqtdir`;
    chomp($base);
    $base =~ s/\//_/g;
    $base =~ s/ /_/g;
    my $basedir = $opts{logdir} . '/' . $base;        # Save fetch details here
    system("mkdir -p $basedir") && die "$Script - unable to create directory '$basedir\n";
    print "Details for '$tdir' will be saved in '$basedir'\n";

    #   Now see if we can get something like a Manifest
    $ENV{ASPERA_SCP_FILEPASS} = $opts{decryptpassword};
    my $cmd = "$opts{ascpshares} -hook -L $basedir $opts{url} $opts{user}:$opts{password} /$opts{rmtdir} $basedir/Manifest";
    if ($opts{verbose}) { print "Executing: $cmd\n"; }
    system($cmd) &&
        die "$Script - Unable to fetch manifest file\n";
    system("ls -l $basedir/Manifest");          # Show file we fetched 

    #   Here is data to be saved.  Assumes destination looks like .../center/dirname
    my $dirname = `basename $fqtdir`;
    chomp($dirname);
    my $center = `dirname $fqtdir`;
    chomp($center);
    $center = `basename $center`;
    chomp($center);
    my %dbdata = (                      # Save key information
        sharedir => $opts{rmtdir},
        destdir => $fqtdir,
        dirname => $dirname,
        center => $center,
        user => $opts{user},
        password => $opts{password},
        decryptpassword => $opts{decryptpassword},
        url => $opts{url},
    );

    #   Get hooks provided by ascpshares
    open(IN,$opts{hooksfile}) ||
        die "$Script - Unable to read hooks file from '$opts{ascpshares} -hook': $!\n";
    while (<IN>) {
        chomp();
        if (/(\S+)=(.+)/) { $dbdata{$1} = $2; }
    }
    close(IN);
    #   Now save those hooks from ascpshares
    if (! exists($dbdata{keyfile})) { die "$Script - missing 'keyfile' in $opts{hooksfile}\n"; }
    if (! exists($dbdata{basecmd})) { die "$Script - missing 'basecmd' in $opts{hooksfile}\n"; }
    system("mv $dbdata{keyfile} $basedir/$opts{keyfile}") &&
        die "$Script - Unable to save keyfile from '$opts{ascpshares} -hook': $!\n";
    $dbdata{keyfile} = "$basedir/$opts{keyfile}";   # Change where keyfile is now
    if ($dbdata{basecmd} !~ /^(.+) -i \S+/) {
        die "$Script - unable to parse ascpshares basecmd: $dbdata{basecmd}\n";
    }
    $dbdata{basecmd} = "$1 -i $basedir/$opts{keyfile} --file-crypt=decrypt";     # Correct command to fetch later
    unlink($opts{hooksfile});
    
    #   Create our master key file so we can fetch more files
    open(OUT, '>' . "$basedir/$opts{database}") ||
        die "$Script - Unable to create file '$basedir/$opts{database}': $!\n";
    foreach my $k (keys %dbdata) { print OUT "$k=$dbdata{$k}\n"; }
    close(OUT);
    print "Saved information for '$base' in $basedir\n";

    #   Now move the Manifest file to the proper place
    $_ = `ls $basedir/Manifest/Manifest*`;
    chomp();
    if ("$_" eq "") { die "$Script - Unable to retrieve '$basedir/Manifest/Manifest*'\n"; }

    my $manfile = "$fqtdir/Manifest.txt";
    system("mv -f $_ $manfile") &&
        die "$Script - Unable to move Manifest.txt to $fqtdir: $!\n";
    system("chmod 0444 $manfile") &&
        die "$Script - Unable to force Manifest.txt as readonly\n";
    rmdir("$basedir/Manifest");
    print "Saved Manifest in $fqtdir\n";

    #   Read the manifest and generate the commands to fetch the file
    my $fetchcmds = "$basedir/fetch-commands.txt";
    open(OUT, '>' . $fetchcmds) ||
        die "$Script - Unable to create file '$fetchcmds': $!\n";
    open(IN, $manfile) ||
        die "$Script - Unable to read file '$manfile': $!\n";
    my $n = 0;
    while (<IN>) {
        if (! /^(\S+)/) { next; }
        my $f = $1;
        $n++;
        print OUT "$0 -dest $dirname fetch $f\n";
    }
    close(IN);
    close(OUT);
    print "Created list of fetch commands in $fetchcmds\n";
}

#==================================================================
# Subroutine:
#   FetchFile($sendfile)
#
#   Fetch a file (maybe as dir/file) as a temporary.
#   When successful, rename it to what it should be.
#==================================================================
sub FetchFile {
    my ($sendfile) = @_;               

    #   Squash the target directory so we have a unique identifier
    my $base = $opts{dest};
    $base =~ s/\//_/g;
    $base =~ s/ /_/g;
    my $basedir = $opts{logdir} . '/' . $base;        # Saved fetch details here
    #   Read database data for this dest
    my $f = "$basedir/$opts{database}";
    open(IN,$f) ||
        die "$Script - Unable to read saved information for '$opts{dest}' from '$f': $!\n";
    my %dbdata = ();
    while (<IN>) {
        if (! /(\S+)=(.+)/) { die "$Script - Unable to parse line in '$f'\n  LINE=$_\n"; }
        $dbdata{$1} = $2;
    }
    close(IN);

    #   Get variables specific to this particular file
    my $rmtfile = "/$dbdata{sharedir}/$sendfile.aspera-env";
    my $cmd = "$opts{asperabroad} $dbdata{user}:$dbdata{password} $rmtfile";
    open(IN, "$cmd |") ||
        die "$Script - Unable to get crucial data for broad file: CMD=$cmd\n";
    my ($user, $host, $token, $srcfile) = split(' ', <IN>);
    close(IN);

    #   Build command to fetch file
    my $realfile = `basename $sendfile`;
    chomp($realfile);
    $ENV{ASPERA_SCP_FILEPASS} = $dbdata{decryptpassword};
    mkdir "$dbdata{destdir}/.$realfile";
    mkdir("$basedir/logs");                 # Aspera logs go here
    my $logdir = "$basedir/logs/$realfile";
    mkdir($logdir);
    $cmd = "$dbdata{basecmd} -W $token -L $logdir $user\@$host:$rmtfile $dbdata{destdir}/.$realfile";
    if ($opts{verbose}) { print "Executing: $cmd\n"; }
    system($cmd) &&
        die "$Script - Unable to fetch file '$sendfile\n";
    system("ls -l $dbdata{destdir}/.$realfile");     # Show tmp file we fetched 
    $cmd = "mv $dbdata{destdir}/.$realfile/$realfile $dbdata{destdir}";
    system($cmd) &&
        die "$Script - Unable to rename newly arrived file:  CMD=$cmd\n";
    system("chmod 0444 $dbdata{destdir}/$realfile");

    rmdir("$dbdata{destdir}/.$realfile");
    system("ls -l $dbdata{destdir}/$realfile");     # Show file we fetched 
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

  topmedcmd.pl show 2199 state_cram        # Show a column
  topmedcmd.pl show 2199 center            # Show center 
  topmedcmd.pl show NWD00234 run           # Show run
  topmedcmd.pl show NWD00234 yaml          # Show everything known about a bamid
 
  topmedcmd.pl where 2199 bam              # Returns real path to bam and host of bam
  topmedcmd.pl where 2199 backup           # Returns path to backups directory and to backup file and host
  topmedcmd.pl where NWD00234 b37          # Returns path to remapped b37 directory and to file
  topmedcmd.pl where 2199 b38              # Returns path to remapped b38 directory and to file

  topmedcmd.pl permit add bai braod 2015oct18   # Stop bai job submissions for a run
  topmedcmd.pl permit remove 12             # Remove a permit control
  topmedcmd.pl permit test bai 4567         # Test if we should submit a bai job for one bam

  topmedcmd.pl -maxlongjobs 35 squeue       # Show what is running
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

=item B<-maxlongjobs N>

Specifies how many of the longest running jobs to show.
This defaults to B<10>.

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

B<mark bamid|nwdid dirname  [verb] [state]>
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

B<set bamid|nwdid columnname value>
Use this to set the value for a column for a particular BAM file.

B<send2ncbi filelist>
Use this to copy data to NCBI with ascp.

B<show arrived>
Use this to show the bamids for all BAMs that are marked arrived.

B<show bamid|nwdid colname|center|run|yaml>
Use this to show information about a particular bamid (or expt_sampleid).
Use 'yaml' to display everything known about the bam of interest.

B<unmark bamid|nwdid [verb]>
Use this to reset the state for a particular BAM file to the default
database value.

B<whatnwdid bamid|nwdid>
Use this to get some details for a particular bam.

B<where bamid|nwdid bam|backup|b37|b38>
If B<bam> was specified, display the path to the real bam file, not one that is symlinked
and the host where the bam exists (or null string).
If B<backup> was specified, display the path to the backup directory 
and the path to the backup file (neither of which may not exist)
and the host where the backup BAM file should exist (or null string).
If B<b37> was specified, display the path to the directory of remapped data for build 37 (or 'none')
and the path to the remapped file (which may not exist).
If B<b38> was specified, display the path to the directory of remapped data for build 38 (or 'none')
and the path to the remapped file (which may not exist).

=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2015-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

