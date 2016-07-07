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
use File::Basename;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#   UW pw = gobuN4je (aspera) and Cah4Meu6 (aspera2)
#--------------------------------------------------------------
our %opts = (
    broadurl => 'https://shares.broadinstitute.org/',
    uwhost => 'aspera.gs.washington.edu',
    center => '',
    database => 'dbfile',                       # As found in $asperadir
    ascp => '/net/mario/cluster/bin/ascp -d -Q -l 1500M -k 1',
    logdir => "/net/topmed/working/topmed-output/aspera",
    keyfile => 'keyfile',                       # Name of keyfile in $asperadir
    hooksfile => '/tmp/ascp-shares.sh.hooks',   # File created by our ascpshares
    ascpshares => '/usr/cluster/monitor/bin/ascp-shares.sh',
    asperabroad => '/usr/cluster/monitor/bin/aspera.broad.sh',
    uwrmthost => 'aspera.gs.washington.edu',
    verbose => 0,
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose center=s
    user=s password=s decryptpassword=s
    target=s rmtdir=s dest=s
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 1 || $opts{help}) {
    warn "$Script [options] -center broad \\\n" .
        "  -user SN0096478 -password FA3XCLEI7JZAMX8 -decryptpassword yk2iCI7zefH3nmu \\\n" .
        "  -rmtdir SN0096478 \\\n" .
        "  manifest targetdir\n" .
        "  or\n" .
        "$Script [options] -center uw \\\n" .
        #  or --user nhlbi-aspera2 --password Cah4Meu6
        "  -user nhlbi-aspera -password gobuN4je \\\n" .
        "  -rmtdir /nhlbi_macrogen_2/arnett_nhlbi_wgs_hypergen_utah/batch_9 \\\n" .
        "  manifest targetdir\n" .
        "  or\n" .
       "$Script [options] -center broad -dest 2016jun23 fetch file\n" .
       "\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

#   Set the URL for this center, unless the user provided it
if (! $opts{center}) {
    die "$Script - Option -center was not provided\n";
}

#--------------------------------------------------------------
#   Fetch the list of files, the manifest
#--------------------------------------------------------------
if ($fcn eq 'manifest') {      # Start getting the manifest. e.g. SN0096478/Manifest
    if ($opts{center} eq 'broad') { GetManifestBroad($ARGV[0]); exit; }
    if ($opts{center} eq 'uw')    { GetManifestUW($ARGV[0]); exit; } 
    die "$Script - Unable to handle Manifest for center '$opts{center}'\n";
}

#--------------------------------------------------------------
#   Fetch a file.
#   If Manifest was done, the dbfile was created with all the important details
#--------------------------------------------------------------
if ($fcn eq 'fetch') {      # Start getting the manifest. Must be like SN0096478/Manifest
    if (! exists($opts{dest})) {
        die "$Script - option -dest was not specified.  Try $Script --help\n";
    }
    if ($opts{center} eq 'broad') { FetchFile($ARGV[0]); exit; }
    if ($opts{center} eq 'uw')    { FetchFile($ARGV[0]); exit; } 
}

die "$Script - Incorrect directive '$fcn', try $Script --help\n";
exit;

#==================================================================
# Subroutine:
#   GetManifestUW($targetdir)
#
#   Fetch the Manifest file. This has the side effect of
#   saving key information for subsequent data being copied to the target.
#   We invoke ascp with an exclude of a single character file/directory
#   and this gets one or more files, possibly called Manifest.txt.
#   If that works, we can get a list of all files to be fetched.
#==================================================================
sub GetManifestUW {
    my ($tdir) = @_;

    #   The broad requires a number of options
    if ((! exists($opts{user})) ||
        (! exists($opts{password})) ||
        (! exists($opts{rmtdir})) ) {
        die "$Script important $opts{center} options were not specified.  Try $Script --help\n";
    }

    # Figure out where these files should go
    if (! -d $tdir) {
        die "$Script - Directory '$tdir' does not exist\n";
    }
    my $runname = basename($tdir);
    my ($fqtdir, $asperadir) = MakeAsperaDir($tdir);

    #   Now see if we can get a mess of md5 files and create a Manifest
    $ENV{ASPERA_SCP_PASS} = $opts{password};
    my $basecmd = "$opts{ascp} --src-base=$opts{rmtdir} -L $asperadir";
    my $cmd = "$basecmd -E '*.bam' $opts{user}\@$opts{uwhost}:$opts{rmtdir} $fqtdir";
    if ($opts{verbose}) { print "Executing: $cmd\n"; }
    system($cmd) &&
        die "$Script - Unable to fetch manifest md5 files for '$tdir'\n  CMD=$cmd\n";
    $cmd = "cat $fqtdir/*.md5 > $fqtdir/Manifest.txt";
    system($cmd) &&
        die "$Script - Unable to create Manifest.txt from md5 files\n  CMD=$cmd\n";
    system("ls -l $fqtdir/Manifest.txt");          # Show file we ceated 
    print "Saved Manifest in $fqtdir\n";

    #   Save data for fetch
    my %datatosave = (
        dirname => $runname,
        center => 'uw',
        password => $opts{password},
        basecmd => $basecmd, 
        destdir => $fqtdir,
        user => $opts{user},
    );
    MakeAsperaData("$fqtdir/Manifest.txt", $asperadir, $runname, \%datatosave);

}

#==================================================================
# Subroutine:
#   GetManifestBroad($targetdir)
#
#   This assumes the $opts{rmtdir} minimally looks like:
#       Manifest.txt
#       x     (set of subdirectories where the data resides)
#       
#   We invoke ascp with an exclude of a single character file/directory
#   and this gets one or more files, possibly called Manifest.txt.
#   If that works, we can get a list of all files to be fetched.
#==================================================================
sub GetManifestBroad {
    my ($tdir) = @_;

    #   The broad requires a number of options
    if ((! exists($opts{user})) ||
        (! exists($opts{password})) ||
        (! exists($opts{rmtdir})) ||
        (! exists($opts{decryptpassword})) ) {
        die "$Script important $opts{center} options were not specified.  Try $Script --help\n";
    }

    # Figure out where these files should go
    if (! -d $tdir) {
        die "$Script - Directory '$tdir' does not exist\n";
    }
    my $runname = basename($tdir);
    my ($fqtdir, $asperadir) = MakeAsperaDir($tdir);

    #   Now see if we can get something like a Manifest
    $ENV{ASPERA_SCP_FILEPASS} = $opts{decryptpassword};
    my $cmd = "$opts{ascpshares} -hook -L $asperadir $opts{broadurl} $opts{user}:$opts{password} /$opts{rmtdir} $asperadir/Manifest";
    if ($opts{verbose}) { print "Executing: $cmd\n"; }
    system($cmd) &&
        die "$Script - Unable to fetch manifest file\n";
    system("ls -l $asperadir/Manifest");          # Show file we fetched 

    #   Here is data to be saved.  Assumes destination looks like .../center/dirname
    my %dbdata = (                      # Save key information
        sharedir => $opts{rmtdir},
        destdir => $fqtdir,
        dirname => $runname,
        center => $opts{center},
        user => $opts{user},
        password => $opts{password},
        decryptpassword => $opts{decryptpassword},
        url => $opts{broadurl},
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
    system("mv $dbdata{keyfile} $asperadir/$opts{keyfile}") &&
        die "$Script - Unable to save keyfile from '$opts{ascpshares} -hook': $!\n";
    $dbdata{keyfile} = "$asperadir/$opts{keyfile}";   # Change where keyfile is now
    if ($dbdata{basecmd} !~ /^(.+) -i \S+/) {
        die "$Script - unable to parse ascpshares basecmd: $dbdata{basecmd}\n";
    }
    $dbdata{basecmd} = "$1 -i $asperadir/$opts{keyfile} --file-crypt=decrypt";     # Correct command to fetch later
    unlink($opts{hooksfile});
    
    #   Now move the Manifest file to the proper place
    $_ = `ls $asperadir/Manifest/Manifest*`;
    chomp();
    if ("$_" eq "") { die "$Script - Unable to retrieve '$asperadir/Manifest/Manifest*'\n"; }

    my $manfile = "$fqtdir/Manifest.txt";
    system("mv -f $_ $manfile") &&
        die "$Script - Unable to move Manifest.txt to $fqtdir: $!\n";
    system("chmod 0444 $manfile") &&
        die "$Script - Unable to force Manifest.txt as readonly\n";
    rmdir("$asperadir/Manifest");
    print "Saved Manifest in $fqtdir\n";

    MakeAsperaData($manfile, $asperadir, $runname, \%dbdata);

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

    #   Saved details for this destinaltion are here
    my $asperadir = $opts{logdir} . '/' . GeteAsperaDest($opts{dest});

    #   Read database data for this dest
    my $f = "$asperadir/$opts{database}";
    open(IN,$f) ||
        die "$Script - Unable to read saved information for '$opts{dest}' from '$f': $!\n";
    my %dbdata = ();
    while (<IN>) {
        if (! /(\S+)=(.+)/) { die "$Script - Unable to parse line in '$f'\n  LINE=$_\n"; }
        $dbdata{$1} = $2;
    }
    close(IN);

    #   Prepare to build the ascp command
    my $realfile = basename($sendfile);
    mkdir("$dbdata{destdir}/.$realfile");
    mkdir("$asperadir/logs");                 # Aspera logs go here
    my $logdir = "$asperadir/logs/$realfile";
    mkdir($logdir);

    my ($rmtfile, $cmd);
    #   BROAD - Build the aspera command to fetch this file
    if ($opts{center} eq 'broad') {
        #   Get variables specific to this particular file
        $rmtfile = "$dbdata{sharedir}/$sendfile.aspera-env";
        $cmd = "$opts{asperabroad} $dbdata{user}:$dbdata{password} $rmtfile";
        open(IN, "$cmd |") ||
            die "$Script - Unable to get data for $opts{center} file: CMD=$cmd\n";
        my ($user, $host, $token, $srcfile) = split(' ', <IN>);
        close(IN);

        #   Build command to fetch file
        $ENV{ASPERA_SCP_FILEPASS} = $dbdata{decryptpassword};
        $cmd = "$dbdata{basecmd} -W $token -L $logdir $user\@$host:$rmtfile $dbdata{destdir}/.$realfile";
    }

    #   UW - Build the aspera command to fetch this file
    if ($opts{center} eq 'uw') {
        $rmtfile = "$dbdata{sharedir}/$sendfile";
        $ENV{ASPERA_SCP_PASS} = $dbdata{password};
        $cmd = "$opts{ascp} --src-base=$dbdata{rmtdir} -L $logdir " .
            "$dbdata{user}\@$opts{uwhost}:$rmtfile $dbdata{destdir}/.$realfile";        
    }

    if (! $cmd) { die "$Script - Unable to build cmd for center '$opts{center}'\n"; }

    #   Run ascp tp get the file. If there was an error, capture some
    #   important information and try again
    if ($opts{verbose}) { print "Executing: $cmd\n"; }
    my $rc = 0;
    foreach my $try (1..10) {
        $rc = system($cmd);
        if (! $rc) { last; }
        my $tmp = "/tmp/$$.topmedaspera";
        system("ls -la $logdir");
        #   Read log and extract key error message
        system("grep ERR $logdir/* > $tmp");
        my $errs = '';
        if (open(IN, $tmp)) {
            while(<IN>) {
                if (! /ERR Sess/) { next; }
                $errs = $_;
                print $errs;
            }
            close(IN);
        }
        unlink($tmp);
        system("echo \"Failed to fetch '$sendfile' try=$try\n$errs\nSee logs at $logdir\"|mail -s err tpg\@hps.com");
        print "Retrying try=$try\n";
        sleep(45*$try);
    }
    if ($rc) { die "$Script - Unable to fetch file '$sendfile\n"; }

    #   File has arrived, make it visible to the outside world
    system("ls -l $dbdata{destdir}/.$realfile");     # Show tmp file we fetched
    $cmd = "echo this should never be\n"; 
    if ($opts{center} eq 'broad') {
        $cmd = "mv $dbdata{destdir}/.$realfile/$realfile $dbdata{destdir}";
    }
    if ($opts{center} eq 'uw') {
        $cmd = "mv $dbdata{destdir}/.$realfile $dbdata{destdir}/$realfile";
    }
    system($cmd) &&
        die "$Script - Unable to rename newly arrived file:  CMD=$cmd\n";
    chmod(0444, "$dbdata{destdir}/$realfile");
    if ($opts{center} eq 'broad') { rmdir("$dbdata{destdir}/.$realfile"); }
    system("ls -l $dbdata{destdir}/$realfile");     # Show file we fetched 
}

#==================================================================
# Subroutine:
#   OLDFetchFile($sendfile)
#
#   Fetch a file (maybe as dir/file) as a temporary.
#   When successful, rename it to what it should be.
#==================================================================
sub OLDFetchFile {
    my ($sendfile) = @_;               

    #   Squash the target directory so we have a unique identifier
    my $asperadir = $opts{logdir} . '/' . GeteAsperaDest($opts{dest});    # Saved fetch details here
    #   Read database data for this dest
    my $f = "$asperadir/$opts{database}";
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
    my $cmd = "$opts{asperabroad} $dbdata{user}:$opts{password} $rmtfile";
    open(IN, "$cmd |") ||
        die "$Script - Unable to get data for $opts{center} file: CMD=$cmd\n";
    my ($user, $host, $token, $srcfile) = split(' ', <IN>);
    close(IN);

    #   Build command to fetch file
    my $realfile = basename($sendfile);
    chomp($realfile);
    $ENV{ASPERA_SCP_FILEPASS} = $dbdata{decryptpassword};
    mkdir "$dbdata{destdir}/.$realfile";
    mkdir("$asperadir/logs");                 # Aspera logs go here
    my $logdir = "$asperadir/logs/$realfile";
    mkdir($logdir);
    $cmd = "$dbdata{basecmd} -W $token -L $logdir $user\@$host:$rmtfile $dbdata{destdir}/.$realfile";
    #   Run ascp tp get the file. If there was an error, capture some
    #   important information and try again
    if ($opts{verbose}) { print "Executing: $cmd\n"; }
    my $rc = 0;
    foreach my $try (1..10) {
        $rc = system($cmd);
        if (! $rc) { last; }
        my $tmp = "/tmp/$$.topmedaspera";
        system("ls -la $logdir");
        #   Read log and extract key error message
        system("grep ERR $logdir/* > $tmp");
        my $errs = '';
        if (open(IN, $tmp)) {
            while(<IN>) {
                if (! /ERR Sess/) { next; }
                $errs = $_;
                print $errs;
            }
            close(IN);
        }
        unlink($tmp);
        system("echo \"Failed to fetch '$sendfile' try=$try\n$errs\nSee logs at $logdir\"|mail -s err tpg\@hps.com");
        print "Retrying try=$try\n";
        sleep(45*$try);
    }
    if ($rc) { die "$Script - Unable to fetch file '$sendfile\n"; }

    system("ls -l $dbdata{destdir}/.$realfile");     # Show tmp file we fetched 
    $cmd = "mv $dbdata{destdir}/.$realfile/$realfile $dbdata{destdir}";
    system($cmd) &&
        die "$Script - Unable to rename newly arrived file:  CMD=$cmd\n";
    chmod(0444, "$dbdata{destdir}/$realfile");

    rmdir("$dbdata{destdir}/.$realfile");
    system("ls -l $dbdata{destdir}/$realfile");     # Show file we fetched 
}


#==================================================================
# Subroutine:
#   (fqtdir, asperadir) = MakeAsperaDir($targetdir)
#
#   This has the side effect of saving key information for subsequent
#   data being copied to the target.
#
#   Returns:
#       fully qualified path to tdir, no symlinks
#       asperadir - directory where we save all the aspera information
#==================================================================
sub MakeAsperaDir {
    my ($tdir) = @_;

    my $fqtdir = abs_path($tdir);           # All our data goes here
    system("mkdir -p $fqtdir") && die "$Script - unable to create directory '$fqtdir'\n";
    system("mkdir -p $opts{logdir}") && die "$Script - unable to create directory '$opts{logdir}'\n";

    #   Squash the target directory so we have a unique identifier
    my $asperadir = $opts{logdir} . '/' . GeteAsperaDest($fqtdir);  # Save fetch details here
    system("mkdir -p $asperadir") && die "$Script - unable to create directory '$asperadir\n";
    print "Details for '$tdir' will be saved in '$asperadir'\n";
    return ($fqtdir, $asperadir);
}


#==================================================================
# Subroutine:
#   ($destkey) = GeteAsperaDest($dir)
#
#   Calculate the destination key for the aspera data we save
#
#   Returns:
#       string
#==================================================================
sub GeteAsperaDest {
    my ($dir) = @_;
   return basename($dir);
}

#==================================================================
# Subroutine:
#   MakeAsperaData($manfile, $asperadir, $runname, $dbdataref)
#
#   Create the 'database' of information used by fetch
#
#==================================================================
sub MakeAsperaData {
    my ($manfile, $asperadir, $runname, $dbdataref) = @_;

    #   Create our master key file so we can fetch more files
    open(OUT, '>' . "$asperadir/$opts{database}") ||
        die "$Script - Unable to create file '$asperadir/$opts{database}': $!\n";
    foreach my $k (keys %$dbdataref) { print OUT "$k=$dbdataref->{$k}\n"; }
    close(OUT);
    print "Saved information for '$runname' in $asperadir\n";

    #   Read the manifest and generate the commands to fetch the file
    my $fetchcmds = "$asperadir/fetch-commands.txt";
    open(OUT, '>' . $fetchcmds) ||
        die "$Script - Unable to create file '$fetchcmds': $!\n";
    open(IN, $manfile) ||
        die "$Script - Unable to read file '$manfile': $!\n";
    my $n = 0;
    while (<IN>) {
        if (! /^(\S+)/) { next; }
        my $f = $1;
        $n++;
        print OUT "$0 -dest $runname fetch $f\n";
    }
    close(IN);
    close(OUT);
    print "Created list of fetch commands in $fetchcmds\n";
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__
topmedaspera.pl [options] -center broad \
  -user SN0096478 -password FA3XCLEI7JZAMX8 -decryptpassword yk2iCI7zefH3nmu \
  -rmtdir SN0096478 \
  manifest targetdir
  or
topmedaspera.pl [options] -center uw \
  -user nhlbi-aspera -password gobuN4je \
  -rmtdir /nhlbi_macrogen_2/arnett_nhlbi_wgs_hypergen_utah/batch_9 \
  manifest targetdir
  or
topmedaspera.pl [options] -dest 2016jun23 fetch file


=head1 NAME

topmedaspera.pl - Automate download of data using ASPERA

=head1 SYNOPSIS

  topmedaspera.pl -center uw -user RMTUSER -pass RMTPW \
    --rmtdir /DIRECTORY/AT/REMOTE/SITE manifest  /incoming/topmed/uw/2016mar33

  topmedaspera.pl --dest 2016mar33 fetch 145474.final.sqz.bam
  
  
  topmedaspera.pl -center broad -user SN0096478 -pass RMTPW \
    decryptpassword yk2iCI7zefH3nmu /
    --rmtdir /SN0096478 manifest  /incoming/topmed/broad/2016jun34

  topmedaspera.pl --dest 2016jun34 fetch a/NWD123456

=head1 DESCRIPTION

This program can automate the downloading of data from another center.
You must specify the center, since each has it's own peculiarities of
how the data looks.

The process consists of two steps. The first is to fetch the manifest
of files to be downloaded. This ultimately results in a Manifest.txt
being created in the target directory. The monitor will notice this and 
then create database entries for all data.

The second step consists of multiple invocations of this script which will
download a single file. Downloading one file at a time (probably run concurrently)
allows the monitor to notice when the file arrives and to begin processing 
the files as they come in.


=head1 OPTIONS

=over 4

=item B<-center NAME>

Specifies a specific center name from which we will download the data, e.g. B<uw>.

=item B<-decryptpassword string>
This must be provided for a BROAD manifest call.

Specifies a descryption string used by the BROAD.

=item B<-dest target>

Specifies a target set of information created when a manifest is downloaded.
This is the target directory where the data is downloaded (e.g. basename targetdirectory).
This is only used by B<fetch>

=item B<-help>

Generates this output.

=item B<-password string>

Specifies the password to be used when accessing the remote site.
This must be provided for a manifest call.

=item B<-rmtdir NAME>

Specifies the local directory where the data will be downloaded.
This must be provided for a manifest call.

=item B<-user string>

Specifies the login account to use when downloading the data.
This must be provided for a manifest call.

=item B<-verbose>

Provided for developers to see additional information.

=back

=head1 PARAMETERS

There are two parameters to this program. 

=over 4

=item B<manifest targetdirectory>

Specifies we are starting a download and need to fetch a list of all the files
to be downloaded - e.g. a Manifest.

This requires you specify the local destination directory where the files 
will be downloaded to. This directory must already exist.

=item B<fetch  file>

Specifies a specific file that will be downloaded. This is in the form of the 
remote site (e.g. B<a/NWD123456>.
The file will be downloaded as 'dot' file, so it will not normally be visible.
Once the file is successfully downloaded, it will be renamed to the basename
of the original file that was fetched (e.g. B<NWD123456>).

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

