#!/usr/bin/env perl
###################################################################
#
# Name: topmedcluster.pl
#
# Description:
#   Use this program to query the cluster queues in a smart way.
#   This program is tightly coupled with /monitor/topmed/index.php
#
#   This was part of topmedcmd.pl, but had to be split out because
#   the web server did not have enough Perl modules installed that
#   we really did not need for this functionality.
#
# ChangeLog:
#   $Log: topmedcluster.pl,v $
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
use My_DB;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
our %opts = (
    realm => '/usr/cluster/topmed/etc/.db_connections/topmed',
    netdir => '/net/topmed',
    bamfiles_table => 'bamfiles',
    runs_table => 'runs',
    squeuecmd => "/usr/cluster/bin/squeue -a -o '%.9i %.15P %.18q %.18j %.8u %.2t %.8M %.6D %R %.9n' -p topmed-working",
    maxlongjobs => 6,
    conf => "/usr/cluster/topmed/etc/topmedthrottle.conf",
    consoledir => 'working/topmed-output',
    squeuefile => '/run/shm/squeue.results',
    verbose => 0,
    squeuefiletime => 30,               # Refresh file if older than this seconds
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose maxlongjobs=i force
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options]";
    warn "$m squeue|summary|newcache\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;
my @squeuelines = ();                   # Global data for functions

#--------------------------------------------------------------
#   Get squeue results, save in memory file every so often
#--------------------------------------------------------------
if ($fcn eq 'newcache')  { CacheFile('new'); exit; }

if ($opts{force}) { CacheFile('new'); }         # Create new cache
else { CacheFile('read'); }                     # Just cache
if ($fcn eq 'squeue')    { SQueue(@ARGV); exit; }
if ($fcn eq 'summary')   { Summary(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   CacheFile(flag)
#
#   Read cachefile into @squeuelines. If 'new' provided, create it new
#==================================================================
sub CacheFile {
    my $flag = $_[0] || 'read';
    my $file = $opts{squeuefile};

    my $in;
    my $filetime;
    if ($flag eq 'read' && -f $file) {
        $filetime = time() - (stat($file))[9];       # How old is file
        if ($filetime > $opts{squeuefiletime}) { $flag = 'new'; }   # Old, force create
    }
    else { $flag = 'new'; }                 # No file, force create of file

    #   If recent file, read and return
    if ($flag ne 'new') {
        open($in, $file) ||
            die "$Script - Unable to open file '$file': $!\n";
        while (<$in>) {
            if (/^#/) { next; }
            push @squeuelines,$_;
        }
        close($in);
        return;
    }

    #   Create a new file, start with squeue data
    $file = "$opts{squeuecmd} |";
    open($in, $file) ||
        die "$Script - Unable to open file '$file': $!\n";
    while (<$in>) {
        if (/^#/) { next; }
        push @squeuelines,$_;
    }   
    close($in);

    #   Now read the conf file so we can get intended max values
    my @maxes = ();
    open($in, $opts{conf}) ||
        die "$Script - Unable to open file '$opts{conf}': $!\n";
    my %hostmaxperaction = ();
    my %systemmaxperaction = ();
    my $OPENAMP = '{';
    my $CLOSEAMP = '}';
    while (<$in>) {
        if (/^#/) { next; }
        if (! /\S/) { next; }
        if (/^job\s+(\S+)\s+$OPENAMP/i) {
            my $action = $1;
            while (<$in>) {
                if (/\s*$CLOSEAMP/) { last; }
                if (/\s*max\s*=\s*(\S+)/) { $systemmaxperaction{$action} = $1; }
                if (/\s*maxperhost\s*=\s*(\S+)/) { $hostmaxperaction{$action} = $1; }
            }
        }
    }
    close($in);

    #   Make summary of maxes per action
    foreach my $action (keys %hostmaxperaction) {
        push @maxes,"action= $action sysmax= $systemmaxperaction{$action} hostmax= $hostmaxperaction{$action}\n";
    }

    #   Write out summary of all data of interest
    my $out;
    open($out, '>'. $opts{squeuefile}) ||
        die "$Script - Unable to create file '$opts{squeuefile}': $!\n";
    foreach (@maxes, @squeuelines) { print $out $_; }
    close($out);
} 

#==================================================================
# Subroutine:
#   Summary()
#
#   Print summary of jobs
#==================================================================
sub Summary {
    #my ($bamid, $set) = @_;
    my %jobtype = ();               # Has of jobtypes queued or running
    foreach my $l (@squeuelines) {
        my @c = split(' ', $l);
        #   Should be none of these, but if so, we want to know about it
        if ($c[4] ne 'topmed') { next; }    # Only interested in my jobs
        if ($c[5] eq 'PD') {        # Queued
            if ($c[3] =~ /\d+-([a-z]+)/) {
                my $act = $1;
                $jobtype{$act}{queued}++;
                $jobtype{$act}{qhosts}{$c[9]}++;
            }
            next;
        }
        if ($c[5] eq 'R') {         # Running
            if ($c[3] =~ /\d+-([a-z]+)/) {
                my $act = $1;
                $jobtype{$act}{running}++;
                $jobtype{$act}{rhosts}{$c[9]}++;
            }
            next;
        }
    }
    foreach my $jtype (sort keys %jobtype) {
        if (! exists($jobtype{$jtype}{queued}))  { $jobtype{$jtype}{queued} = 0; }
        if (! exists($jobtype{$jtype}{running})) { $jobtype{$jtype}{running} = 0; }
        my $r = '';
        foreach my $h (sort keys %{$jobtype{$jtype}{rhosts}}) {
            $r .= "$h $jobtype{$jtype}{rhosts}{$h} ";
        }
        my $q = '';
        foreach my $h (sort keys %{$jobtype{$jtype}{qhosts}}) {
            $q .= "$h $jobtype{$jtype}{qhosts}{$h} ";
        }
        print "$jtype queued=$jobtype{$jtype}{queued} running=$jobtype{$jtype}{running} qhosts=$q rhosts=$r\n";
    }
}

#==================================================================
# Subroutine:
#   SQueue()
#
#   Print summary of topmed-related jobs
#==================================================================
sub SQueue {
    #my ($bamid, $set) = @_;
    my $cmd = $opts{squeuecmd};
    my %partitions = ();            # Unique partition names
    my %qos = ();                   # Unique QOS names
    my %qosheld = ();               # Count of held jobs in QOS
    my %running = ();
    my %mosixrunning = ();          # Counts of jobs by nomosix users
    my %queued = ();
    my %nottopmed = ();             # Not my user list
    my %jobtype = ();               # If one QOS, we want to know types of jobs queued

    DBConnect($opts{realm});
    foreach my $l (@squeuelines) {
        my @c = split(' ', $l);
        #   Should be none of these, but if so, we want to know about it
        if ($c[1] eq 'nomosix' && $c[8] =~ /topmed/) {    # nomosix on topmed node
            $mosixrunning{$c[4]}++;    # Count of running jobs for this user
            next;
        }
        if ($c[4] ne 'topmed') {            # Save partition and not topmed user
            $nottopmed{$c[1]}{$c[4]}++;
            next;
        }
        $partitions{$c[1]} = 1;
        if ($c[5] eq 'PD') {        # Queued
            push @{$queued{$c[1]}{data}},$l;
            $queued{$c[1]}{count}++;
            $qos{$c[2]}{queued}++;
            if ($l =~ /held/i) { $queued{$c[1]}{held}++; }     # Job held somehow
            if ($c[3] =~ /\d+-([a-z]+)/) { $jobtype{$1}++; }   # Count of types of jobs
            next;
        }
        if ($c[5] eq 'R') {         # Running
            push @{$running{$c[1]}{data}},$l;
            $running{$c[1]}{count}++;
            $qos{$c[2]}{running}++;
            next;
        }
    }
    
    #   Show summary of partitions
    print "Partition Summary\n";
    foreach my $p (sort keys %partitions) {
        if (! defined($queued{$p}{count}))  { $queued{$p}{count} = 0; }
        if (! defined($running{$p}{count})) { $running{$p}{count} = 0; }
        if (! defined($queued{$p}{held}))   { $queued{$p}{held} = 0; }
        my $s = '';
        foreach my $u (sort keys %{$nottopmed{$p}}) {
            $s .= "$u [${$nottopmed{$p}}{$u}] ";
        }
        my $k = scalar(keys %{$nottopmed{$p}}) || 0;
        printf("  %-18s %3d running / %3d queued / %3d held / %3d foreign user: $s \n",
            $p, $running{$p}{count}, $queued{$p}{count}, $queued{$p}{held}, $k);
        if (%jobtype) {
            print '    queued: ';
            foreach my $jtype ( sort keys %jobtype) {
                print $jtype . ' [' . $jobtype{$jtype} . ']  ';
            }
            print "\n";
        }
    }
    print "\n";

    #   Show summary of QOS   Only one QOS, no need for a summary
    #print "QOS Summary\n";
    #foreach my $q (sort keys %qos) {
    #    my $s = $qos{$q}{held} || '';
    #    if ($s) { $s .= ' jobs held'; }
    #    if (! defined($qos{$q}{running}))  { $qos{$q}{running} = 0; }
    #    if (! defined($qos{$q}{queued}))   { $qos{$q}{queued} = 0; }
    #    my $qq = $q;                    # Patch normal to be something meaningful?
    #    if ($q eq 'normal') { $qq = 'default_qos'; }
    #    printf("  %-18s %3d running / %3d queued   %s\n",
    #        $qq, $qos{$q}{running}, $qos{$q}{queued}, $s); 
    #}
    #print "\n";

    #   Show summary of running on partitions
    print "Running jobs per host\n";
    foreach my $p (sort keys %partitions) {
        my $format = '    %-9s  %s';
        my %hosts = ();
        my %total = ();                     # Get counts of number of each jobname
        foreach my $l (@{$running{$p}{data}}) {
            my @c = split(' ', $l);
            if ($c[3] =~ /\d+-([a-z]+)/) { $c[3] = $1; }
            #$c[3] =~ s/^NWD\d+/NWD/;       # This is remapping job
            $hosts{$c[8]}{$c[3]}++;         # Increment $hosts{topmed2}{cram}
        }
        print "   Partition $p:\n";
        printf ($format . "\n", 'Host', 'Job types and count');
        foreach my $h (sort keys %hosts) {
            my $s = '';
            foreach my $jobname (sort keys %{$hosts{$h}}) {
                $s .= sprintf('  %-10s', $jobname . "  [" . $hosts{$h}{$jobname} . '] ');
                $total{$jobname} += $hosts{$h}{$jobname};
            }
            printf ($format . "\n", $h, $s);
        }
        #   Show total number of jobs per action
        my $s = '';
        foreach my $jobname (sort keys %total) {
            $s .= sprintf('  %-10s', $jobname . " [" . $total{$jobname} . '] ');
        }
        printf ($format . "\n", 'Totals:', $s);
    }
    print "\n";
    if (%mosixrunning) {
        my $format = "    %-8s  %s\n";
        print "  Partition nomosix jobs on topmed hosts:\n";
        printf ($format, 'User', 'Job count');
        foreach my $u (sort keys %mosixrunning) {
            printf($format, $u, "  [" . $mosixrunning{$u} . '] ');
        }
        print "\n";
    }
    
    #   Show longest running jobs
    print "Longest running jobs per partition\n";
    foreach my $p (sort keys %partitions) {
        if (! defined($running{$p}{data})) { next; }
        print "  $p:\n";
        my %longjobs = ();
        foreach my $l (@{$running{$p}{data}}) {
            my @c = split(' ', $l);
            if ($c[3] =~ /^NWD/) { next; }  # Ignore jobs I did not submit
            if ($c[3] =~ /^nhlbi/) { next; }
            my $t = $c[6];                  # Normalize all times to dd-hh:mm:ss
            my $days = '00-';
            if ($t =~ /^(\d+\-)(.+)/) { $days = $1; $t = $2; }
            if ($t =~ /^(\d+):(\d+):(\d+)$/) { $t = sprintf('%02d:%02d:%02d', $1, $2, $3); }
            elsif ($t =~ /^(\d+):(\d+)$/) { $t = sprintf('00:%02d:%02d', $1, $2); }
            elsif ($t =~ /^:(\d+)$/) { $t = sprintf('00:00:%02d', $1); }
            $longjobs{$t} .= '  ' . $l . "\n";      # Descriptions of longest running job
        }
        my $i = $opts{maxlongjobs};
        foreach my $t (reverse sort keys %longjobs) {
            foreach (split("\n", $longjobs{$t})) {
                if (! $_) { next; }
                print ReFormatPartitionData($_);
                $i--;
                if ($i <= 0) { last; }
            }
            if ($i <= 0) { last; }
        }
    }
    print "\n";

    #   Get summary of IO data for key nodes
    my @iostatinfo = ();
    if (opendir(my $dh, "$opts{netdir}/$opts{consoledir}")) {
        @iostatinfo = sort sorthost grep { /iostat.topmed\d*.txt/ } readdir($dh);
        closedir $dh;
        my $in;
        my $prefix = ' ' x 12;
        my $date = $prefix;
        my $hdr = 0;
        my $iohdr = '   %user %iowait  %idle';
        my $iodata = '';
        print "    I/O activity for /incoming and /working filesystems per host\n";
        foreach my $f (@iostatinfo) {
            my $host = $prefix;
            if (open($in,"$opts{netdir}/$opts{consoledir}/$f")) {
                if ($f =~ /(topmed\d*)/) { $host = sprintf('%-12s', $1); }
                #   Linux 4.2.0-42-generic (topmed) 	02/10/2017 	_x86_64_	(120 CPU)
                #
                #   avg-cpu:  %user   %nice %system %iowait  %steal   %idle
                #              0.35    0.00    1.66    0.55    0.00   97.44
                #
                #   Device:            tps    MB_read/s    MB_wrtn/s    MB_read    MB_wrtn
                #   md20              0.40         0.00         0.01          0          0
                #   md21           7782.60       213.60        13.17       1068         65
                while (<$in>) {
                    if (/Linux.+(\d\d\/\d\d\/20\d\d)/) { $date = sprintf('%-12s', $1); next; }
                    if (/Device/ && (! $hdr)) {
                        print $date . ReFormatDeviceData($_) . $iohdr . "\n";
                        $hdr++;
                        next;
                    }
                    if (/avg-cpu/) {
                        $_ = <$in>;             # Get actual numbers
                        $iodata = ReFormatIOData($_);
                        next;
                    }
                    if (/md2\d/) {
                        print $host . ReFormatDeviceData($_) . $iodata ."\n";
                        $host = $prefix;
                        $iodata = '';
                        next;
                    }
                }
                close($in);
            }
        }
        print "\n";
    }
    return;
}

#==================================================================
# Subroutine:
#   sorthost    Used to sort topmed hosts in numnerical order
#
#   Return reformatted device data
#==================================================================
sub sorthost {
    my $cmp1 = 0;
    my $cmp2 = 0;
    if ($a =~ /topmed(\d+)/) { $cmp1 = sprintf('%02d',$1); }
    if ($b =~ /topmed(\d+)/) { $cmp2 = sprintf('%02d',$1); }
    return ($cmp1 <=> $cmp2);
}

#==================================================================
# Subroutine:
#   ReFormatIOData($str)
#
#   Return reformatted device data
#==================================================================
sub ReFormatIOData {
    my ($str) = @_;
    if ((! defined($str)) || $str eq '') { return ''; }
    my @c = split(' ', $str);
    if ($#c < 5) { return ''; }
    return sprintf('   %5s %7s %6s', $c[0], $c[3], $c[5]);
}

#==================================================================
# Subroutine:
#   ReFormatDeviceData($str)
#
#   Return reformatted device data
#==================================================================
sub ReFormatDeviceData {
    my ($str) = @_;
    if ((! defined($str)) || $str eq '') { return ''; }
    my @c = split(' ', $str);
    if ($c[0] eq 'md20') { $c[0] = '/incoming'; }
    else {
        if ($c[0] eq 'md21') { $c[0] = '/working'; }
    }
    if ($#c < 3) { return ''; }
    return sprintf('%-9s %10s %10s', $c[0], $c[2], $c[3]);
}

#==================================================================
# Subroutine:
#   ReFormatPartitionData($str)
#
#   Return reformatted partition data
#==================================================================
sub ReFormatPartitionData {
    my ($str) = @_;
    my $fmt = '    %-9s %-9s %-12s %-10s %-9s %s  %s';
    my @c = split(' ', $str);
    my $nwdid = '?';
    my $dir = '?';
    if ($c[3] =~ /^(\d+)/) {
        my $bamid=$1;
        my $sth = DoSQL("SELECT runid,expt_sampleid FROM $opts{bamfiles_table} WHERE bamid=$bamid");
        if ($sth) {
            my $href = $sth->fetchrow_hashref;
            if (exists($href->{expt_sampleid})) { $nwdid = $href->{expt_sampleid}; }
            if (exists($href->{runid})) {
                $sth = DoSQL("SELECT dirname FROM $opts{runs_table} WHERE runid=$href->{runid}");
                if ($sth) {
                    $href = $sth->fetchrow_hashref;
                    $dir = $href->{dirname};
                }
            }
        }
    }
    return sprintf($fmt, $c[0], $nwdid, $c[2], $c[3], $c[8], $c[6], $dir)  . "\n";
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmedcluster.pl - Query for cluster system information


=head1 SYNOPSIS

  topmedcluster.pl squeue           # Show summary of queues
  topmedcluster.pl -force squeue    # Create new cache and show summary of queues
  topmedcluster.pl summary          # Show summary of jobtypes
  topmedcluster.pl newcache         # Create new cache file, nothing else


=head1 DESCRIPTION

This program supports simple commands to show information about the cluster systems.


=head1 OPTIONS

=over 4

=item B<-force>

Forces the creation of a new cache file.

=item B<-help>

Generates this output.

=item B<-maxlongjobs N>

Specifies the number of long running jobs to show. This defaults to B<6>.

=item B<-verbose>

Provided for developers to see additional information.

=back


=head1 PARAMETERS

Parameters to this program are grouped into several groups which are used
to deal with specific sets of information in the monitor databases.

B<summary>
Use this to provide a short summary. No longer used I think.

B<squeue>
Use this to provide a smart summary of SLURM queues.

B<newcache>
Use this to force creation of a new cache file. No output.


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2016-2017 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

