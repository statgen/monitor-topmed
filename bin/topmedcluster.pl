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
    squeuecmd => "/usr/cluster/bin/squeue -a -o '%.9i %.15P %.18q %.14j %.8u %.2t %.8M %.6D %R %.9n' -p topmed-working",    
    verbose => 0,
    maxlongjobs => 6,
    consoledir => 'working/topmed-output',
);

Getopt::Long::GetOptions( \%opts,qw(
    help verbose maxlongjobs=i
    )) || die "$Script - Failed to parse options\n";

#   Simple help if requested
if ($#ARGV < 0 || $opts{help}) {
    my $m = "$Script [options]";
    warn "$m squeue\n" .
        "More details available by entering: perldoc $0\n\n";
    if ($opts{help}) { system("perldoc $0"); }
    exit 1;
}
my $fcn = shift @ARGV;

#--------------------------------------------------------------
#   Execute the command provided
#--------------------------------------------------------------
if ($fcn eq 'squeue')    { SQueue(@ARGV); exit; }

die "$Script  - Invalid function '$fcn'\n";
exit;

#==================================================================
# Subroutine:
#   SQueue()
#
#   Print summary of topmed-related jobs
#==================================================================
sub SQueue {
    #my ($bamid, $set) = @_;
    my $cmd = $opts{squeuecmd} . '| grep topmed';
    my %partitions = ();            # Unique partition names
    my %qos = ();                   # Unique QOS names
    my %qosheld = ();               # Count of held jobs in QOS
    my %running = ();
    my %mosixrunning = ();          # Counts of jobs by nomosix users
    my %queued = ();
    my %nottopmed = ();             # Not my user list

    DBConnect($opts{realm});
    foreach my $l (split("\n", `$cmd`)) {
        my @c = split(' ', $l);
        if ($c[1] eq 'nomosix' && $c[8] =~ /topmed/) {    # nomosix on topmed node
            $mosixrunning{$c[4]}++;    # Count of running jobs for this user
            next;
        }
        if ($c[1] !~ /^topmed/) { next; }    # Only want my partition of interest
        if ($c[4] ne 'topmed' && $c[4] ne 'schelcj') {  # Save partition and not topmed user
            $nottopmed{$c[1]}{$c[4]} = 1;
            next;
        }
        $partitions{$c[1]} = 1;
        if ($c[5] eq 'PD') {        # Queued
            push @{$queued{$c[1]}{data}},$l;
            $queued{$c[1]}{count}++;
            $qos{$c[2]}{queued}++;
            if ($l =~ /held state/ || $l =~ /JobHeldUser/) { $qos{$c[2]}{held}++; }
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
        my $s = join(' ',sort keys %{$nottopmed{$p}});
        my $k = scalar(keys %{$nottopmed{$p}});
        printf("  %-18s %3d running / %3d queued / %3d foreign user: $s \n",
            $p, $running{$p}{count}, $queued{$p}{count}, $k);
    }
    print "\n";

    #   Show summary of QOS
    print "QOS Summary\n";
    foreach my $q (sort keys %qos) {
        my $s = $qos{$q}{held} || '';
        if ($s) { $s .= ' jobs held'; }
        if (! defined($qos{$q}{running}))  { $qos{$q}{running} = 0; }
        if (! defined($qos{$q}{queued}))   { $qos{$q}{queued} = 0; }
        my $qq = $q;                    # Patch normal to be something meaningful?
        if ($q eq 'normal') { $qq = 'default_qos'; }
        printf("  %-18s %3d running / %3d queued   %s\n",
            $qq, $qos{$q}{running}, $qos{$q}{queued}, $s); 
    }
    print "\n";

    #   Show summary of running on partitions
    print "Running jobs per host\n";
    my $format = '    %-11s  %s';
    foreach my $p (sort keys %partitions) {
        my %hosts = ();
        foreach my $l (@{$running{$p}{data}}) {
            my @c = split(' ', $l);
            $c[3] =~ s/\d+\-//;             # Remove bamid from jobname
            $c[3] =~ s/^NWD\d+/NWD/;        # This is remapping job
            $hosts{$c[8]}{$c[3]}++;         # Increment $hosts{topmed2}{cram}
        }
        print "   Partition $p:\n";
        my $format = '    %-11s  %s';
        printf ($format . "\n", 'Host', 'Job types and count');
        foreach my $h (sort keys %hosts) {
            my $s = '';
            foreach my $jobname (sort keys %{$hosts{$h}}) {
                $s .= sprintf('  %-12s', $jobname . "  [" . $hosts{$h}{$jobname} . '] ');
            }
            printf ($format . "\n", $h, $s);
        }
    }
    print "\n";
    if (%mosixrunning) {
        $format = "    %-8s  %s\n";
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

  topmedcluster.pl squeue       # Show summary of queues


=head1 DESCRIPTION

This program supports simple commands to show information about the cluster systems.


=head1 OPTIONS

=over 4

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

B<squeue>
Use this to provide a smart summary of SLURM queues.


=head1 EXIT

If no fatal errors are detected, the program exits with a
return code of 0. Any error will set a non-zero return code.

=head1 AUTHOR

Written by Terry Gliedt I<E<lt>tpg@umich.eduE<gt>> in 2017-2016 and is
is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; See http://www.gnu.org/copyleft/gpl.html

=cut

