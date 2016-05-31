#!/usr/bin/perl -I/usr/cluster/lib/perl5/site_perl
###################################################################
#
# Name: topmed_daemon.pl
#
# Description:
#   Use this program to issue commans on the topmed server to
#   get information for the web server which cannot see files
#   on the topmed machines.
#
#   Note: error codes and a decent example 'are shown at
#       http://search.cpan.org/~ether/libwww-perl-6.13/lib/LWP/Simple.pm
#       http://stackoverflow.com/questions/14425220/perl-http-server
#
# ChangeLog:
#   $Log: topmed_daemon.pl,v $
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
###################################################################
use strict;
use warnings;
use FindBin qw($Script);
use Getopt::Long;

use HTTP::Daemon;
use HTTP::Response;
use HTTP::Status;
use URI::Escape;

#--------------------------------------------------------------
#   Initialization - Sort out the options and parameters
#--------------------------------------------------------------
my %opts = (
    verbose => 0,
    restart => 0,
    restart_cmd => '/home/topmed/restart_daemon.sh -submit',
    stop => 0,
    port => 48109,
    fromips => '127.0.0.1 192.168.1.56 192.168.2.200 192.168.1.8',
    squeue_unused => '/usr/cluster/bin/squeue --format "%.9i %.10P %.12j %.8u %.2t %.10M %R"',
    df => '/bin/df -h /net/topmed/incoming /net/topmed/working /net/topmed2/incoming /net/topmed2/working  /net/topmed3/incoming /net/topmed3/working /net/topmed4/working /net/topmed4/incoming /net/topmed5/incoming /net/topmed5/working   /net/topmed6/incoming /net/topmed6/working',
    squeue => '/usr/cluster/monitor/bin/topmedcmd.pl squeue', 
    errorcheck => '/usr/cluster/monitor/bin/topmed_monitor.pl check',
    shownode => '/usr/cluster/bin/scontrol show node',
);
$opts{mycommand} = "$0 "  . join(' ', @ARGV);    # To support restart

Getopt::Long::GetOptions( \%opts,qw(
    help verbose port=n fromips=s
    )) || die "Failed to parse options\n";

#   Simple help if requested
if ($opts{help}) {
    warn "$Script [options]\n" .
        "This is a daemon to run SLURM squeue commands.\n";
    exit 1;
}
my %okips = ();
foreach (split(' ', $opts{fromips})) {
    $okips{$_} = 1;
}

#--------------------------------------------------------------
#   Launch daemon, get requests, process them
#--------------------------------------------------------------
my $daemon = HTTP::Daemon->new(LocalPort => $opts{port}) ||
    die "Unable to start daemon: $!\n";
if ($opts{verbose}) {
    print "Contact this server at URL: " . $daemon->url . "\n" .
    "Accepting connections from: $opts{fromips}\n";
}

#   Wait for something to show up, process it
while (! $opts{stop}) {
    my $conn = $daemon->accept;
    while (my $req = $conn->get_request) {
        my $ip = $conn->peerhost();
        #   We only accept queries from the right places
        if (! defined($ip)) { next; }
        if (! $okips{$ip}) {
            if ($opts{verbose}) { print "Cannot accept query from '$ip'  URI=" . $req->uri . "\n"; }
            next;
        }
        #   If this is a GET, parse out the pieces
        if ($req->method ne 'GET') {
            $conn->send_error(RC_FORBIDDEN);
            next;
        }
        my $uri = $req->uri;
        if (! defined($uri)) {
            $conn->send_error(RC_FORBIDDEN);
            if ($opts{verbose}) { print "URI path was undefined, how'd you do that?\n"; }
            next;
        }
        if ($opts{verbose}) { print "$ip =>'$uri'\n"; }

        my @uriparts = $uri =~ m|(?:([^:/?#]+):)?(?://([^/?#]*))?([^?#]*)(?:\?([^#]*))?(?:#(.*))?|;
        ProcessGet($conn, $uriparts[2], $uriparts[3]);
        if ($opts{stop}) { last; }
    }
    $conn->close;
    undef($conn);
}
if (! $opts{restart}) { exit; }

#   Restart turns out to be far more complex than I expected
#   I can restart with 'at' so for debugging, that's sufficient
#   This does not seem to behave. Maybe 'at' requires a terminal
system($opts{restart_cmd}) && print "Unable to restart\n";
exit;


#==================================================================
# Subroutine:
#   ProcessGet - Process GET request
#
# Arguments:
#   conn - connection object
#   path - URI path, eg /squeue
#   query - URI path after /, eg ?-p topmed
#==================================================================
sub ProcessGet {
    my ($conn, $path, $query) = @_;
    if (! defined($path)) { $path = ''; }
    if (! defined($query)) { $query = ''; }
    $query = uri_unescape($query);

    #   Handle squeue/queue
    if ($path =~ /squeue.(\S+)/) {
        my $part = Clean($1);
        my $cmd = $opts{squeue};
        if ($opts{verbose}) { print "  CMD=$cmd\n"; }
        my $s = `$cmd  2>&1`;
        SendText($conn, $s);
        return;
    }

    #   Handle shownode/node
    if ($path =~ /shownode.(\S+)/) {
        my $node = Clean($1);
        my $cmd = "$opts{shownode} $node";
        if ($opts{verbose}) { print "  CMD=$cmd\n"; }
        my $s = `$cmd | grep State=`;
        my @c = split(' ',$s);
        if ($opts{verbose}) { print "  RESULTS=$c[0]\n"; }
        SendText($conn, "$node $c[0]");
        return;
    }

    #   Handle stop and restart
    if ($path =~ /(restart|stop)\/(\S+)/) {
        my $part = $2;
        my $fcn = $1;
        my $s = '';
        if ($part ne 'please') {
            $s = "Say please and thank you\n";
        }
        else {
            if ($fcn eq 'restart') {
                $s .= "$0 will restart\n";
                $opts{restart} = 1;
                $fcn = 'stop';
            }
            if ($fcn eq 'stop') {
                $s .= "$0 will stop\n";
                $opts{stop} = 1;
            }
        }
        if ($opts{verbose}) { print "Restart/Stop: $s\n"; }
        SendText($conn, $s);
        return;
    }

    #   Handle df
    if ($path =~ '/df') {
        my $cmd = "$opts{df}";
        if ($opts{verbose}) { print "  CMD=$cmd\n"; }
        my $s = `$cmd  2>&1`;
        SendText($conn, $s);
        return;
    }

    #   Handle errorcheck
    if ($path =~ '/errorcheck') {
        my $cmd = "$opts{errorcheck}";
        if ($opts{verbose}) { print "  CMD=$cmd\n"; }
        my $s = `$cmd  2>&1`;
        SendText($conn, $s);
        return;
    }

    #   Handle logs
    if ($path =~ '/logs') {
        my $d = '/home/topmed/output';
        my $s = '';
        if (! chdir($d)) { $s = "Cannot CD to '$d': $!\n"; }
        else {
            foreach my $f (qw(
                topmed_init.log
                topmed_failures.log
                topmed_monitor_arrive.log
                topmed_monitor_b37.log
                topmed_monitor_bai.log
                topmed_monitor_cram.log
                topmed_monitor_expt.log
                topmed_monitor_ncbi.log
                topmed_monitor_orig.log
                topmed_monitor_phs.log
                topmed_monitor_qplot.log
                topmed_monitor_verify.log
                )) {
                $s .= "<b>Showing $f</b>\n" . `tail -6 $f`;
            }
            #   Show more details for logs we are really interested in at the moment
            foreach my $f (qw(
                topmed_monitor_ncbi.log
                )) {
                  $s .= "<pre><b>More Details on $f</b>\n" . `tail -60 $f` . "</pre>$f\n";
            }
        }
        SendText($conn, $s);
        return;
    }

    #   Nothing we knew about
    $conn->send_error(RC_NOT_ACCEPTABLE, "Invalid URI path '$path'");
    return;
}

#==================================================================
# Subroutine:
#   Clean($str)
#   Remove nasty chars from a string
#
# Arguments:
#   str - string
#
# Returns:
#   str
#==================================================================
sub Clean {
    my ($str) = @_;
    $str =~ s/\n/NL/g;
    $str =~ s/\t/ /g;
    $str =~ s/\r/CR/g;
    $str =~ s/;/SEMI/g;
    $str =~ s/\|/BAR/g;
    $str =~ s/\&/AMP/g;
    return $str;
}

#==================================================================
# Subroutine:
#   SendText($conn, $str)
#   Send $str as plain text on connection
#
# Arguments:
#   conn - connection object
#   str - string
#==================================================================
sub SendText {
    my ($conn, $str) = @_;

    $conn->send_response(
        HTTP::Response->new(RC_OK, undef,
            [
                'Content-Type' => 'plain/text',
                'Cache-Control' => 'no-store, no-cache, must-revalidate, post-check=0, pre-check=0',
                'Pragma' => 'no-cache',
                'Expires' => 'Thu, 01 Dec 1994 16:00:00 GMT',
            ],
            $str,
        )
    );
}

#==================================================================
#   Perldoc Documentation
#==================================================================
__END__

=head1 NAME

topmed_daemon.pl - Issue commands to the topmed server

=head1 SYNOPSIS

  topmed_daemon.pl

=head1 DESCRIPTION


Use this program to issue a SLURM squeue command and return
the results. This program is expected to run as a daemon 
usually with xinetd or inetd

This is necessary because our webservers are created with
a different distribution than the rest of our machines.
This means that the slurm binaries (like squeue) cannot
be run on the web server.

When we have all of our infrastructure machines on the same OS level,
then this should be removed from the topmed PHP code and done
more rationally. Until then, this is an embarrassing hack, but useful.

=head1 OPTIONS

=over 4

=item B<-fromips str>

Specifies a white space delimited string of IP addresses 
which will be allowed to issue requests to this server.

=item B<-help>

Generates this output.

=item B<-port N>

Specifies the port to listen on for requests.

=item B<-verbose>

Provided for developers to see additional information.

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

