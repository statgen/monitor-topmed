<?php
/*#################################################################
#
# Name: local_config.php
#
# Description:
#   Use this to customize monitors. Set $MON to the type
#   of monitor before including this. Allows this code to
#   provide the local config values needed for different monitors.
#
#  This i set up to use Twitter Bootstrap
#
# Copyright (C) 2010-2015 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/
$VERSION = '0.3';

//  Be sure I can see syntax errors
//ini_set('error_reporting', E_ALL ^ E_NOTICE);
//ini_set('error_reporting', E_ALL);
error_reporting(E_ALL);

//  These guys can do anything for any monitor action
$MGRS = array('tpg', 'scaron', 'schelcj');

//  Define fields for specific monitor functions
if ($MON == 'coredump') {
    $HDR['home'] = '/monitor/coredump';     // Banner URL
    $HDR['title'] = 'CoreDump Monitor';     // Banner title
    $HDR['logo'] = '/monitor/images/cdlogo_banner.png';     // Banner logo
    $HDR['local'] = "    <link href='/monitor/css/COREstyle.css' rel='stylesheet'>\n" .
        "    <script src='/monitor/js/COREjs.js'></script>\n";
    $HDR['footer'] = "<a class='footer' href=\"javascript:popup('cleanlog.php?fcn=rmlogs',640,480)\">Clean logs</a><br/>";

    $LDB['realm'] = '/usr/cluster/monitor/etc/.db_connections/monitor';
    $LDB['coredump'] = 'status';        // Table for coredump
    $LDB['events'] = 'events';          // Table for events for a run
    $LDB['hb'] = 'heartbeat';           // Table for heartbeat when things execute
    $LDB['statuses'] = 'statuses';      // Table for messages about some status

    //  These guys can make a few changes in the project information
    $PEDITOR = array('');
    //  These guys can request activities
    $MD5MGRS = array('tpg', 'scaron', 'rswillia', 'acoredump', 'gruberjd', 'hengw', 'svrieze', 'omergk');
}

if ($MON == 'topmed') {
    $HDR['home'] = '/monitor/topmed';        // Banner URL
    $HDR['title'] = 'NHLBI TOPMed - Data Tracking';   // Banner title
    $HDR['logo'] = '/monitor/images/nhlbilogo_banner.png';     // Banner logo
    $HDR['local'] = "    <link href='/monitor/css/TOPMEDstyle.css' rel='stylesheet'>\n" .
        "    <script src='/monitor/js/TOPMEDjs.js'></script>\n";

    $LDB['realm'] = '/usr/cluster/monitor/etc/.db_connections/topmed';
    $LDB['centers'] = 'centers';            // SQL table names
    $LDB['runs'] = 'runs';
    $LDB['bamfiles'] = 'bamfiles';
    $LDB['pulls'] = 'requestfiles';

    $FILES['topdir'] = '/net/topmed/incoming/topmed';

    //  SLURM partitions we might use
    $SLURM['qlocal'] = 'topmed-incoming topmed2-incoming';
    $SLURM['qclients'] = 'topmed nomosix';

    //  These guys can request activities
    $REQMGRS = array('tblackw');
}


if ($MON == 'backups') {
    $HDR['home'] = '/monitor/backups';      // Banner URL
    $HDR['title'] = 'RSNAP - Last Successful Daily Backup';       // Banner title
    $HDR['logo'] = '/monitor/images/backlogo_banner.png';     // Banner logo
    $HDR['local'] = "    <link href='/monitor/css/RSNAPstyle.css' rel='stylesheet'>\n" .
        "    <script src='/monitor/js/RSNAPjs.js'></script>\n";

    $TOPDIR = '/net/statgen/statgen/monitor/backups';
    $CLUSTERDIR = '/usr/cluster/adm/rsnap';

    $FILES['groupfile'] = $CLUSTERDIR . '/host.list';
    $FILES['calcsecs'] = $TOPDIR . '/calcsecs.pl';

    //  Defaults for connecting to database
    $LDB['realm'] = $CLUSTERDIR .  '/rsnap.realm';
    $LDB['table']= 'rsnap';
}

//================================================================
//     Common for all monitor programs
//================================================================
//  Here are common META tags
$HDR['meta'] = <<<END
    <meta name='description' content='Monitor CSG Data'/>
    <meta name='robots' content='index,follow'/>
    <meta name='resource-type' content='document'/>
    <meta http-equiv='expires' content='0'/>
    <meta name='author' content='Terry Gliedt and others, University of Michigan'/>
    <meta name='copyright' content='Copyright (c) 2009- by the University of Michigan'/>
    <meta name='keywords' content='track monitor backups genonmic data delivery'/>

END;
//  Add Twitter Bootstrap to the local CSS/JS
$tb = <<<END
    <link href='/monitor/css/bootstrap.min.css' rel='stylesheet'>
    <script src='https://oss.maxcdn.com/html5shiv/3.7.2/html5shiv.min.js'></script>
    <script src='https://oss.maxcdn.com/respond/1.4.2/respond.min.js'></script>
    <script src='https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js'></script>
    <script src='/monitor/js/bootstrap.min.js'></script>

END;
$HDR['local'] = $tb . $HDR['local'];

//  If footer already defined, insert it in the standard footer
if (isset($HDR['footer'])) { $f = $HDR['footer']; }
else { $f = ''; }
$HDR['footer'] = <<<END
Copyright (c) 2010- University of Michigan<br/>
Report problems to <a class='footer' href='mailto:tpg@umich.edu'>tpg@umich.edu</a><br/>
$f
<br/>Last Revision: Apr 2015 <br/>

END;


?>
