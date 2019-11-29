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
#  This is set up to use Twitter Bootstrap
#
# Copyright (C) 2010-2015 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/
$VERSION = '1.0';

//  Be sure I can see syntax errors
//ini_set('error_reporting', E_ALL ^ E_NOTICE);
//ini_set('error_reporting', E_ALL);
error_reporting(E_ALL);

//  These guys can do anything for any monitor action
$MGRS = array('tpg');
//  These guys can request activities
$REQMGRS = array('tblackw');

$HDR['home'] = '/topmed';        // Banner URL
$HDR['title'] = 'NHLBI TOPMed - Data Tracking';   // Banner title
$HDR['logo'] = '/topmed/images/nhlbilogo_banner.png';     // Banner logo
$HDR['local'] = "    <link href='/topmed/css/TOPMEDstyle.css' rel='stylesheet'>\n" .
    "    <script src='/topmed/js/TOPMEDjs.js'></script>\n";
$LDB['bindir'] = '/usr/cluster/topmed/bin';
$LDB['SLURMNODES'] = '/net/topmed/working/topmed-output/slurm.nodes';
$LDB['centers'] = 'centers';
$LDB['permissions'] = 'permissions';

//	Figure out what kind of data we are to display
if (strstr($_SERVER['PHP_SELF'],'rnaseq')) {
	$LDB['datatype'] = 'rnaseq';
	$LDB['realm'] = $LDB['bindir'] . '/../etc/.db_connections/rnaseq';
	$LDB['datatype'] = 'rnaseq';
	$LDB['projects'] = 'tx_projects';   	// Nicknames for tables and pkeys
	$LDB['projects_pkey'] = 'rnaprojectid';
	$LDB['samples'] = 'tx_samples';
	$LDB['samples_pkey'] = 'txseqid';
	$LDB['files'] = 'tx_files';
	$LDB['files_pkey'] = 'fileid';
}
if (strstr($_SERVER['PHP_SELF'],'genome') || strstr($_SERVER['PHP_SELF'],'index')) {
	$LDB['datatype'] = 'genome';
	$LDB['realm'] = $LDB['bindir'] . '/../etc/.db_connections/topmed';
	$LDB['runs'] = 'runs';
	$LDB['runs_pkey'] = 'runid';
	$LDB['samples'] = 'bamfiles';
	$LDB['samples_pkey'] = 'bamid';
	$LDB['stepstats'] = 'stepstats';
}

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
    <link href='/topmed/css/bootstrap.min.css' rel='stylesheet'>
    <script src='https://oss.maxcdn.com/html5shiv/3.7.2/html5shiv.min.js'></script>
    <script src='https://oss.maxcdn.com/respond/1.4.2/respond.min.js'></script>
    <script src='https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js'></script>
    <script src='/topmed/js/bootstrap.min.js'></script>

END;
$HDR['local'] = $tb . $HDR['local'];

//  If footer already defined, insert it in the standard footer
if (isset($HDR['footer'])) { $f = $HDR['footer']; }
else { $f = ''; }
$HDR['footer'] = <<<END
Copyright (c) 2009- University of Michigan<br/>
Report problems to <a class='footer' href='mailto:tpg@umich.edu'>tpg@umich.edu</a><br/>
$f
<br/>Last Revision: June 2019 <br/>

END;


?>
