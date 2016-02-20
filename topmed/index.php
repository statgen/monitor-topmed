<?php
/*#################################################################
#
# Name: index.php
#
# Description:
#   Display information for tracking NHLBI TopMed data
#
# Copyright (C) 2015 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/
$MON='topmed'; include_once 'local_config.php';
include_once 'common.php';
include_once 'header.php';
include_once 'DBMySQL.php';
include_once "edit.php";

//print "<!-- _POST=\n"; print_r($_POST); print " -->\n";

// Hack to tell if we are on statgen or not
if ( file_exists('/exports/monitor/index.html')) { $oldsite = 1; }
else { $oldsite = 0; }

$qurl =  $_SERVER['SCRIPT_NAME'] . "?fcn=queue'";
$STATUSLETTERS =  "<i><b>A</b>=File Arrived, <b>5</b>=MD5 Verified, <b>C</b>=BAM=>CRAM, <b>I</b>=BAI created<br/>" .
    "<b>Q</b>=qplot run, <b>7</b>=Remapped Build=37, <b>8</b>=Remapped Build=38, <b>X</b>=EXPT=>NCBI<br/>" .
    "<b>S</b>=Secondary BAM=>NCBI, <b>P=</b>B37=>Primary BAM=>NCBI, <b>T</b>=b38=>Tertiary BAM=>NCBI";

$SHOWQUEUES = "STATUS: &nbsp;&nbsp;&nbsp;" .
    "<a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=showqlocal' " .
    "onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=showqlocal\",680,720); " .
    "return false;'>Local Queue</a> &nbsp;&nbsp;&nbsp;" .

    //"<a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=showqrmt' " .
    //"onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=showqrmt\",680,720); " .
    //"return false;'>Remote Queue</a> &nbsp;&nbsp;&nbsp;" .

    "<a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=df' " .
    "onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=df\",680,720); " .
    "return false;'>Disk Usage</a> &nbsp;&nbsp;&nbsp;";

    //"<a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=errorcheck' " .
    //"onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=errorcheck\",680,720); " .
    //"return false;'>Error Check</a> &nbsp;&nbsp;&nbsp;" .

    //"<a href='/topmed/loadedfiles.txt' target='blank'> " .
    //" NCBI Log</a> &nbsp;&nbsp;&nbsp;" .

    if ($oldsite) { 
        $SHOWQUEUES .= "<a href='/monitor/topmed/plot.php' target='plots'> " .
            " Plots</a> &nbsp;&nbsp;&nbsp;" .
            "<a href='http://gcsdev.sph.umich.edu/topmed/' target='newcode'> New Site </a> &nbsp;&nbsp;&nbsp;";
    }
    else {
        $SHOWQUEUES .= "<a href='/topmed/plot.php' target='plots'> " .
            " Plots</a> &nbsp;&nbsp;&nbsp;";
    }

    $SHOWQUEUES .="<a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=logs' " .
    "onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=logs\",680,720); " .
    "return false;'>Tail Logs</a>";

$BAMNOTE = "<p><b>Note:</b><br>" .
    "<i>Colors:</i> " .
    "<span class='cancelled'> cancelled </span>&nbsp;" .
    "<span class='completed'> completed </span>&nbsp;" .
    "<span class='delivered'> delivered </span>&nbsp;" .
    "<span class='failed'> failed </span>&nbsp;" .
    "<span class='requested'> requested </span>&nbsp;" .
    "<span class='started'> started </span>&nbsp;" . 
    "<span class='submitted'> submitted </span>&nbsp;" .
    "<br>" .
    "<i>Status:</i> state of tasks ($STATUSLETTERS)<br>" .
    "</p>\n";
$RUNNOTE = "<p><b>Note:</b><br>" .
    "<i>Status of <b>all BAMs</b> in run:</i>&nbsp;" .
    "<span class='done'> finished </span>&nbsp;" .
    "<span class='processing'> being processed </span>&nbsp;" .
    "<span class='unknown'> unstarted, etc </span>&nbsp;" .
    "<span class='failed'> at least one failure </span>&nbsp;" .
    "<br>" .
    "<i>Steps:</i> $STATUSLETTERS<br>" .
    "</p>\n";

$quickcols = array(                     // Map of status column to topmedcmd verb
    'state_arrive'   => 'arrived',
    'state_md5ver'   => 'md5ver',
//    'state_backup'   => 'backedup',
    'state_bai'      => 'bai',
    'state_qplot'    => 'qplot',
    'state_cram'     => 'cramed',
    'state_b37'      => 'mapping37',
    'state_b38'      => 'mapping38',
    'state_ncbiexpt' => 'sendexpt',
    'state_ncbiorig' => 'sendorig',
    'state_ncbib37'  => 'sendb37',
    'state_ncbib38'  => 'sendb38'
);
$quickletter = array(                   // Map of status column to letter we see
    'state_arrive'   => 'A',
    'state_md5ver'   => '5',
    'state_bai'      => 'I',
    'state_qplot'    => 'Q',
    'state_cram'     => 'C',
    'state_b37'      => '7',
    'state_b38'      => '8',
    'state_ncbiexpt' => 'X',
    'state_ncbiorig' => 'S',
    'state_ncbib37'  => 'P',
    'state_ncbib38'  => 'T'
);
$validfunctions = array('all', 'verify', 'bai', 'qplot', 'cram', 'sexpt', 'sorig', 'sb37', 'sb38');
$NOTSET = 0;                // Not set
$REQUESTED = 1;             // Task requested
$SUBMITTED = 2;             // Task submitted to be run
$STARTED   = 3;             // Task started
$DELIVERED = 19;            // Data delivered, but not confirmed
$COMPLETED = 20;            // Task completed successfully
$CANCELLED = 89;            // Task cancelled
$FAILED    = 99;            // Task failed
$state2str = array(         // Values here are class for SPAN tag
    $NOTSET => 'notset',
    $REQUESTED => 'requested',
    $SUBMITTED => 'submitted',
    $STARTED => 'started',
    $DELIVERED => 'delivered',
    $COMPLETED => 'completed',
    $CANCELLED => 'cancelled',
    $FAILED => 'failed'
);

$NODELIST = array('topmed', 'topmed2', 'topmed3', 'topmed4');
$TOPMEDJOBNAMES = array('verify', 'bai', 'qplot', 'cram', 'expt', 'orig', 'b37', 'b38');
$LOGFILES = array('topmed_init.log',
    'topmed_monitor_arrive.log',
    'topmed_monitor_bai.log',
    'topmed_monitor_cram.log',
    'topmed_monitor_expt.log',
    'topmed_monitor_ncbi.log',
    'topmed_monitor_orig.log',
    'topmed_monitor_phs.log',
    'topmed_monitor_qplot.log',
    'topmed_monitor_verify.log',
    'topmed_monitor_b37.log' );

//  These columns are state values to be converted to people readable strings
//  See DateState() for possible values
$conv2dat = array_keys($quickcols);

//  Handy javascript fragments
$JS['VSPACER'] = "<br/><br/><br/><br/>";
$JS['SPACER'] = "&nbsp;&nbsp;&nbsp;&nbsp;";
$JS['CLOSE'] = "<p align='right'><font size='-1'><a href='javascript:window.close()'>Close</a>" . $JS['SPACER'];
$JS['BACK'] = "<p align='right'><font size='-1'><a href='javascript:history.back()'>Back</a>" . $JS['SPACER'];

//  Fixed values for status in requestfiles database
$REQ_STATUS['QUEUED'] = 'request_queued';

//print "<!-- _POST=\n"; print_r($_POST); print " -->\n";

//-------------------------------------------------------------------
//  Get parameters passed in via normal invocation
//-------------------------------------------------------------------
if (! isset($_SERVER['REMOTE_USER'])) { $_SERVER['REMOTE_USER'] = 'none'; }
$iammgr = 0;
if (in_array($_SERVER['REMOTE_USER'], $MGRS)) { $iammgr = 1; }
if (in_array($_SERVER['REMOTE_USER'], $REQMGRS)) { $iammgr = 1; }
//  If a manager, allow them to control job submission
if ($iammgr) {
$SHOWQUEUES .= " &nbsp;&nbsp;&nbsp; <a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=control' " .
    "onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=control\",680,720); " .
    "return false;'>Control Jobs</a>";
}
print doheader($HDR['title'], 1);

if ($iammgr) {
    $s = "<b>See TOPMed monitor docs " .
        "<a href='https://statgen.sph.umich.edu/wiki/NHLBI_automation_steps' target='_blank'>" .
        "here</a>.</b>\n";
    if (! $oldsite) {
        $s .= "<font color='red'>New Web Site, Not Everything Works</font></br>\n";
    }
}
else {
    $s = '';
    if (! $oldsite) {
        $s .= "<font color='red'>New Web Site, Not Everything Works</font></br>\n";
    }
}
print "<p class='intro'>The <a href='http://www.nhlbi.nih.gov/'>NHLBI</a> provides " .
    "science-based, plain-language information related to heart, lung " .
    "and blood diseases and conditions and sleep disorders.\n" .
    "Details about this data tracking are available from Tom Blackwell " .
    "(<a href='mailto:tblackw@umich.edu'>tblackw@umich.edu</a>). $s</p>\n";

//  Real parameters for each form, default is ''
$parmcols = array('fcn', 'maxdir', 'desc', 'center',
    'run', 'runid', 'bamid', 'centerid', 'fetchpath', 'hostname', 'col',
    'op', 'id');
extract (isolate_parms($parmcols));
if (! $center) { $center = 'all'; }
if (! $fcn)    { $fcn = 'runs'; }
if (! $maxdir) { $maxdir = '150'; }

DB_Connect($LDB['realm']);
GetCenters();                   // Get maps to identify centers
$HTML = '';

if ($fcn == 'queue') {              // This cannot work until statgen can run squeue :-(
    $q = 'topmed-incoming';
    print "<h3 align='center'>SLURM Queue:  $q</h3>\n" .
        "<div class='indent'><table width='80%'><pre>\n";
    print shell_exec("/usr/cluster/bin/squeue -p $q");
    print "<p align='right'><font size='-1'><a href='javascript:window.close()'>Close</a>&nbsp;&nbsp;&nbsp;</p>\n";
    print dofooter($HDR['footer']);
    exit;
}
if ($fcn == 'runs') {
    print ViewRuns($center, $maxdir, $iammgr);
    print dofooter($HDR['footer']);
    exit;
}
if ($fcn == 'rundetail') {
    print ViewRunDetail($runid);
    print dofooter($HDR['footer']);
    exit;
}
if ($fcn == 'bams') {
    print ViewBams($runid, $maxdir, $iammgr);
    print dofooter($HDR['footer']);
    exit;
}
if ($fcn == 'bamdetail') {
    print ViewBamDetail($bamid);
    print dofooter($HDR['footer']);
    exit;
}
if ($fcn == 'reqshow') {
    $a = 1;
    if ($iammgr) { $a = 0; }
    print ViewPullQueue($centerid, $a);
    print dofooter($HDR['footer']);
    exit;
}
if ($fcn == 'permit') {
    $h = HandlePermit($op, $center, $run, $id);
    print ControlJobs($h);
    print dofooter($HDR['footer']);
    exit;
}
if ($fcn == 'control') {
    print ControlJobs('');
    print dofooter($HDR['footer']);
    exit;
}

if ($fcn == 'showqlocal') {
    print "<center>$SHOWQUEUES &nbsp;&nbsp;&nbsp;</center>\n";
    if ($oldsite) {
        foreach ($NODELIST as $n) {
            $cmd = "/usr/cluster/monitor/bin/slurm_query.sh -squeue ${n}-incoming";
            print `$cmd`;
        }
        $cmd = "/usr/cluster/monitor/bin/slurm_query.sh -squeue topmed-working";
        print `$cmd`;
        exit;
    }
    foreach ($NODELIST as $n) { print ShowSLURMIcoming($n); }
    exit;
}

if ($fcn == 'df') {
    print "<center>$SHOWQUEUES &nbsp;&nbsp;&nbsp;</center>\n";
    if ($oldsite) {
        $cmd = "/usr/cluster/monitor/bin/slurm_query.sh -df ignored";
        print `$cmd`;
        exit;
    }
    $cmd = "df -h";
    foreach ($NODELIST as $n) { $cmd .= " /net/$n/incoming /net/$n/working";  }
    print "<pre>" . `$cmd` . "</pre>\n";
    exit;
}

if ($fcn == 'logs') {
    print "<center>$SHOWQUEUES &nbsp;&nbsp;&nbsp;</center>\n";
    if ($oldsite) {
        $cmd = "/usr/cluster/monitor/bin/slurm_query.sh -logs ignored";
        print `$cmd`;
        exit;
    }
    $d = '/net/topmed/working/topmed-output';
    if (! chdir($d)) {
        print Emsg("Unable to find logs in '$d'", 1) . "<br/>\n";
        exit;
    }
    $s = '<pre>';
    foreach ($LOGFILES as $f) { $s .= "<b>Showing $f</b>\n" . `tail -6 $f`; }
    print "$s</pre>\n";
    exit;
}

//  Certain users can do more
if ($iammgr) {
    if ($fcn == 'setcancel') {
        $sql = 'UPDATE ' . $LDB['bamfiles'] . " SET $col=$CANCELLED WHERE bamid=$bamid";
        $result = SQL_Query($sql);
        Msg("Cancelled '$col' for BAM '$bamid'  SQL=$sql");
        print $JS['CLOSE'];
        print dofooter($HDR['footer']);
        exit;
    }
    if ($fcn == 'setrequest') {
        $sql = 'UPDATE ' . $LDB['bamfiles'] . " SET $col=$REQUESTED WHERE bamid=$bamid";
        $result = SQL_Query($sql);
        Msg("Requested '$col' for BAM '$bamid'");
        print $JS['CLOSE'];
        print dofooter($HDR['footer']);
        exit;
    }
    if ($fcn == 'reqfetch') {
        print PromptFetch($centerid);
        print dofooter($HDR['footer']);
        exit;
    }
    if ($fcn == 'reqsave') {
        print QueueFetch($centerid, $hostname, $fetchpath);
        print dofooter($HDR['footer']);
        exit;
    }
    
    if ($fcn == 'editrun') {
        $sql = 'SELECT * FROM ' . $LDB['runs'] . " WHERE runid='$runid'";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);

        $HTML .= "<h3 align='center'>Edit Entry for " . $row['dirname'] . "</h3>\n" .
            '<center>' . Emsg("Just because you CAN change a field does not mean you should", 1) . "<br/></center>\n" .
            "<div class='indent'>\n";
        $HTML .= Edit($LDB['runs'], $row) . "</div>\n";
        print $HTML;
        print dofooter($HDR['footer']);
        exit;
    }
    if ($fcn == 'editbam') {
        $sql = 'SELECT * FROM ' . $LDB['bamfiles'] . " WHERE bamid='$bamid'";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);

        $HTML .= "<h3 align='center'>Edit Entry for " . $row['bamname'] . "</h3>\n" .
            '<center>' . Emsg("Just because you CAN change a field does not mean you should", 1) . "<br/></center>\n" .
            "<div class='indent'>\n";
        $HTML .= Edit($LDB['bamfiles'], $row) . "</div>\n";
        print $HTML;
        print dofooter($HDR['footer']);
        exit;
    }
    if ($fcn == 'modify') {
        $id = $runid;
        $t = $LDB['runs'];
        if ($bamid) { $id = $bamid; $t = $LDB['bamfiles']; }
        if (Modify($t, $id)) { $HTML .= Msg("Record '$id' modified", 1); }
        else { $HTML .= Emsg("Failed to modify record '$id'", 1); }
        print $HTML .
            "<p align='right'><font size='-1'><a href='javascript:window.close()'>Close</a>&nbsp;&nbsp;&nbsp;</p>\n";
        print dofooter($HDR['footer']);
        exit;
    }
}

//  What was that?
Emsg("Unknown directive '$fcn'.");
Nice_Exit("How'd you do that?");
exit;

/*---------------------------------------------------------------
#   html = PromptFetch($cid)
#   Prompt for a fetch request of remote data
---------------------------------------------------------------*/
function PromptFetch($cid) {
    global $CENTERID2NAME, $JS;

    //  Verify the center is known
    if (! isset($CENTERID2NAME[$cid])) {
        return Emsg("Center '$cid' is unknown. How'd you do that?", 1) . $JS['CLOSE'];
    }
    $center = strtoupper($CENTERID2NAME[$cid]);
    $html = "<h3 align='center'>Request Pull of a Run From '$center'</h3>\n" .
        "<p>Filling out this form will initiate a pull of your data " .
        "by NHLBI TopMed at the University of Michigan. Be careful to provide the correct " .
        "hostname and path to your Run.</p>\n";

    //  Now generate form so user can queue a fetch of remote data
    $html .= "<form action='" . $_SERVER['PHP_SELF'] . "?fcn=reqsave' method='post'>\n" .
        "<input type='hidden' name='centerid' value='$cid'>\n";
    $html .= "<b>Hostname</b>&nbsp;&nbsp;" .
        "<input type='text' name='hostname' size='24'>" .
        "<br/><font size='-1' color='green'>Fully qualified hostname to fetch data from</font></td>" .
        "<br/><br/>" .
        "<b>Fetch Path</b>&nbsp;&nbsp;</td>" .
        "<input type='text' name='fetchpath' size='70'>" .
        "<br/><font size='-1' color='green'>Fully qualified path to Run directory on host</font></td>" .
        "<br/><br/>" .
        "<input type='submit' value=' Queue Request to Fetch Data '>\n" .
        "<input type='button' name='close' value=' Close ' onclick='javascript:window.close()'>\n" .
        "</form>\n";
    $html .= $JS['VSPACER'] . $JS['CLOSE'];
    return $html;
}

/*---------------------------------------------------------------
#   html = QueueFetch($cid, $host, $path)
#   Verify and queue a request to fetch remote data
---------------------------------------------------------------*/
function QueueFetch($cid, $host, $path) {
    global $LDB, $CENTERID2NAME, $JS, $REQ_STATUS;

    //  Verify the center is known
    if (! isset($CENTERID2NAME[$cid])) {
        return Emsg("Center '$cid' is unknown. How'd you do that?", 1) . $JS['CLOSE'];
    }
    $center = strtoupper($CENTERID2NAME[$cid]);
    $html = "<h3 align='center'>Queue Pull of a Run From '$center'</h3>\n";

    //  Do basic sanity check on input
    if (! preg_match('/^\//', $path)) {
        return Emsg("Path '$path' is not fully qualified (must start with '/'.", 1) . $JS['BACK'];

    }
    $hostip = gethostbyname($host);
    if ($hostip == $host) {
        return Emsg("Hostname '$host' cannot be resolved. Invalid host", 1) . $JS['BACK'];
    }

    //  Make sure this request has not been made before
    $sql = 'SELECT * FROM ' . $LDB['pulls'] . " WHERE centerid=$cid AND " .
        "fetchpath='$path'";
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);
    if ($numrows) {             // Yikes, this already exists
        return Emsg("Duplicate request: '$path' already queued. Nothing done.", 1) .
         $JS['VSPACER'] . $JS['BACK'];
    }
    $d = date('Y/m/d');
    $sql = 'INSERT INTO ' . $LDB['pulls'] .
        ' (centerid,hostname,fetchpath,daterequestor,iprequestor,ipuser,status) ' .
        "VALUES($cid,'$host','$path','$d','" . $_SERVER['REMOTE_ADDR'] . "','" .
        $_SERVER['REMOTE_USER'] . "','" . $REQ_STATUS['QUEUED'] . "')";
    $result = SQL_Query($sql);
    $html .= "<p>Request for '$path' from host '$host' has been queued</p>\n" .
         $JS['VSPACER'] . $JS['CLOSE'];
    return $html;
}

/*---------------------------------------------------------------
#   html = ViewPullQueue($cid, $onlyactive)
#   Verify and queue a request to fetch remote data
---------------------------------------------------------------*/
function ViewPullQueue($cid, $onlyactive=1) {
    global $LDB, $CENTERID2NAME, $JS, $REQ_STATUS;

    //  Verify the center is known
    if (! isset($CENTERID2NAME[$cid])) {
        return Emsg("Center '$cid' is unknown. How'd you do that?", 1) . $JS['CLOSE'];
    }
    $center = strtoupper($CENTERID2NAME[$cid]);
    if ($onlyactive) { $a = 'Active'; }
    else { $a = 'All'; }
    $html = "<h3 align='center'>Queue of $a Pull Runs From '$center'</h3>\n";

    //  Make sure this request has not been made before
    $sql = 'SELECT * FROM ' . $LDB['pulls'] . " WHERE centerid=$cid";
    if ($onlyactive) { $sql .= "AND status='" . $REQ_STATUS['QUEUED'] . "'"; }
    $html .= ShowTable($LDB['pulls'], $sql);
    return $html;
}

/*---------------------------------------------------------------
#   html = ViewRuns($center, $maxdirs, $iammgr)
#   Show summary of directories of runs
---------------------------------------------------------------*/
function ViewRuns($center, $maxdirs, $iammgr) {
    global $LDB, $CENTERS, $CENTERID2NAME, $CENTERNAME2ID, $RUNNOTE, $SHOWQUEUES;
    $hdrcols  = array('dirname', 'status', 'bamcount', 'dateinit', 'datecomplete');

    //  Generate HTML header for page
    //  Get list of centers doing:  select distinct(project) from status;
    $html = "<h3 align='center'>Summary of Last $maxdirs Directories</h3>\n" .
        "<center>" .
        "<a href='" . $_SERVER['SCRIPT_NAME'] . "?maxdirs=25'>Show Last 25</a>" .
        "&nbsp;&nbsp;&nbsp; or &nbsp;&nbsp;&nbsp;" .
        "<a href='" . $_SERVER['SCRIPT_NAME'] . "?maxdir=75'>Show Last 75</a>" .
        "&nbsp;&nbsp;&nbsp; or &nbsp;&nbsp;&nbsp;" .
          "<a href='" . $_SERVER['SCRIPT_NAME'] . "?" . $_SERVER['QUERY_STRING'] .
        "'>Refresh</a><br>\n Or Choose a Center: " .
        "&nbsp;&nbsp;<a href='" . $_SERVER['SCRIPT_NAME'] .
        "?center=all&amp;maxdir=$maxdirs'>All</a>&nbsp;&nbsp;\n";
    foreach ($CENTERS as $c) {
        $html .= "&nbsp;&nbsp;<a href='" . $_SERVER['SCRIPT_NAME'] .
            "?center=$c&amp;maxdir=$maxdirs'>" . strtoupper($c) . "</a>&nbsp;&nbsp;\n";
    }
    $html .= "<br/>$SHOWQUEUES</center>\n";

    $centers2show = array();                // Get list of centers for this query
    if ($center == 'all') { $centers2show = $CENTERS; }
    else { array_push($centers2show, $center); }

    //  For each center, show details from database ($rows)
    foreach ($centers2show as $centr) {
        $cid = $CENTERNAME2ID[$centr];
        $html .= "<br><div class='indent'><b>" . strtoupper($centr);
        if ($iammgr) {
            $u = $_SERVER['SCRIPT_NAME'] ."?fcn=reqfetch&amp;centerid=$cid";
            $html .= "&nbsp;&nbsp;&nbsp;<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
                "<font color='blue' size='-1'>Queue Pull</font></a>";
            $u = $_SERVER['SCRIPT_NAME'] ."?fcn=reqshow&amp;centerid=$cid";
            $html .= "&nbsp;&nbsp;&nbsp;<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
                "<font color='blue' size='-1'>Show Queue</font></a>";
        }
        $html .= "</b></div>\n";
        //  Build start of table for each center
        $html .= "<table align='center' width='100%' border='1'><tr>\n";
        foreach ($hdrcols as $c) {
            $html .= "<th class='heading'>" . ucfirst($c) . "</th>\n";
        }
        if ($iammgr) { $html .= "<th>&nbsp;</th>"; }
        $html .= "</tr>\n";

        //  Walk through database getting data for this center
        $sql = 'SELECT * FROM ' . $LDB['runs'] . " WHERE centerid=$cid";
        if ($maxdirs) { $sql .= " LIMIT $maxdirs"; }
        $result = SQL_Query($sql);
        $numrows = SQL_NumRows($result);
        $rows = array();            // Save DB info for later display
        $centers2show = array();    // Get list of centers for this query
        for ($i=0; $i<$numrows; $i++) {
            $row = SQL_Fetch($result);
            $rows[$row['runid']] = $row;
        }

        reset($rows);
        foreach ($rows as $id => $row) {
            if ($row['centerid'] != $cid) { continue; }
            //print "<!-- CID=$cid ROW=\n"; print_r($row); print " -->\n";
            //  Show data for this center
            $html .= "<tr>\n";
            reset($hdrcols);
            foreach ($hdrcols as $c) {
                $d = $row[$c];
                if ((! isset($d)) || ($d == '')) { $d = '&nbsp;'; }
                if ($c == 'dirname') {
                    $u = $_SERVER['SCRIPT_NAME'] . "?fcn=bams&amp;runid=" . $row['runid'];
                    $d = "<a href='$u'>$d</a>";
                }
                if ($c == 'status') { $d = CalcRunStatus($d); }
                //if ($c == 'dateinit' && (! preg_match('/\D/', $d))) { $d = date('Y/m/d H:i', $d); }
                if ($c == 'dateinit') { $d = date('Y/m/d H:i', $d); }
                if ($c == 'datecomplete' && $d != '&nbsp;') { $d = date('Y/m/d H:i', $d); }
                $html .= "<td align='center'>$d</td>\n";
            }
            
            $html .= "<td align='center'>";
            $u = $_SERVER['SCRIPT_NAME'] ."?fcn=rundetail&amp;runid=" . $row['runid'];
            $html .= "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
                    "<font color='green' size='-1'>Details</font></a>&nbsp;&nbsp;&nbsp;";
            if ($iammgr) {
                $u = $_SERVER['SCRIPT_NAME'] ."?fcn=editrun&amp;runid=" . $row['runid'];
                $html .= "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
                    "<font color='red' size='-1'>Edit</font></a>";
            }
            $html .= "</td></tr>\n";
        }
        $html .= "</table>\n";
    }
    $html .= "<div class='indent'>\n$RUNNOTE</div>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = ViewRunDetail($runid)
#   Show details for a particular run
---------------------------------------------------------------*/
function ViewRunDetail($runid) {
    global $LDB, $CENTERS, $CENTERID2NAME, $CENTERNAME2ID;

    //  Show all we have for this directory
    $sql = 'SELECT * FROM ' . $LDB['runs'] . " WHERE runid=$runid";
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);
    if ($numrows != 1) {
        Emsg("Unable to find detail for run '$runid'");
        return '';
    }
    $row = SQL_Fetch($result);
    $dir = $row['dirname'];

    //  Now generate the details section for the directory
    $html = "<h3 align='center'>Details for '$dir'</h3>\n" .
        "<div class='indent'><table width='80%'>";
    $collist = array_keys($row);
    foreach ($collist as $c) {
        if ($c == 'runid') { continue; }
        $d = $row[$c];
        if ($c == 'centerid') { $d = "$d (" . $CENTERID2NAME[$d] . ')'; }
        if ((! isset($d)) || (is_null($d)) || ($d == '') || ($d == 'NULL')) { $d = 'n/a'; }
        //  If the data has a newline in it do more than just show it
        $p = strpos($d, "\n");
        if ($p) { $d = "<pre>$d</pre>"; }
        $html .= "<tr><td align='left' valign='top' width='15%'><b>" . ucfirst($c) . ":</b></td>" .
            "<td>&nbsp;&nbsp;&nbsp;</td>";
        $html .= "<td align='left' colspan='4'>$d</td>";
        $html .= "</tr>\n";
    }
    $html .= "</table>\n" .
        "<p align='right'><font size='-1'><a href='javascript:window.close()'>Close</a>&nbsp;&nbsp;&nbsp;</p>\n" .
        "</div>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = ViewBams($runid, $maxdirs, $iammgr)
#   Show list of BAM files for a particular run
---------------------------------------------------------------*/
function ViewBams($runid, $maxdirs, $iammgr) {
    global $LDB, $HDR, $CENTERS, $CENTERID2NAME, $CENTERNAME2ID, $FILES, $BAMNOTE, $SHOWQUEUES;
    $hdrcols  = array('bamname', 'QUIKSTAT', 'bamsize', 'piname');
    $html = '';
    $maxdirs = 0;                   // For now, show all BAMs

    //  Figure out RUN for this BAM
    $sql = 'SELECT dirname,centerid,bamcount FROM ' . $LDB['runs'] . " WHERE runid=$runid";
    $result = SQL_Query($sql, 0);
    $e = DB_CheckError($result);
    if ($e) { return EMsg("Surprise: Unable to find RUN id=$runid: $e", 1); }
    $row = SQL_Fetch($result);
    $runname = $row['dirname'];
    $bamcount = $row['bamcount'];
    $centername = $CENTERID2NAME[$row['centerid']];
    $center = strtoupper($centername);

    //  Generate HTML header for page
    $sql = 'SELECT * FROM ' . $LDB['bamfiles'] . " WHERE runid=$runid";
    if ($maxdirs) { $sql .= " LIMIT $maxdirs"; }
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);

    //  Show details for each bam
    $url = $HDR['home'] . "/index.php?center=$centername&amp;maxdir=50";
    $html .= "<h3 align='center'>$bamcount BAM Files for '$runname' [$runid] in center " .
        "<a href='$url'>$center</a></h3>\n";
    $html .= "<p align='center'/>$SHOWQUEUES</p>\n";

    $html .= "<table align='center' width='100%' border='1'><tr>\n";
    foreach ($hdrcols as $c) {
        $html .= "<th class='heading'>" . ucfirst($c) . "</th>\n";
    }
    if ($iammgr) { $html .= "<th>&nbsp;</th>"; }
    $html .= "</tr>\n";

    for ($i=0; $i<$numrows; $i++) {
        $row = SQL_Fetch($result);
        reset($hdrcols);
        foreach ($hdrcols as $c) {
            if ($c == 'QUIKSTAT') { $d = QuickStatus($row); }
            else { $d = $row[$c]; }
            if ((! isset($d)) || ($d == '')) { $d = '&nbsp;'; }
            if ($c == 'bamname') {
                $u = $_SERVER['SCRIPT_NAME'] . "?fcn=bamdetail&amp;bamid=" . $row['bamid'];
                $d = "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>$d</a>";
            }
            if ($c == 'bamsize') { $d = ShortForm($d); }
            $html .= "<td align='center'>$d</td>\n";
        }           
        if ($iammgr) {
            $u = $_SERVER['SCRIPT_NAME'] ."?fcn=editbam&amp;bamid=" . $row['bamid'];
            $html .= "<td>&nbsp;<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
                "<font color='red' size='-1'>Edit</font></a></td>";
        }
        $html .= "</tr>\n";
    }
    $html .= "</table>\n" .
        "<div class='indent'>\n$BAMNOTE</div>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = ShowSLURMIcoming($host)
#   Show summary of SLURM queues for this host
---------------------------------------------------------------*/
function ShowSLURMIcoming($host) {
    global $TOPMEDJOBNAMES;
    $JOBSTOSHOW = 10;

    //  Get state for this machine
    $cmd = "/usr/cluster/bin/scontrol show node $host | grep State=";
    $results = `$cmd`;
    $nodestate = 'State=Node not known';
    $nodestate = 'State=Fix scontrol, Sean';

    #   This should work too, returning something like 'State=MIXED'
    #   except gcsdev does not work. Somehow the 48109 port is blocked maybe?
    #$results = `/usr/bin/wget -o /dev/null -O $tmpfile.0 $url/shownode/$h`;

    if (preg_match('/(State=\S+)/', $results, $m)) { $nodestate = $m[1]; }
    if (preg_match('/DRAIN/', $nodestate)) { $nodestate = "<font color=red>$nodestate</font>"; }

    //  Get number jobs queued
    $tmpfile = '/tmp/tempfile.topmed';
    $partition = "${host}-incoming";
    $cmd = "/usr/cluster/bin/squeue --format '%A %.14j %.3t %.10M %R' -p $partition > $tmpfile";
    system($cmd);
    $cmd = "wc -l $tmpfile";
    $results = `$cmd`;
    $queuesize = 'empty';
    if (preg_match('/(\d+) /', $results, $m)) { $queuesize = $m[1] - 1; }
    print "<p>Partition <b>$partition</b> has $queuesize jobs queued &nbsp;&nbsp;&nbsp;&nbsp;  $nodestate\n";

    //  Show summary of topmed-related jobs of interest
    if ($queuesize) {
        reset($TOPMEDJOBNAMES);
        foreach ($TOPMEDJOBNAMES as $t) {
            $cmd = "grep $t $tmpfile | wc -l";
            if (preg_match('/(\d+)/', `$cmd`, $m)) { $queued = $m[1]; }
            else { $queued = '?'; }
            if ($queued == 0) { continue; }
            $cmd = "grep $t $tmpfile | grep ' R ' | wc -l";
            if (preg_match('/(\d+)/', `$cmd`, $m)) { $r = $m[1]; }
            else { $r = '?'; }
            print "<br/>&nbsp;&nbsp;&nbsp; $t jobs: $queued queued ($r running)\n";
        }
    
        //  Show a few queued jobs
        $cmd = "tail -$JOBSTOSHOW $tmpfile";
        print "<br/>Last $JOBSTOSHOW queued jobs are:</p><pre>";
        print `$cmd`;
        print "</pre>\n";
  }
  unlink($tmpfile);
}

/*---------------------------------------------------------------
#   html = ViewBamDetail($bamid)
#   Show details for a particular run
---------------------------------------------------------------*/
function ViewBamDetail($bamid) {
    global $LDB, $JS, $conv2dat;

    //  Show all we have for this directory
    $sql = 'SELECT * FROM ' . $LDB['bamfiles'] . " WHERE bamid=$bamid";
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);
    if ($numrows != 1) {
        Emsg("Unable to find detail for BAM '$bamid'");
        return '';
    }
    $row = SQL_Fetch($result);

    //  Now generate the details section for the file
    $html = "<h3 align='center'>Details for '" . $row['bamname'] . "'</h3>\n" .
        "<div class='indent'><table width='80%'>";
    $collist = array_keys($row);
    $bamid = '?';
    foreach ($collist as $c) {
        if ($c == 'bamid') { $bamid = $row[$c]; continue; }
        $d = $row[$c];
        if ($c == 'dateinit') { $d = date('Y/m/d H:i', $d); }
        if ($c == 'datecomplete' && $d != '&nbsp;') { $d = date('Y/m/d H:i', $d); }
        if (in_array($c, $conv2dat)) {  // Special field needs formatting
            $val = DateState($d);
            if ($val != 'notset') {
                $d = "<span class='$val'>$val</span>";

                $url = $_SERVER['PHP_SELF'] . "?bamid=$bamid&amp;col=$c&amp;fcn=setrequest";
                $d .= $JS['SPACER'] . 
                    "<span class='requestcancel'>" .
                    "<a href='$url' onclick='javascript:confirmpopup(\"Request $c\"," .
                    "\"$url\", 640,480); return false;'>Request</a></span>" .
                    " / ";
                $url = $_SERVER['PHP_SELF'] . "?bamid=$bamid&amp;col=$c&amp;fcn=setcancel";
                $d .= "<span class='requestcancel'>" .
                    "<a href='$url' onclick='javascript:confirmpopup(\"Cancel $c\"," .
                    "\"$url\", 640,480); return false;'>Cancel</a></span>";
            }
        }
        if ((! isset($d)) || (is_null($d)) || ($d == '') || ($d == 'NULL')) { $d = 'n/a'; }
        //  If the data has a newline in it do more than just show it
        $p = strpos($d, "\n");
        if ($p) { $d = "<pre>$d</pre>"; }
        $html .= "<tr><td align='left' valign='top' width='15%'><b>" . ucfirst($c) . ":</b></td>" .
            "<td>&nbsp;&nbsp;&nbsp;</td>";
        $html .= "<td align='left' colspan='4'>$d</td>";
        $html .= "</tr>\n";
    }
    $html .= "</table>\n" .
        "<p align='right'><font size='-1'><a href='javascript:window.close()'>Close</a>&nbsp;&nbsp;&nbsp;</p>\n" .
        "</div>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = ControlJobs($h)
#   Show details for a particular run
---------------------------------------------------------------*/
function ControlJobs($h) {
    global $LDB, $SHOWQUEUES, $CENTERS, $CENTERNAME2ID, $validfunctions;
    $url = $_SERVER['PHP_SELF'] . "?fcn=permit";    
    $html = '';

    //  Dump table for each permission for now
    $html .= "<h3 align='center'>Control Job Submission</h3>\n";
    $html .= "<center>These controls do not affect jobs already submitted<br>" .
        "$SHOWQUEUES</center><br/>\n";
    if ($h) { $html .= "<div class='indent'><span class='surprise'>$h</span></div>\n"; }

    //  Prompt for fields which can control permissions
    $html .= "<h4>Use this to stop future job submissions</h4>\n" .
        "<form action='$url' method='post'>\n" .
        "<input type='hidden' name='fcn' value='permit'>\n" .
        "<table align='left' width='80%' border='0'>\n" .
        "<tr>" .
        "<td align='right'><b>Center</b></td>" .
        "<td>&nbsp;</td>" .
        "<td><input type='text' name='center' length='8' value='all'></td>" .
        "<td>&nbsp;</td>" .
        "<td><font color='green'> broad, nygc, uw etc.  </font></td>" .
        "</tr>\n" .
        "<tr>" .
        "<td align='right'><b>Name of Run</b></td>" .
        "<td>&nbsp;</td>" .
        "<td><input type='text' name='run' length='16' value='all'></td>" .
        "<td>&nbsp;</td>" .
        "<td><font color='green'> 2015aug22.weiss.06, 2015oct14 etc. </font></td>" .
        "</tr>\n" .
        "<tr>" .
        "<td align='right'><b>Operation</b></td>" .
        "<td>&nbsp;</td>";
    $html .= "<td><select name='op'>";
    foreach ($validfunctions as $fcn) { $html .= "<option value='$fcn'>$fcn</option>"; }
    $html .= "</select></td>" .
        "<td><font color='green'>&nbsp; </font></td>" .
        "<td>&nbsp;</td>" .
        "</tr>\n" .
        "<tr>" .
        "<td><input type='submit' value=' Stop Submissions '></td>" .
        "<td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td>" .
        "</tr>\n" .
        "</table>\n</form>\n";
    
    //  This shows what is in effect and allows on to delete an entry
    $hdrcols  = array('operation', 'centername', 'dirname', 'centerid', 'runid');
    $sql = 'SELECT * FROM ' . $LDB['permissions'];
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);

    $html .= "<br><br><br><br><br><br><br><br>\n" .
        "<h4>These job submission controls are in effect</h4>\n" .
        "<p>Use '<b>Remove</b>' to remove a control and allow job submissions to resume.</p>\n" .
        "<table align='center' width='90%' border='1'>\n" .
        "<tr>\n";
    foreach ($hdrcols as $c) {
        $html .= "<th class='heading'>" . ucfirst($c) . "</th>\n";
    }
    $html .= "<th class='heading'>&nbsp;</th></tr>\n" .
        "<tr>";

    for ($i=0; $i<$numrows; $i++) {
        $row = SQL_Fetch($result);
        reset($hdrcols);
        foreach ($hdrcols as $c) {
            $d = $row[$c];
            if ((! isset($d)) || ($d == '')) { $d = 'all'; }
            $html .= "<td align='center'>$d</td>\n";
        }
        $html .= "<td align='center'><a href='$url&amp;op=del&amp;id=$row[id]'>" .
            "<font color='green'>Remove</font></a></td>" .
            "</tr>\n";
    }
    $html .= "</table>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = HandlePermit($op, $center, $run, $id)
#   Update permissions table, return HTML about results
---------------------------------------------------------------*/
function HandlePermit($op, $center, $run, $id) {
    global $validfunctions, $LDB;
    $html = "Nothing honey";
    if ($op == '') { return ''; }       // Avoid misleading err msgs

    //  Delete a permission
    if ($op == 'del') {
        $cmd = escapeshellcmd($LDB['bindir'] . "/topmedcmd.pl permit remove $id");
        $s = `$cmd 2>&1`;
        return "<pre>$s</pre>\n";
        //return "<pre>cmd=$cmd\n$s</pre>\n";
    }
    
    //  Set a permission
    if (! in_array($op, $validfunctions)) {
        return Emsg("Invalid operation '$op' - How'd you do that?", 0);
    }
    $cmd = escapeshellcmd($LDB['bindir'] . "/topmedcmd.pl permit add $op $center $run");
    $s = `$cmd 2>&1`;
    return "<pre>$s</pre>\n";
    //return "<pre>cmd=$cmd\n$s</pre>\n";
}

/*---------------------------------------------------------------
#   GetCenters()
#   Set globals with details about all centers
---------------------------------------------------------------*/
function GetCenters() {
    global $LDB, $CENTERID2NAME, $CENTERNAME2ID, $CENTERS;

    $sql = 'SELECT * FROM ' . $LDB['centers'] . ' ORDER BY centername ASC';
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);
    $CENTERID2NAME = array();               // Hash of centerid to names
    $CENTERNAME2ID = array();               // Hash of names to centerid
    $CENTERS = array();                     // Array of names
    for ($j=0; $j<$numrows; $j++) {
        $row = SQL_Fetch($result);
        $i = $row['centerid'];
        $n = $row['centername'];
        $CENTERID2NAME[$i] = $n;
        $CENTERNAME2ID[$n] = $i;
        array_push ($CENTERS,$n);
    }
}

/*---------------------------------------------------------------
# href = QuickStatus($r)
#   $r = row of data for this BAM, contains dates for tasks
#       datearrived,datemd5ver,datebackup,datecp2ncbi
#   Generate a short summary of the state for this directory
#   Parameter is the row from the database
---------------------------------------------------------------*/
function QuickStatus($r) {
    global $quickcols, $quickletter;
    $h = '';
    $span='notset';
    $cols = array_keys($quickcols);
    foreach ($cols as $c) {
        $val = DateState($r[$c]);
        $h .= "<span class='$val'>&nbsp;" . $quickletter[$c] . "&nbsp;</span>";
    }
    return $h;
}

/*---------------------------------------------------------------
# href = DateState($t)
#   Return state for a particular time. See values at top of this pgm
#   returns array (state)
---------------------------------------------------------------*/
function DateState($t) {
    global $state2str;
    $state = 'notset';
    if (in_array($t, $state2str)) { $state = $state2str[$t]; }
    $state = $state2str[$t];
    return $state;
}

/*---------------------------------------------------------------
# html = CalcRunStatus($str)
#   Convert status string into short hand status
#   str should look like: A=done,5=processing,B=unknown etc
#   returns string of html
---------------------------------------------------------------*/
function CalcRunStatus($str) {
    global $quickletter;
    //return $str;          // To see original state
    $h = '';
    $cols = array_values($quickletter);
    foreach ($cols as $c) {
        $span='notset';
        if (preg_match("/$c=([^,]+),/", $str, $m)) { $span = $m[1]; }
        $h .= "<span class='$span'> $c </span>";
    }
    return $h;
}

/*---------------------------------------------------------------
# string = ShortForm(number)
#   Convert a number into a shorter form in K, M or G
#   returns string
---------------------------------------------------------------*/
function ShortForm($n) {
    if (! isset($n))  { return ''; }
    if (! is_numeric($n)) { return $n; }
    if ($n < 102400)  { return $n; }
    if ($n < 1048576) { return sprintf('%.2fK', $n/1024); }
    if ($n < 1073741824) { return sprintf('%.2fM', $n/1048576); }
    return sprintf('%4.2fG', $n/1073741824);
}

/*---------------------------------------------------------------
# href = ShowTable($table, $sql, $tfoot)
#   Generate an HTML dump of a database table
#   returns HTML
---------------------------------------------------------------*/
function ShowTable($table, $sql='', $tfoot=0) {

    $h = "<font size='-1'> <table align='center' border='1' width='90%'>\n";
    //  Builder header for columns in table
    $desc = SQL_Desc($table);
    if (count($desc) <= 0) { return Emsg("Unable to figure out columns in '$table'", 1); }
    $colarray = array();
    foreach ($desc as $c => $val) {
        if (substr($c,0,1) == '_') { continue; }
        array_push($colarray, $c);
    }
    $th = '';
    foreach ($colarray as $index => $c) {
        if (substr($c,0,1) == '_') { continue; }
        $th .= "<th>" . ucfirst($c) . "</th>";
    }
    $h .= "<thead><tr>$th</tr>\n</thead>\n";
    if ($tfoot) { $h .= "<tfoot><tr>$th</tr>\n</tfoot>\n"; }
    $h .= "<tbody>\n";
    
    //  Now generate all the rows in this SQL
    if (! $sql) { $sql = 'SELECT * FROM ' . $table; }
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);
    if ($numrows <= 0) { return Emsg("No data returned from SQL: $sql", 1); }
    for ($i=0; $i<$numrows; $i++) {
        $row = SQL_Fetch($result);
        $h .= "<tr align='center'>";
        foreach ($colarray as $index => $c) {
            $h .= "<td>" . $row[$c] . "</td>";
        }
        $h .= "</tr>\n";
    }
    $h .= "</tbody>\n</table></font>\n";
    return $h;
}


?>
