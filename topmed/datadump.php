<?php
/*#################################################################
#
# Name: datadump.php
#
# Description:
#   Dump data from the database in a convenient format
#
# Copyright (C) 2016 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/
$MON='topmed'; include_once 'local_config.php';
include_once 'common.php';
include_once 'DBMySQL.php';
include_once "download.php";

//print "<!-- _POST=\n"; print_r($_POST); print " -->\n";
$NOTSET = 0;                // Not set
$REQUESTED = 1;             // Task requested
$SUBMITTED = 2;             // Task submitted to be run
$STARTED   = 3;             // Task started
$DELIVERED = 19;            // Data delivered, but not confirmed
$COMPLETED = 20;            // Task completed successfully
$CANCELLED = 89;            // Task cancelled
$FAILEDCHECKSUM = 98;       // Task failed, because checksum at NCBI bad
$FAILED    = 99;            // Task failed
$state2str = array(         // Values here are class for SPAN tag
    $NOTSET => 'notset',
    $REQUESTED => 'requested',
    $SUBMITTED => 'submitted',
    $STARTED => 'started',
    $DELIVERED => 'delivered',
    $COMPLETED => 'completed',
    $CANCELLED => 'cancelled',
    $FAILEDCHECKSUM => 'failedchecksum',
    $FAILED => 'failed'
);

//-------------------------------------------------------------------
//  You could also do this in one giant select like this
//
//  select c.centername,r.dirname,b.expt_sampleid,b.datayear from bamfiles as b join runs as r on  b.runid = r.runid  join centers as c on c.centerid = r.centerid;

//-------------------------------------------------------------------

DB_Connect($LDB['realm']);
GetCenters();                   // Get maps to identify centers

$dumppart1 = "center\trun\trunid\t";
$dumpcols = array("expt_sampleid",
    "datayear",
    "piname",
    "state_ncbiorig",
    "state_ncbib37",
    "state_ncbib38",
    "bamid",
    "dateinit",
    "build",
    "bamflagstat",
    "cramflagstat",
    "bamsize"
);
$dumpstring = $dumppart1 . implode("\t",$dumpcols) . "\n";

//  Dump all nwdid, run and center
$sql = 'SELECT * FROM ' . $LDB['runs'] . ' ORDER BY centerid';
$runresult = SQL_Query($sql);
$runnumrows = SQL_NumRows($runresult);
for ($r=0; $r<$runnumrows; $r++) {
    $row = SQL_Fetch($runresult);
    $dirname = $row['dirname'];
    $runid = $row['runid'];
    $cid = $row['centerid'];
    $centername = $CENTERID2NAME[$cid];
    $sql = 'SELECT * FROM ' . $LDB['bamfiles'] . " WHERE runid=$runid";
    $bamresult = SQL_Query($sql);
    $bamnumrows = SQL_NumRows($bamresult);
    for ($b=0; $b<$bamnumrows; $b++) {
        $row = SQL_Fetch($bamresult);
        if ($row['expt_sampleid'] && $row['datayear']) {
            $dumpstring .= $CENTERID2NAME[$cid] .  "\t" . $dirname .  "\t" . $runid;
            foreach ($dumpcols as $c) {
                if (substr($c, 0, 6) == 'state_') { $dumpstring .= "\t" . $state2str{$row[$c]}; }
                else { $dumpstring .= "\t" . $row[$c]; }
            }
            $dumpstring .= "\n";
        }
    }
}
//print "<pre>\n$dumpstring</pre>\n"; exit;
Download($fn='bamfile_datadump.txt', $dumpstring);
exit;


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

?>
