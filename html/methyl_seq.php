<?php
/*#################################################################
#
# Name: methyl_seq.php
#
# Description:
#   Code to support methylation processing
#
# Copyright (C) 2019 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/

/*---------------------------------------------------------------
#   MethylFunctions()    Most useful data is in $PARMS
#   Handle functions related to methylation processing. Might not return
---------------------------------------------------------------*/
function MethylFunctions($fcn) {
    global $HDR, $LDB, $GLOBS, $PARMS;

	if ($fcn == 'projects') {
		print $infotext;
		print ViewProjects($PARMS['center'], $PARMS['maxdir']);        
		print "<center>" . GetChooseLines() . "<br/>" . $GLOBS['links'] . "</center>\n";
		print dofooter($HDR['footer']);
		exit;
	}
	if ($fcn == 'batches') {
		print ViewBatches($PARMS['id'], $PARMS['maxdir']);
		print "<center>" . GetChooseLines() . "<br/>" . $GLOBS['links'] . "</center>\n";
		print dofooter($HDR['footer']);
		exit;
	}
	if ($fcn == 'samples') {
		print ViewSamples($PARMS['id'], $PARMS['maxdir']);
		print "<center>" . GetChooseLines() . "<br/>" . $GLOBS['links'] . "</center>\n";
		print dofooter($HDR['footer']);
		exit;
	}
	if ($fcn == 'detail') {
		$tablenick = $PARMS['table'];
		$id = $PARMS['id'];
		$t = $LDB[$tablenick];
		$pkey = $tablenick . '_pkey';
		$pkey = $LDB[$pkey];
		$HTML .= "<h3 align='center'>View Database Entry [table=$t entry=$id]</h3>\n" .
			"<div class='indent'><table width='80%'>\n";
		$sql = "SELECT * FROM $t WHERE $pkey='$id'";
    	$HTML .= ViewSampleDetail($tablenick, $id, $sql);
		print $HTML . "</div>\n";
		print dofooter($HDR['footer']);
		exit;
	}
}

/*---------------------------------------------------------------
#   html = ViewProjects()
#   Show summary of directories of runs for all datayears
---------------------------------------------------------------*/
function ViewProjects($center, $maxdir) {
    global $GLOBS;
    $hdrcols  = array('dirname', 'status', 'count');	// Columns for this type of data
    //  Generate HTML header for page
    //  Get list of centers doing:  select distinct(project) from status;
    $html = "<h3 align='center'>Methylation Data Projects</h3>\n" .
        "<center>" . GetChooseLines() . "<br/>" . $GLOBS['links'] . "</center>\n";
    $yearstart = 2019;
    $yearstop = 2018;
    $centers2show = array();                // Get list of centers for this query
    if ($center) { array_push($centers2show, $center); }
    else { $centers2show = $GLOBS['centers']; }

    //  Show data for center by datayear
    for ($datayear=$yearstart; $datayear>$yearstop; $datayear--) {
        //  For each center show details from database ($rows)
        foreach ($centers2show as $centr) {
            $cid = $GLOBS['centername2id'][$centr];
            $h = ShowProjectYear($cid, $maxdir, $datayear, $GLOBS[iammgr]);
            if (! $h) { continue; }
            $html .= "<br><div class='indent'><b>" . strtoupper($centr) .
                ", Year $datayear</b></div>$h\n";
        }
        $html .= "<br>";
        if ($datayear > 1) { $html .= "<hr width='80%' class='separator'>\n"; }
    }
    
    $html .= "<div class='indent'>\n" . $GLOBS['statusruns'] . "</div>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = ShowProjectYear() {
#   Show summary of directories of runs
---------------------------------------------------------------*/
function ShowProjectYear($cid, $maxdir, $datayear) {
    global $LDB, $GLOBS;
    $hdrcols  = array('dirname', 'status', 'count');
    $tablenick = 'projects';
    $projtable = $tablenick;
    $projpkey = $projtable . '_pkey';
    $projtable = $LDB[$projtable];
    $projpkey = $LDB[$projpkey];

    //  Walk through database getting data for this center
    $sql = "SELECT * FROM $projtable WHERE centerid=$cid AND datayear=$datayear " .
    	"ORDER BY $projpkey DESC";
    if ($maxdir) { $sql .= " LIMIT $maxdir"; }
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);
    $rows = array();            // Save DB info for later display
    $centers2show = array();    // Get list of projects for this query
    for ($i=0; $i<$numrows; $i++) {
        $row = SQL_Fetch($result);
        $rows[$row[$projpkey]] = $row;
    }
    if (! count($rows)) { return ''; }        // Nothing here

    //  Build start of table for each center
    $html = "<table align='center' width='100%' border='1'><tr>\n";
    foreach ($hdrcols as $c) {
        $html .= "<th class='heading'>" . ucfirst($c) . "</th>\n";
    }
    if ($GLOBS['iammgr']) { $html .= "<th>&nbsp;</th>"; }
    $html .= "</tr>\n";

    reset($rows);
    foreach ($rows as $id => $row) {
        //  Show data for this project
        $html .= "<tr>\n";
        reset($hdrcols);
        foreach ($hdrcols as $c) {
            $d = $row[$c];
            if ((! isset($d)) || ($d == '')) { $d = '&nbsp;'; }
            if ($c == 'dirname') {
                $u = $_SERVER['SCRIPT_NAME'] . "?fcn=batches&amp;id=" . $row[$projpkey];
                $d = "<a href='$u'>$d</a>";
            }
            if ($c == 'status') { $d = CalcRunStatus($d); }
            $html .= "<td align='center'>$d</td>\n";
        }
            
        $html .= "<td align='center'>";
        $u = $_SERVER['SCRIPT_NAME'] ."?fcn=detail&amp;table=$tablenick&amp;id=" . $row[$projpkey];
        $html .= "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
            "<font color='green' size='-2'>Details</font></a>&nbsp;";
        if ($GLOBS['iammgr']) {
            $u = $_SERVER['SCRIPT_NAME'] ."?fcn=edit&amp;table=$tablenick&amp;id=" . $row[$projpkey];
            $html .= "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
                "<font color='red' size='-2'>Edit</font></a>";
        }
        $html .= "</td></tr>\n";
    }
    $html .= "</table>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = ViewSamples()
#   Show list of RNA files for a particular project
---------------------------------------------------------------*/
function ViewSamples($id, $maxdir) {
    global $LDB, $GLOBS, $PARMS;
    $hdrcols  = array('expt_sampleid', 'QUIKSTAT', 'type', 'version', 'rowcol');

	$samplestable = 'samples';
    $samplespkeynick = $samplestable . '_pkey';
 	$samplestable = $LDB[$samplestable];
    $samplespkey = $LDB[$samplespkeynick];

	$batchestable = 'batches';
    $batchespkeynick = $batchestable . '_pkey';
 	$batchestable = $LDB[$batchestable];
    $batchespkey = $LDB[$batchespkeynick];
	$projectstable = 'projects';
    $projectspkeynick = $projectstable . '_pkey';
 	$projectstable = $LDB[$projectstable];
    $projectspkey = $LDB[$projectspkeynick];


    //  From parent table get columns of interest
    $sql = "SELECT $projectspkey,batchname,count FROM $batchestable WHERE $batchespkey=$id";
    $result = SQL_Query($sql, 0);
    $row = SQL_Fetch($result);
    $projectspkey = $row[$projectspkey];
    $batchname = $row['batchname'];
    $count = $row['count'];
    //  From grandparent table get center
    $sql = "SELECT centerid FROM $projectstable WHERE $projectspkey=$projectspkey";
    $result = SQL_Query($sql, 0);
    $row = SQL_Fetch($result);
    $centerid = $row['centerid'];
    $sql = "SELECT centername FROM " . $LDB['centers'] . " WHERE " . $LDB['centers_pkey'] . "=$centerid";
    $result = SQL_Query($sql, 0);
    $row = SQL_Fetch($result);
    $centername = $row['centername'];

    $sql = "SELECT * FROM $samplestable WHERE $batchespkey=$id";
    $url = $_SERVER['SCRIPT_NAME'] . "?center=$center&amp;maxdir=$maxdir";
    $hdr = "<h3 align='center'>$count Samples for '$dirname' in center " .
        "<a href='$url'>$center</a></h3>\n";

	return ShowSamples($sql, $hdrcols, 'samples', $hdr);
}

/*---------------------------------------------------------------
#   html = ViewBatches()
#   Show list of batches for a particular methylation project
---------------------------------------------------------------*/
function ViewBatches($id, $sample) {
    global $LDB, $GLOBS, $PARMS;
    $hdrcols  = array('batchname', 'count', 'file0' );
    $projtable = 'projects';
    $projpkey = $projtable . '_pkey';
    $projtable = $LDB[$projtable];
    $projpkey = $LDB[$projpkey];
    $batchestable = 'batches';
    $tablenick = $batchestable;
    $batchespkey = $batchestable . '_pkey';
    $batchestable = $LDB[$batchestable];
    $batchespkey = $LDB[$batchespkey];
    $samplestable = 'samples';
    $samplespkey = $samplestable . '_pkey';
    $samplespkey = $LDB[$samplespkey];
    $center = $PARMS['center'];
    $html = '';
 
    //  Generate HTML header for page
    $sql = "SELECT dirname FROM $projtable WHERE $projpkey=$id";
    $result = SQL_Query($sql);
    $row = SQL_Fetch($result);
    $proj = $row['dirname'];
    $sql = "SELECT * FROM $batchestable WHERE $projpkey=$id ORDER BY $batchespkey";
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);

    //  Show details for each file
    $url = $_SERVER['SCRIPT_NAME'] . "?center=$centername&amp;maxdir=" . $PARMS['maxdir'];
    $html .= "<h3 align='center'>Batches in project '$proj' [$id] in center " .
        "<a href='$url'>$center</a></h3>\n";
    $html .= "<p align='center'/>" . $GLOBS['links'] . "</p>\n";

    $html .= "<table align='center' width='100%' border='1'><tr>\n";
    foreach ($hdrcols as $c) {
        $html .= "<th class='heading'>" . ucfirst($c) . "</th>\n";
    }
    if ($GLOBS['iammgr']) { $html .= "<th>&nbsp;</th>"; }
    $html .= "</tr>\n";

    for ($i=0; $i<$numrows; $i++) {
        $row = SQL_Fetch($result);
        reset($hdrcols);
        foreach ($hdrcols as $c) {
            $d = $row[$c];
        	if ($c == 'batchname') {
                $u = $_SERVER['SCRIPT_NAME'] . "?fcn=samples&amp;id=" . $row[$batchespkey];
                $d = "<a href='$u'>$d</a>";
            }
            if ((! isset($d)) || ($d == '')) { $d = '&nbsp;'; }
            $html .= "<td align='center'>$d</td>\n";
        }
        $html .= "<td align='center'>";
        $u = $_SERVER['SCRIPT_NAME'] ."?fcn=detail&amp;table=$tablenick&amp;id=" . $row[$batchespkey];
        $html .= "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
            "<font color='green' size='-2'>Details</font></a>&nbsp;";
        if ($GLOBS['iammgr']) {
            $u = $_SERVER['SCRIPT_NAME'] ."?fcn=edit&amp;table=$tablenick&amp;id=" . $row[$batchespkey];
            $html .= "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
                "<font color='red' size='-2'>Edit</font></a>";
        }
        $html .= "</td></tr>\n";
    }
    $html .= "</table>\n";
    return $html;
}

?>
