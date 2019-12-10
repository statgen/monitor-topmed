<?php
/*#################################################################
#
# Name: rna_seq.php
#
# Description:
#   Code to support RNA sequence processing
#
# Copyright (C) 2019 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/

/*---------------------------------------------------------------
#   RNAFunctions()    Most useful data is in $PARMS
#   Handle functions related to RNA processing. Might not return
---------------------------------------------------------------*/
function RNAFunctions($fcn) {
    global $HDR, $LDB, $GLOBS, $PARMS;

	if ($fcn == 'projects') {
		print $infotext;
		print ViewProjects($PARMS['center'], $PARMS['maxdir']);        
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

	if ($fcn == 'files') {
		print ViewFiles($PARMS['id']);
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
    $html = "<h3 align='center'>RNA Sequence Data Projects</h3>\n" .
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
                $u = $_SERVER['SCRIPT_NAME'] . "?fcn=samples&amp;id=" . $row[$projpkey];
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
    $hdrcols  = array('expt_sampleid', 'QUIKSTAT', 'count', 'fileprefix', 'rnasubject');

	$samplestable = 'samples';
    $samplespkeynick = $samplestable . '_pkey';
 	$samplestable = $LDB[$samplestable];
    $samplespkey = $LDB[$samplespkeynick];
	$projectstable = 'projects';
    $projectspkeynick = $projectstable . '_pkey';
 	$projectstable = $LDB[$projectstable];
    $projectspkey = $LDB[$projectspkeynick];
    //  Get parent for this sample
    $sql = "SELECT $projectspkey FROM $samplestable WHERE $samplespkey=$id";
    $result = SQL_Query($sql, 0);
    $row = SQL_Fetch($result);
    $projectid = $row[$projectspkey];
    //  From parent table get columns of interest and id of center
    $sql = "SELECT centerid,dirname,count FROM $projectstable WHERE $projectspkey=$projectid";
    $result = SQL_Query($sql, 0);
    $row = SQL_Fetch($result);
    $centerid = $row['centerid'];
    $dirname = $row['dirname'];
    $count = $row['count'];
    $sql = "SELECT centername FROM " . $LDB['centers'] . " WHERE " . $LDB['centers_pkey'] . "=$centerid";
    $result = SQL_Query($sql, 0);
    $row = SQL_Fetch($result);
    $center = $row['centername'];

    $sql = "SELECT * FROM $samplestable WHERE $projectspkey=$id";
    $url = $_SERVER['SCRIPT_NAME'] . "?center=$center&amp;maxdir=$maxdir";
    $hdr = "<h3 align='center'>$count Samples for '$dirname' in center " .
        "<a href='$url'>$center</a></h3>\n";

	return ShowSamples($sql, $hdrcols, 'samples', $hdr);
}

/*---------------------------------------------------------------
#   html = ViewFiles()
#   Show list of files for a particular RNA sample
---------------------------------------------------------------*/
function ViewFiles($id, $sample) {
    global $LDB, $GLOBS, $PARMS;
    $hdrcols  = array('intar', 'checksum', 'filename' );
    $filestable = 'files';
    $filespkey = $filestable . '_pkey';
    $filestable = $LDB[$filestable];
    $filespkey = $LDB[$filespkey];
    $samplestable = 'samples';
    $samplespkey = $samplestable . '_pkey';
    $samplespkey = $LDB[$samplespkey];
    $samplestable = $LDB[$samplestable];
    $center = $PARMS['center'];
    $html = '';
 
    //  Generate HTML header for page
    $sql = "SELECT expt_sampleid FROM $samplestable WHERE $samplespkey=$id";
    $result = SQL_Query($sql);
    $row = SQL_Fetch($result);
    $sample = $row['expt_sampleid'];

    $sql = "SELECT * FROM $filestable WHERE $samplespkey=$id ORDER BY $filespkey";
    $url = $_SERVER['SCRIPT_NAME'] . "?center=$center&amp;maxdir=$maxdir";
    $hdr = "<h3 align='center'>Files Associated with Sample '$sample' [$id] in center " .
        "<a href='$url'>$center</a></h3>\n";

	return ShowSamples($sql, $hdrcols, 'samples', $hdr, 0);
}

?>
