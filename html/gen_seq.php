<?php
/*#################################################################
#
# Name: gen_seq.php
#
# Description:
#   Code to support Genome sequence processing
#
# Copyright (C) 2019 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/

/*---------------------------------------------------------------
#   GENFunctions()
#   Handle functions related to genome processing. Might not return
---------------------------------------------------------------*/
function GENFunctions($fcn) {
    global $HDR, $LDB, $GLOBS, $PARMS;

	if ($fcn == 'runs') {
		if ($GLOBS['iammgr']) {
    		$s = "See TOPMed monitor docs " .
        		"<a href='https://statgen.sph.umich.edu/wiki/NHLBI_automation_steps' target='_blank'>here</a>.\n";
		}
		else { $s = ''; }
		print "<p class='intro'>The <a href='http://www.nhlbi.nih.gov/'>NHLBI</a> provides " .
    		"science-based, plain-language information related to heart, lung " .
    		"and blood diseases and conditions and sleep disorders.\n" .
    		"Details about this data tracking are available from Tom Blackwell " .
    		"(<a href='mailto:tblackw@umich.edu'>tblackw@umich.edu</a>). $s</p>\n";

		print ViewRuns($PARMS['center'], $GLOBS['maxdir']);        
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
	if ($fcn == 'bams') {
		print ViewBams($PARMS['id'], $GLOBS['maxdir']);
		print "<center>" . GetChooseLines() . "<br/>" . $GLOBS['links'] . "</center>\n";
		print dofooter($HDR['footer']);
		exit;
	}
}

/*---------------------------------------------------------------
#   html = ViewRuns()
#   Show summary of directories of runs for all datayears
---------------------------------------------------------------*/
function ViewRuns($center, $maxdir) {
    global $LDB, $GLOBS, $PARMS;

    //  Generate HTML header for page
    //  Get list of centers doing:  select distinct(project) from status;
    $html = "<h3 align='center'>Summary of Genome Sequence Data</h3>\n" .
        "<center>" . GetChooseLines() . "<br/>" . $GLOBS['links'] . "</center>\n";
    $yearstart = 5;
    $yearstop = 1;
    $centers2show = array();                // Get list of centers for this query
    //	For some reason we're using center to spcify year. That's dumb
    if ($center) {
		if (preg_match('/year(\d)/',$center,$m)) {
			$centers2show = $GLOBS['centers'];
			$yearstart = $m[1];
			$yearstop = $m[1];
			if ($m[1] > '1') { $yearstop = $m[1] - 1; }
		}
		else { array_push($centers2show, $center); }
    }
    else { $centers2show = $GLOBS['centers']; }

    //  Show data for center by datayear
    $lineshown = 0;
    for ($datayear=$yearstart; $datayear>=$yearstop; $datayear--) {
        //  For each center show details from database ($rows)
        foreach ($centers2show as $centr) {
            $cid = $GLOBS['centername2id'][$centr];
            $h = ShowRunYear($cid, $maxdir, $datayear, $GLOBS[iammgr]);
            if (! $h) { continue; }
            $html .= "<br><div class='indent'><b>" . strtoupper($centr) .
                ", Year $datayear</b></div>$h\n";
            $lineshown++;
        }
        if ($lineshown) {
        	$html .= "<br>";
        	if ($datayear > 1) { $html .= "<hr width='80%' class='separator'>\n"; }
        }
    }
    $html .= "<div class='indent'>\n" . $GLOBS['statusruns'] . "</div>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = ShowRunYear() {
#   Show summary of directories of runs
---------------------------------------------------------------*/
function ShowRunYear($cid, $maxdirs, $datayear, $iammgr) {
    global $LDB, $GLOBS, $PARMS;
    $hdrcols  = array('dirname', 'status', 'count', 'build');

    $tablenick = 'runs';
    $projtable = $tablenick;
    $projpkey = $projtable . '_pkey';
    $projtable = $LDB[$projtable];
    $projpkey = $LDB[$projpkey];
    $samplesnick = 'samples';
    $samplestable = $samplesnick;
	$samplespkey = $samplestable . '_pkey';
	$samplestable = $LDB[$samplestable];
    $samplespkey = $LDB[$samplespkey];

    //  Walk through database getting data for this center
    $sql = "SELECT * FROM $projtable WHERE centerid=$cid " .
    	"AND datayear=$datayear ORDER BY $projpkey DESC";
    if ($maxdirs) { $sql .= " LIMIT $maxdirs"; }
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);
    $rows = array();            // Save DB info for later display
    $centers2show = array();    // Get list of runs for this query
    for ($i=0; $i<$numrows; $i++) {
        $row = SQL_Fetch($result);
        //  Get build for run. This should be done much smarter in SQL
        $buildsql = "SELECT build FROM $samplestable WHERE $projpkey=" . $row[$projpkey] . " LIMIT 1";
        $buildresult = SQL_Query($buildsql);
        $buildrow = SQL_Fetch($buildresult);
        $row['build'] = $buildrow['build'];
        $rows[$row[$projpkey]] = $row;
    }
    if (! count($rows)) { return ''; }        // Nothing here

    //  Build start of table for each center
    $html = "<table align='center' width='100%' border='1'><tr>\n";
    foreach ($hdrcols as $c) {
        $html .= "<th class='heading'>" . ucfirst($c) . "</th>\n";
    }
    if ($iammgr) { $html .= "<th>&nbsp;</th>"; }
    $html .= "</tr>\n";

    reset($rows);
    foreach ($rows as $id => $row) {
        //  Show data for this run
        $html .= "<tr>\n";
        reset($hdrcols);
        foreach ($hdrcols as $c) {
            $d = $row[$c];
            if ((! isset($d)) || ($d == '')) { $d = '&nbsp;'; }
            if ($c == 'dirname') {
                $u = $_SERVER['SCRIPT_NAME'] . "?fcn=bams&amp;id=" . $row['runid'];
                $d = "<a href='$u'>$d</a>";
            }
            if ($c == 'status') { $d = CalcRunStatus($d); }
            $html .= "<td align='center'>$d</td>\n";
        }
            
        $html .= "<td align='center'>";
        $u = $_SERVER['SCRIPT_NAME'] ."?fcn=detail&amp;table=$tablenick&amp;id=" . $row[$projpkey];
        $html .= "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
            "<font color='green' size='-2'>Details</font></a>&nbsp;";
        if ($iammgr) {
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
#   html = ViewBams($runid, $maxdirs)
#   Show list of BAM files for a particular run
---------------------------------------------------------------*/
function ViewBams($id, $maxdir) {
    global $LDB, $GLOBS, $PARMS;
    $hdrcols  = array('bamname', 'QUIKSTAT', 'bamsize', 'piname');

	$samplestable = 'samples';
    $samplespkeynick = $samplestable . '_pkey';
 	$samplestable = $LDB[$samplestable];
    $samplespkey = $LDB[$samplespkeynick];
	$runstable = 'runs';
    $runspkeynick = $runstable . '_pkey';
 	$runstable = $LDB[$runstable];
    $runspkey = $LDB[$runspkeynick];

    //  Get columns of interest and id of center
    $sql = "SELECT centerid,dirname,count FROM $runstable WHERE $runspkey=$id";
    $result = SQL_Query($sql, 0);
    $row = SQL_Fetch($result);
    $centerid = $row['centerid'];
    $dirname = $row['dirname'];
    $count = $row['count'];
    $sql = "SELECT centername FROM " . $LDB['centers'] . " WHERE " . $LDB['centers_pkey'] . "=$centerid";
    $result = SQL_Query($sql, 0);
    $row = SQL_Fetch($result);
    $center = $row['centername'];

    $sql = "SELECT * FROM $samplestable WHERE $runspkey=$id";
    $url = $_SERVER['SCRIPT_NAME'] . "?center=$center&amp;maxdir=$maxdir";
    $hdr = "<h3 align='center'>$count Samples for '$dirname' [$id] in center " .
        "<a href='$url'>$center</a></h3>\n";

	return ShowSamples($sql, $hdrcols, 'samples', $hdr);
}

?>
