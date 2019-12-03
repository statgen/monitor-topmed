<?php
/*#################################################################R
#
# Name: rnaseq.php
#
# Description:
#   Display information for tracking NHLBI TopMed RNA data
#
# Copyright (C) 2015-2019 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/
include_once 'local_config.php';
if ($LDB['datatype'] == 'genome') { include_once 'gen_seq.php'; }
if ($LDB['datatype'] == 'rnaseq') { include_once 'rna_seq.php'; }
include_once 'managejobs.php';
include_once 'iammgrfcns.php';
include_once 'common.php';
include_once 'header.php';
include_once 'DBMySQL.php';
include_once "edit.php";
include_once "quickletters.php";

//print "<!-- _SERVER=\n"; print_r($_SERVER); print " -->\n";
//print "<!-- GLOBS=\n"; print_r($GLOBS); print "\n PARMS=\n"; print_r($PARMS); print "\n LDB=\n"; print_r($LDB); print " -->\n";

$GLOBS['links'] = "LINKS: " .
    "<a onclick='javascript:window.location.reload()'><img src='refresh.png' alt='refresh'></a> &nbsp;" .
    "<a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=queue' " .
    "onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=queue\",680,720); " .
    "return false;'>Local_Queue</a> &nbsp;" .
	"<a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=logs' " .
    "onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=logs\",680,720); " .
    "return false;'>Tail_Logs</a> &nbsp;";

if ($LDB['datatype'] == 'rnaseq') {
	$$GLOBS['validfunctions'] = array('all', 'verify', 'backup', 'awscopy');
}
if ($LDB['datatype'] == 'genome') {
	$$GLOBS['validfunctions'] = array('all', 'verify', 'cram', 'backup', 'qplot',
    	'gcepush', 'gcepull', 'bcf', 'gcecopy', 'gcecpbcf', 'gcecleanup', 'awscopy', 'fix');
}

//-------------------------------------------------------------------
//  Get parameters passed in via normal invocation
//-------------------------------------------------------------------
putenv('PROJECT=topmed');           # Set env variable so pgms play nice
if (! isset($_SERVER['REMOTE_USER'])) { $_SERVER['REMOTE_USER'] = 'none'; }

$GLOBS['iammgr'] = 0;
if (in_array($_SERVER['REMOTE_USER'], $MGRS)) { $GLOBS['iammgr'] = 1; }
if (in_array($_SERVER['REMOTE_USER'], $REQMGRS)) { $GLOBS['iammgr'] = 1; }
//  If a manager, allow them to control job submission
if ($GLOBS['iammgr']) {
    $GLOBS['links'] .= "<a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=control' " .
        "onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=control\",680,720); " .
        "return false;'>Control_Jobs</a> &nbsp;" .
        "<a href='" . $_SERVER['SCRIPT_NAME'] . "?fcn=restart' " .
        "onclick='javascript:popup2(\"" . $_SERVER['SCRIPT_NAME'] . "?fcn=restart\",680,720); " .
        "return false;'>Restart_Jobs</a> &nbsp;" .       
        // "<a href='https://www-inpsyght1.sph.umich.edu/inpsyght/' target='_blank'>Inpsyght</a> &nbsp;" .
        //"<a href='/topmed/origindex.php' target='_blank'>Orig_Index</a> &nbsp;" .
        "<a href='/topmed/rnaseq.php' target='_blank'>RNA_Seq</a> &nbsp;" .
        "<a href='/topmed/genome.php' target='_blank'>Genome_Seq</a> &nbsp;";
    if ($LDB['datatype'] == 'genome') {
         $GLOBS['links'] .= "<a href='http://nhlbi.sph.umich.edu/report/monitor.php' " .
        "onclick='javascript:popup2(\"http://nhlbi.sph.umich.edu/report/monitor.php\",680,720); " .
        "return false;'>ReMapping</a> &nbsp;" .
        "<a href='https://console.cloud.google.com/storage/browser/' target='_blank'>Browser</a> &nbsp;";
    }
}

print doheader($HDR['title'], 1);

//  Real parameters for each form, default is ''
//	Removed these Npv 30: 'desc', 'sample', 'pkey', 'runid', 'bamid', 'centerid', 'fetchpath', 'hostname',
$parmcols = array('fcn', 'maxdir', 's', 'col', 'center', 'datayear',   
    'table', 'op', 'id', 'samplestate', 'run');
//	Capture $parmcols in hash $PARMS so we can access them from anywhere easily
$PARMS = isolate_parms($parmcols);
//extract (isolate_parms($parmcols));

//	Set defaults
if ($LDB['datatype'] == 'genome') {
	if ($PARMS['center'] == '') { $PARMS['center'] = 'year4'; }
	if ($PARMS['fcn'] == '')    { $PARMS['fcn'] = 'runs'; }
	if ($PARMS['maxdir'] == '') { $PARMS['maxdir'] = 0; }
}
if ($LDB['datatype'] == 'rnaseq') {
	if ($PARMS['datayear'] == '')   { $PARMS['datayear'] = '2019'; }
	if ($PARMS['center'] == '')     { $PARMS['center'] = 'uw'; }
	if ($PARMS['fcn'] == '' || $PARMS['fcn'] == 'runs')     { $PARMS['fcn'] = 'projects'; }
}

DB_Connect($LDB['realm']);
GetCenters();                       // Get maps to identify centers
$HTML = '';

//	Figure out what we're supposed to do   Useful data is in $PARMS and $GLOBS
$fcn = $PARMS['fcn'];
if ($LDB['datatype'] == 'rnaseq') {
	RNAFunctions($fcn);
}

if ($LDB['datatype'] == 'genome') {
	GENFunctions($fcn);
}
JOBFunctions($fcn);		// Handle JOBs
if ($GLOBS['iammgr']) { ImmgrFunctions($fcn);  }	// Handle functions for managers

//  What was that?
Emsg("Unknown directive '" . $fcn . "'.");
Nice_Exit("How'd you do that?");
exit;

/*---------------------------------------------------------------
#   $html = GetChooseLines()
#   Set globals with details about all centers
---------------------------------------------------------------*/
function GetChooseLines() {
    global $LDB, $GLOBS, $PARMS;
	$maxdir = $PARMS['maxdir'];
	$html = 'Choose: '; 		// Ain't consistency great?
	if ($LDB['datatype'] == 'genome') {
		$html .= "&nbsp;&nbsp;<a href='" . $_SERVER['SCRIPT_NAME'] .
			"?center=year5&amp;maxdir=$maxdir'>Year 5</a>&nbsp;&nbsp;\n" .
			"<a href='" . $_SERVER['SCRIPT_NAME'] .
			"?center=year4&amp;maxdir=$maxdir'>Year 4</a>&nbsp;&nbsp;\n" .
			"<a href='" . $_SERVER['SCRIPT_NAME'] .
			"?center=year3&amp;maxdir=$maxdir'>Year 3</a>&nbsp;&nbsp;\n" .
			"<a href='" . $_SERVER['SCRIPT_NAME'] .
			"?center=year2&amp;maxdir=$maxdir'>Year 2</a>&nbsp;&nbsp;\n" .
			"<a href='" . $_SERVER['SCRIPT_NAME'] .
			"?center=year1&amp;maxdir=$maxdir'>Year 1</a>&nbsp;&nbsp;\n";
	}
	if ($LDB['datatype'] == 'rnaseq') {
		$html .= "&nbsp;&nbsp;<a href='" . $_SERVER['SCRIPT_NAME'] .
        "?datayear=2019&amp;maxdir=$maxdir'>2019</a>&nbsp;&nbsp;\n";
    }
    foreach ($GLOBS['centers'] as $c) {
        $html .= "<a href='" . $_SERVER['SCRIPT_NAME'] .
            "?center=$c&amp;maxdir=$maxdir'>" . strtoupper($c) . "</a>&nbsp;&nbsp;\n";
    }
    return $html;
}

/*---------------------------------------------------------------
#   GetCenters()
#   Set globals with details about all centers
---------------------------------------------------------------*/
function GetCenters() {
    global $LDB, $GLOBS;

	$GLOBS['centers'] = array();
	$GLOBS['centerid2name'] = array();    	// Hash of centerid to names
	$GLOBS['centername2id'] = array();		// Hash of names to centerid
	//	General case, get centers from database
    $sql = 'SELECT * FROM ' . $LDB['centers'];
	$sql .= " WHERE datatype LIKE '%" . $LDB['datatype'] . "%'";
	$sql .= ' ORDER BY centername ASC';
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);
    for ($j=0; $j<$numrows; $j++) {
        $row = SQL_Fetch($result);
        $i = $row['centerid'];
        $n = $row['centername'];
        $GLOBS['centername2id'][$n] = $i;
        $GLOBS['centerid2name'][$i] = $n;
        array_push ($GLOBS['centers'],$n);
    }
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

/*---------------------------------------------------------------
#   html = ShowSamples()
#   Show list of samples files for a particular set of data
---------------------------------------------------------------*/
function ShowSamples($id, $hdrcols, $tablenick, $parenttablenick) {
    global $LDB, $GLOBS, $PARMS;
    $html = '';
    $maxdir = 0;                    // For now, show all samples
	$samplestable = $tablenick;
    $samplespkey = $samplestable . '_pkey';
    $projtable = $parenttablenick;
	$projpkey = $parenttablenick . '_pkey';
    $samplestable = $LDB[$samplestable];
    $samplespkey = $LDB[$samplespkey];
    $projtable = $LDB[$projtable];
    $projpkey = $LDB[$projpkey];

    //  Figure out project for this sample
    $sql = "SELECT dirname,centerid,count FROM $projtable WHERE $projpkey=$id";
    $result = SQL_Query($sql, 0);
    $e = DB_CheckError($result);
    if ($e) { return EMsg("Surprise: Unable to find $projtable $projpkey=$id: $e", 1); }
    $row = SQL_Fetch($result);
    $projectname = $row['dirname'];
    $count = $row['count'];
    $centername = $GLOBS['centerid2name'][$row['centerid']];
	$center = strtoupper($centername);

    //  Get all samples from this project
    $sql = "SELECT * FROM $samplestable WHERE $projpkey=$id";
    if ($maxdir) { $sql .= " LIMIT $maxdir"; }
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);

    //  Show details for each sample
    $url = $_SERVER['SCRIPT_NAME'] . "?center=$centername&amp;maxdir=$maxdir";
    $html .= "<h3 align='center'>$count Samples for '$projectname' [$id] in center " .
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
            if ($c == 'QUIKSTAT') {
                $u = $_SERVER['SCRIPT_NAME'] . "?fcn=showout&amp;table=$tablenick&amp;id=" .
                    $row[$samplespkey] . "&amp;samplestate=XX";
                $u = "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>XX</a>";
                $d = QuickStatus($row, $u);     // Pass on URL for failing states
            }
            else { $d = $row[$c]; }
            if ((! isset($d)) || ($d == '')) { $d = '&nbsp;'; }
            $html .= "<td align='center'>$d</td>\n";
        }           

        $html .= "<td align='center'>";
		if ($LDB['datatype'] == 'rnaseq') {
         	$u = $_SERVER['SCRIPT_NAME'] . "?fcn=files&amp;id=" . $row[$samplespkey];
         	$html .= "<a href='$u' onclick='javascript:popup2(\"$u\",800,720); return false;'>" .
            "<font color='green' size='-2'>Files</font></a>&nbsp;";
     	}
        $u = $_SERVER['SCRIPT_NAME'] ."?fcn=detail&amp;table=$tablenick&amp;id=" . $row[$samplespkey];
        $html .= "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
            "<font color='green' size='-2'>Details</font></a>&nbsp;";
        if ($GLOBS['iammgr']) {
            $u = $_SERVER['SCRIPT_NAME'] ."?fcn=edit&amp;table=$tablenick&amp;id=" . $row[$samplespkey];
            $html .= "<a href='$u' onclick='javascript:popup2(\"$u\",680,720); return false;'>" .
                "<font color='red' size='-2'>Edit</font></a>";
        }
        $html .= "</td></tr>\n";
    }
    $html .= "</table>\n" .
        "<div class='indent'>\n" . $GLOBS['statuscolor'] . "</div>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = ViewSampleDetail($tablenick, $id, $sql)
#   Show details for a particular entry
---------------------------------------------------------------*/
function ViewSampleDetail($tablenick, $id, $sql) {
	global $quickcols;

	$result = SQL_Query($sql);
	$row = SQL_Fetch($result);
	$numrows = SQL_NumRows($result);
	if ($numrows != 1) {
		Emsg("Unable to find detail for this entry");
		return '';
	}
	$conv2dat = array_keys($quickcols);		// Job management columns

	$collist = array_keys($row);
    foreach ($collist as $c) {
        $d = $row[$c];
        if ((! isset($d)) || (is_null($d)) || ($d == '') || ($d == 'NULL')) { $d = 'n/a'; }
        if ($c == 'dateinit') { $d = date('Y/m/d H:i', $d); }
        if ($c == 'datecomplete' && $d != '&nbsp;') { $d = date('Y/m/d H:i', $d); }
        if (in_array($c, $conv2dat)) {  // Special field needs formatting
            $val = DateState($d);
            if ($val != '') {
                $d = "<span class='$val'>$val</span>";

                $url = $_SERVER['PHP_SELF'] . "?id=$id&amp;table=$tablenick&amp;col=$c&amp;fcn=setrequest";
                $d .= "&nbsp;&nbsp;&nbsp;&nbsp;" . 
                    "<span class='requestcancel'>" .
                    "<a href='$url' onclick='javascript:confirmpopup(\"Request $c\"," .
                    "\"$url\", 640,480); return false;'>Request</a></span>" .
                    " / ";
                $url = $_SERVER['PHP_SELF'] . "?id=$id&amp;table=$tablenick&amp;col=$c&amp;fcn=setcancel";
                $d .= "<span class='requestcancel'>" .
                    "<a href='$url' onclick='javascript:confirmpopup(\"Cancel $c\"," .
                    "\"$url\", 640,480); return false;'>Cancel</a></span>";
            }
        }
        //  If the data has a newline in it do more than just show it
        $p = strpos($d, "\n");
        if ($p) { $d = "<pre>$d</pre>"; }
        $html .= "<tr><td align='left' valign='top' width='15%'><b>" . ucfirst($c) . ":</b></td>" .
            "<td>&nbsp;&nbsp;&nbsp;</td>";
        $html .= "<td align='left' colspan='4'>$d</td>";
        $html .= "</tr>\n";
    }
    $html .= "</table></div><br/>\n";
    $html .= "<font size='-1'><table border='0' align='center' width='80%'><tr>" .
		"<td align='left'>&nbsp;&nbsp;&nbsp;<a href='javascript:window.location.reload()'>Reload</a></td>" .
		"<td align='right'><a href='javascript:window.close()'>Close</a>&nbsp;&nbsp;&nbsp;</td>" .
		"</tr></table></font>\n";
    return $html;
}

?>
