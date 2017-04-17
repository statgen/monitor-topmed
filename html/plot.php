<?php
/*#################################################################
#
# Name: plot.php
#
# Description:
#   Generate plots of interesting for tracking NHLBI TopMed data
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
include_once 'phplot.php';

//print "<!-- _POST=\n"; print_r($_POST); print " -->\n";

$NOTSET = 0;                // Not set
$REQUESTED = 1;             // Task requested
$SUBMITTED = 2;             // Task submitted to be run
$STARTED   = 3;             // Task started
$DELIVERED = 19;            // Data delivered, but not confirmed
$COMPLETED = 20;            // Task completed successfully
$CANCELLED = 89;            // Task cancelled
$FAILED    = 99;            // Task failed

$YMAX = 240;                // Globals for screen size
$XMAX = 600;

$COLORS = array('DarkGreen', 'red', 'orange', 'SlateBlue', 'gray');

//-------------------------------------------------------------------
//  Get parameters passed in via normal invocation
//-------------------------------------------------------------------
print doheader($HDR['title'], 1);

//  Real parameters for each form, default is ''
$parmcols = array('fcn', 'lastdays', 'datayear');
extract (isolate_parms($parmcols));


DB_Connect($LDB['realm']);
if (! $fcn)    { $fcn = 'plot'; }
if (! $lastdays)   { $lastdays = 20; }
if (! $datayear)   { $datayear = 1; }

$totalcompletedbams = '';
$sql = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . " WHERE datayear=$datayear AND state_ncbiorig=$COMPLETED";
$result = SQL_Query($sql);
$row = SQL_Fetch($result);
$totalcompletedbams .= $row['count(*)'] . '/';
$sql = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . " WHERE datayear=$datayear AND state_ncbib37=$COMPLETED";
$result = SQL_Query($sql);
$row = SQL_Fetch($result);
$totalcompletedbams .= $row['count(*)'] . '/';
$sql = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . " WHERE datayear=$datayear AND state_ncbib38=$COMPLETED";
$result = SQL_Query($sql);
$row = SQL_Fetch($result);
$totalcompletedbams .= $row['count(*)'];

$sql = 'SELECT * FROM ' . $LDB['stepstats'];
$result = SQL_Query($sql);
$numrows = SQL_NumRows($result);
$sqldata = array();                         // Save all SQL data
for ($i=0; $i<$numrows; $i++) {
    $row = SQL_Fetch($result);
    array_push($sqldata, $row);
}

$sql = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . " WHERE datayear=$datayear AND state_verify=$COMPLETED";
$result = SQL_Query($sql);
$row = SQL_Fetch($result);
$totalbamcount = $row['count(*)'];


//-------------------------------------------------------------------
//  Generate all interesting plots
//-------------------------------------------------------------------
if ($fcn == 'plot') {
    $NCBIBAMDATE = date('Y/m/d', strtotime("-$lastdays days"));
    //print "NCBIBAMDATE=$NCBIBAMDATE<br>\n";
    //-------------------------------------------------------------------
    //  Details about steps for processing each BAM (non-NCBI)
    //-------------------------------------------------------------------
    $reshowurl =  $_SERVER['SCRIPT_NAME'] . "?datayear=$datayear&amp;lastdays=$lastdays";
    print "<table width='80%' align='center' border='0'> <tr>\n" .
        "<td align='left'>" .
          "<form action='" . $_SERVER['PHP_SELF'] . "' method='post'>\n" .
          "<b><input type='submit' value=' Generate Plots: '>\n" .
          "Last <input type='text' name='lastdays' value='$lastdays' size='2'>\n" .
          "days</td>" .
        "<td><b>Project:</b>" .
          "<select name='datayear'>" .
          "<option value='2'>Year 2</option><option value='1'>Year 1</option></select>\n" . 
          "</td>" .
        "<td align='right'>" .
        "<a href='$reshowurl'>Reshow Plots</a>" .
        "</td></tr></form></table></b></br>\n";   
        
    // Plot totals of all states for each step
    print "<p><b>Legend: <font size='-1'>" .
        "<font color='$COLORS[0]'>Completed</font>, " .
        "<font color='$COLORS[1]'>Failed/Canceled</font>, " .
        "<font color='$COLORS[2]'>Not Loaded</font>, " .
        "<font color='$COLORS[3]'>In Process</font> (submitted, running), or " .
        "<font color='$COLORS[4]'>Not Started</font> " .
        "</font></b></p>\n";
    $legend = array('bcf', 'gce38post', 'gce38pull', 'gce38push',
        'b38', 'b37', 'cram', 'qplot', 'bai');    // Reversed
    $title = "Current Counts for Each Step [$totalbamcount Verified BAMs]";
    $plotdata = array();
    $s = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . " WHERE datayear=$datayear AND";
    foreach ($legend as &$c) {  
        $d = array();
        array_push($d, $c);
        $sql = $s . "  state_$c=$COMPLETED";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        array_push($d, $row['count(*)']);

        $sql = $s . " state_$c=$FAILED OR state_$c=$CANCELLED";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        array_push($d, $row['count(*)']);

        $sql = $s . " state_$c=$DELIVERED";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        array_push($d, $row['count(*)']);

        $sql = $s . " state_$c=$SUBMITTED OR state_$c=$STARTED OR state_$c=$REQUESTED";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        array_push($d, $row['count(*)']);

        $sql = $s . " state_$c=$NOTSET";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        array_push($d, $row['count(*)']);

        //print "<pre>$c Data=\n"; print_r($d); print "\n</pre>\n";
        array_push($plotdata, $d);
    }
    $legend = array();
    $xtmp = $XMAX;
    $XMAX = 900;                // Extra wide so we can read numbers
    MakePlot($plotdata, $title, $legend, '', '', 'stackedbars', 'text-data-yx');
    $XMAX = $xtmp;

    $legend = array('cram', 'qplot', 'verify', 'bcf');
    $title = "Daily Count of Steps Completed";
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['count_cram']);
        array_push($d, $row['count_qplot']);
        array_push($d, $row['count_bai']);
        array_push($d, $row['count_verify']);
        array_push($d, $row['count_bcf']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend);

    $title = "Ave Completion Time/Step";
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['avetime_cram']);
        array_push($d, $row['avetime_qplot']);
        array_push($d, $row['avetime_bai']);
        array_push($d, $row['avetime_verify']);
        array_push($d, $row['avetime_bcf']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend, 'Seconds');

    //-------------------------------------------------------------------
    //  Details about steps sending data to NCBI
    //-------------------------------------------------------------------
    print "<h4>Google Cloud Activity</h4>\n";
    $title = "Daily Count of Samples PUSHed to GCE";
    $legend = array('push', 'pull', 'post');
    $title = "Daily Count of Steps Completed";
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['count_gcepush']);
        array_push($d, $row['count_gcepull']);
        array_push($d, $row['count_gcepost']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend);

    $title = "Ave Time for GCE-Related Steps";
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['avetime_gcepush']);
        array_push($d, $row['avetime_gcepull']);
        array_push($d, $row['avetime_gcepost']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend, 'Seconds', 'y');

    print "<p align='right'><a href='$reshowurl'>Reshow Plots</a>\n";
    print dofooter($HDR['footer']);
    exit;
}

//  What was that?
Emsg("Unknown directive '$fcn'.");
Nice_Exit("How'd you do that?");
exit;

/*    DEAD CODE

    //-------------------------------------------------------------------
    //  Details about steps sending data to NCBI
    //-------------------------------------------------------------------
    print "<h4>Sending BAMs to NCBI</h4>\n" .
        "<p><font size='-1'>" .
        "The following describe the number of various BAMs and average send times " .
        "for the three types of BAMs:  <b>orig</b> are the original BAMs (BROAD " .
        "original BAMs are sent as a CRAM). <b>b37</b> are BAMs remapped using " .
        "build 37 and are created from CRAMs and <b>b38</b> is similar except using " .
        "build 38." .
        "</font></p>\n";
    $title = "Daily Count of BAMs Sent to NCBI";
    $legend = array('orig', 'b37', 'b38');
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['count_orig']);
        array_push($d, $row['count_b37']);
        array_push($d, $row['count_b38']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend, '', 'y');

    $title = "Ave Time to Send BAM to NCBI";
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['avetime_orig']);
        array_push($d, $row['avetime_b37']);
        array_push($d, $row['avetime_b38']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend, 'Seconds', 'y');

    $title = "Daily Count of BAMs loaded at NCBI   Totals=$totalcompletedbams";
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['loadedorigbamcount']);
        array_push($d, $row['loadedb37bamcount']);
        array_push($d, $row['loadedb38bamcount']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend, '', 'y');

    print "<p><font size='-1'>" .
        "Failure rates sending and retrieving data from GCE." .
        "</font></p>\n";
    $title = "Daily Count of Errors for GCE Data";
    $legend = array('orig', 'origchecksum', 'b37', 'b37checksum', 'b38', 'b38checksum');
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['errorigcount']);
        array_push($d, $row['errckorigcount']);
        array_push($d, $row['errb37count']);
        array_push($d, $row['errckb37count']);
        array_push($d, $row['errb38count']);
        array_push($d, $row['errckb38count']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend, '', 'y');

*/

/*---------------------------------------------------------------
#   MakePlot($plotdata, $title, $legend, $ytitle, $ypoints, $type, $datatype)
#
#   Generate a plot in the current HTMNL stream
#   $ytitle could be the title onn the Y axis
#   $ypoints is a boolean if the Y values should be annotated on the plot
#   $type should be lines or bars
#   $datatype should be text-data or text-data-yx
---------------------------------------------------------------*/
function MakePlot($plotdata, $title, $legend, $ytitle='', $ypoints='', $type='lines', $datatype='text-data') {
    global $JS, $YMAX, $XMAX, $COLORS;

    $plot = new PHPlot($XMAX, $YMAX);
    $plot->SetFailureImage(False);  // No error images
    $plot->SetPrintImage(False);    // No automatic output
    $plot->SetImageBorderType('plain');
    $plot->SetPlotType($type);
    $plot->SetDataType($datatype);
    $plot->SetDataValues($plotdata);
    $plot->SetLineWidths(3);
    if ($datatype == 'text-data-yx') {
        $plot->SetYTickPos('none');
        //$plot->SetXTickPos('none');
        //$plot->SetXTickLabelPos('none');
        //$plot->SetXDataLabelPos('plotin');
        //$plot->SetLegendPixels($xmax-85, 5);
        $plot->SetDataColors($COLORS);
        if ($ypoints) { $plot->SetXDataLabelPos('plotin'); }
        else { $plot->SetXDataLabelPos('plotstack'); }
    }
    else {
        $plot->SetLegend($legend);
        $plot->SetXLabelAngle(90);
        $plot->SetLegendPixels(45, 25);
        //$plot->SetLegendPosition(0, 0, 'plot', 0, 0, 5, 5);
        if ($ypoints) { $plot->SetYDataLabelPos('plotin'); }
    }
    $plot->SetTitle($title);
    if ($ytitle) {
        $plot->SetYTitle($ytitle);
        // With Y data labels, we don't need Y ticks or their labels, so turn them off.
        //$plot->SetYTickLabelPos('none');
        //$plot->SetYTickPos('none');
    }
    $plot->DrawGraph();
    print "<img src=\""; print $plot->EncodeImage(); print "\" alt='$title'>\n";
    print "<br/><br/>\n";
}

?>
