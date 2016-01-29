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
$MON='topmed'; include_once '../local_config.php';
include_once '../common.php';
include_once '../header.php';
include_once '../DBMySQL.php';
include_once "../edit.php";
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

$LASTNDAYS = 20;            // How many days of data to plot

//  Handy javascript fragments
$JS['VSPACER'] = "<br/><br/>";
$JS['SPACER'] = "&nbsp;&nbsp;&nbsp;&nbsp;";
$JS['CLOSE'] = "<p align='right'><font size='-1'><a href='javascript:window.close()'>Close</a>" . $JS['SPACER'];
$JS['BACK'] = "<p align='right'><font size='-1'><a href='javascript:history.back()'>Back</a>" . $JS['SPACER'];

//-------------------------------------------------------------------
//  Get parameters passed in via normal invocation
//-------------------------------------------------------------------
//print doheader($HDR['title'], 1);

//  Real parameters for each form, default is ''
$parmcols = array('fcn');
extract (isolate_parms($parmcols));

DB_Connect($LDB['realm']);
if (! $fcn)    { $fcn = 'whatever'; }

$totalcompletedbams = '';
$sql = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . " WHERE state_ncbiorig=$COMPLETED";
$result = SQL_Query($sql);
$row = SQL_Fetch($result);
$totalcompletedbams .= $row['count(*)'] . '/';
$sql = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . " WHERE state_ncbib37=$COMPLETED";
$result = SQL_Query($sql);
$row = SQL_Fetch($result);
$totalcompletedbams .= $row['count(*)'] . '/';
$sql = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . " WHERE state_ncbib38=$COMPLETED";
$result = SQL_Query($sql);
$row = SQL_Fetch($result);
$totalcompletedbams .= $row['count(*)'];

$totalbamcount = 0;
$sql = 'SELECT * FROM ' . $LDB['stepstats'];
$result = SQL_Query($sql);
$numrows = SQL_NumRows($result);
$sqldata = array();                         // Save all SQL data
for ($i=0; $i<$numrows; $i++) {
    $row = SQL_Fetch($result);
    $totalbamcount = $row['bamcount'];
    array_push($sqldata, $row);
}
//  We have saved all SQL data in $sqldata


if ($fcn == 'whatever') {
    print doheader($HDR['title'], 1);
    print "<h2 align='center'>TopMed Activity</h2>\n";
    $yyyy = date('Y');
    $mm = date('m');
    $dd = sprintf('%02d', date('d') - $LASTNDAYS);
    if ($dd < 0) {
        $dd = sprintf('%02d', 30 - date('d'));
        $mm = sprintf('%02d', date('m') - 1);
    }
    if ($mm == '00') {
        $mm = '12';
        $yyyy--;
    }
    $NCBIBAMDATE = $yyyy . '/' . $mm . '/' . $dd;   // No BAMs sent to NCBI before this

    //-------------------------------------------------------------------
    //  Details about steps for processing each BAM (non-NCBI)
    //-------------------------------------------------------------------
    print "<h4>Processing Steps Before Sending to NCBI</h4>\n" .
        "<p>The following describe the various of steps to process a BAM.</p>\n";
    $legend = array('bams');
    $title = "Daily Count of Verified BAMs  Max=$totalbamcount";
    $plotdata = array();
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['bamcount']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend);

/*
    // Plot counts of steps not completed 
    $legend = array('ncbib38', 'ncbib37', 'ncbiorig', 'ncbiexpt', 'b38', 'b37', 'cram', 'qplot', 'bai', 'md5ver', 'arrive');    // Reversed
    $title = "Steps Not Completed for $totalbamcount BAMs";
    $plotdata = array(); 
    foreach ($legend as &$c) {  
        $sql = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . " WHERE state_$c!=$COMPLETED";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        $d = array();
        array_push($d, $c);
        array_push($d, $row['count(*)']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend, '', 'y', 'bars', 'text-data-yx');
*/

    // Plot totals of all states for each step
    $legend = array('ncbib38', 'ncbib37', 'ncbiorig', 'ncbiexpt', 'b38', 'b37', 'cram', 'qplot', 'bai');    // Reversed
    $title = "Counts for Each Step [$totalbamcount Verified BAMs]";
    $plotdata = array();
    $s = 'SELECT count(*) FROM ' . $LDB['bamfiles'] . ' ';
    foreach ($legend as &$c) {  
        $d = array();
        array_push($d, $c);
        $sql = $s . " WHERE state_$c=$COMPLETED";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        array_push($d, $row['count(*)']);

        $sql = $s . " WHERE state_$c=$FAILED OR state_$c=$CANCELLED";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        array_push($d, $row['count(*)']);

        $sql = $s . " WHERE state_$c=$SUBMITTED OR state_$c=$STARTED OR state_$c=$DELIVERED";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        array_push($d, $row['count(*)']);

        $sql = $s . " WHERE state_$c=$NOTSET OR state_$c=$REQUESTED";
        $result = SQL_Query($sql);
        $row = SQL_Fetch($result);
        array_push($d, $row['count(*)']);

        //print "<pre>$c Data=\n"; print_r($d); print "\n</pre>\n";
        array_push($plotdata, $d);
    }
    $legend = array();
    MakePlot($plotdata, $title, $legend, '', 'y', 'stackedbars', 'text-data-yx');
    print "<p><font color='DarkGreen'><b>Completed</b></font>, " .
        "<font color='red'><b>Failed/Canceled</b></font>, " .
        "<font color='SlateBlue'><b>In Process</b></font> (submitted, running), or " .
        "<font color='gray'><b>Not Started</b></font> " .
        "</br><br></p>\n";

    $legend = array('cram', 'qplot', 'bai', 'md5ver');
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
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend, 'Seconds');

    //-------------------------------------------------------------------
    //  Details about steps sending data to NCBI
    //-------------------------------------------------------------------
    print "<h4>Sending BAMs to NCBI</h4>\n" .
        "<p>The following describe the number of various BAMs and verage send times " .
        "for the three types of BAMs  <b>orig</b> are the original BAMs (BROAD " .
        "original BAMs are recreated from CRAM. <b>b37</b> BAMs are remapped using " .
        "build 37 and are created from CRAMs. <b>b38</b> is similar except using " .
        "build 38.</p>\n";
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

    $title = "Daily Count of Errors When Sending BAMs to NCBI";
    $legend = array('totalerrs');
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['errcount']);
        array_push($plotdata, $d);
    }
    MakePlot($plotdata, $title, $legend, '', 'y');

    exit;
}

//  What was that?
Emsg("Unknown directive '$fcn'.");
Nice_Exit("How'd you do that?");
exit;

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
    global $JS;
    $ymax = 240;
    $xmax = 600;

    $plot = new PHPlot($xmax, $ymax);
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
        $plot->SetDataColors(array('DarkGreen', 'red', 'SlateBlue', 'gray'));
        if ($ypoints) { $plot->SetXDataLabelPos('plotin'); }
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
    print $JS['VSPACER'] . "\n";
}

?>
