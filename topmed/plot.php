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


$HTML = '';

$sql = 'SELECT * FROM ' . $LDB['stepstats'];
$result = SQL_Query($sql);
$numrows = SQL_NumRows($result);
$sqldata = array();                         // Save all SQL data
for ($i=0; $i<$numrows; $i++) {
    $row = SQL_Fetch($result);
    array_push($sqldata, $row);
}
//  We have saved all SQL data in $sqldata


if ($fcn == 'whatever') {
    print doheader($HDR['title'], 1);
    print "<h2 align='center'>TopMed Activity</h2>\n";
    $NCBIBAMDATE = '2016/01/01';            // No BAMs sent to NCBI before this

/*
  yyyymmdd CHAR(10) NOT NULL,
  count_expt        INT DEFAULT 0,
  avetime_expt      INT DEFAULT 0,
  ncbicount_expt    INT DEFAULT 0,
  count_orig        INT DEFAULT 0,
  avetime_orig      INT DEFAULT 0,
  ncbicount_orig    INT DEFAULT 0,
  count_b37         INT DEFAULT 0,
  avetime_b37       INT DEFAULT 0,
  ncbicount_b37     INT DEFAULT 0,
  count_b38         INT DEFAULT 0,
  avetime_b38       INT DEFAULT 0,
  ncbicount_b38     INT DEFAULT 0,
*/

    //-------------------------------------------------------------------
    //  Details about steps for processing each BAM (non-NCBI)
    //-------------------------------------------------------------------
    print "<h4>Processing Steps Before Sending to NCBI</h4>\n" .
        "<p>The following describe the various of steps completed per day " .
        "and the average time per step.</p>\n";
    $legend = array('verify', 'bai', 'qplot', 'cram');
    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['count_verify']);
        array_push($d, $row['count_bai']);
        array_push($d, $row['count_qplot']);
        array_push($d, $row['count_cram']);
        array_push($plotdata, $d);
    }
    $title = "Count of Steps Completed per BAM";
    MakePlot($plotdata, $title, $legend);

    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['avetime_verify']);
        array_push($d, $row['avetime_bai']);
        array_push($d, $row['avetime_qplot']);
        array_push($d, $row['avetime_cram']);
        array_push($plotdata, $d);
    }
    $title = "Ave Seconds for Steps Completed per BAM";
    MakePlot($plotdata, $title, $legend);

    //-------------------------------------------------------------------
    //  Details about steps sending data to NCBI
    //-------------------------------------------------------------------
    print "<h4>Sending BAMs to NCBI</h4>\n" .
        "<p>The following describe the number of various BAMs and verage send times " .
        "for the three types of BAMs  <b>orig</b> are the original BAMs (BROAD " .
        "original BAMs are recreated from CRAM. <b>b37</b> BAMs are remapped using " .
        "build 37 and are created from CRAMs. <b>b38</b> is similar except using " .
        "build 38.</p>\n";
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
    $title = "Count of BAMs Sent to NCBI";
    MakePlot($plotdata, $title, $legend);

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
    $title = "Ave Time to Send BAM to NCBI";
    MakePlot($plotdata, $title, $legend);

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
    $title = "Count of loaded BAMs at NCBI";
    MakePlot($plotdata, $title, $legend);

    $plotdata = array(); 
    for ($i=0; $i<$numrows; $i++) {
        $row = $sqldata[$i];
        if ($row['yyyymmdd'] < $NCBIBAMDATE) { continue; }
        $d = array();
        array_push($d, substr($row['yyyymmdd'],5,5));
        array_push($d, $row['errorigcount']);
        array_push($d, $row['errb37count']);
        array_push($d, $row['errb38count']);
        array_push($plotdata, $d);
    }
    $title = "Count of Errors Sending BAMs to NCBI";
    MakePlot($plotdata, $title, $legend);

    exit;
}

//  What was that?
Emsg("Unknown directive '$fcn'.");
Nice_Exit("How'd you do that?");
exit;

/*---------------------------------------------------------------
#   MakePlot($plotdata, $title, $legend)
#
#   Generate a plot in the current HTMNL stream
---------------------------------------------------------------*/
function MakePlot($plotdata, $title, $legend) {
    global $JS;

    $plot = new PHPlot(600, 240);
    $plot->SetFailureImage(False);  // No error images
    $plot->SetPrintImage(False);    // No automatic output
    $plot->SetImageBorderType('plain');
    $plot->SetPlotType('lines');
    $plot->SetDataType('text-data');
    $plot->SetDataValues($plotdata);
    $plot->SetXLabelAngle(90);
    $plot->SetLineWidths(3);
    $plot->SetLegend($legend);
    $plot->SetLegendPixels(45, 25);
    $plot->SetTitle($title);
    $plot->DrawGraph();
    print "<img src=\""; print $plot->EncodeImage(); print "\" alt='$title'>\n";
    print $JS['VSPACER'] . "\n";
}

?>
