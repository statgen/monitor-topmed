<?php
/*#################################################################
#
# Name: quickletters.php
#
# Description:
#   Code to support mapping of short code letters to actions
#
# Copyright (C) 2019 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/
$statusletters = "<i><b>a</b>=File Arrived, <b>5</b>=MD5 Verified, <b>B</b>=Local Backup Done, " .
    //"<b>A</b></b>=Upload data to AWS,<br/>" .
    "<b>F</b>=FIX";
if ($LDB['datatype'] == 'genome') {
	$statusletters = "<i><b>a</b>=File Arrived, <b>5</b>=MD5 Verified, <b>B</b>=Local Backup of CRAM, <b>C</b>=BAM=>CRAM, <b>Q</b>=qplot run,<br/>" .
    	"<b>7</b>=Remapped Build=37, " .
    	"<b>s</b>=Push Build=38 to GCE, <b>r</b>=Pull Build=38 from GCE, <b>8</b>=Remapped Build=38," .
    	"<br/><b>V</b>=Completed BCF/VT 38,<b>G</b></b>=Upload CRAM to GCE," .
    	"<b>g</b></b>=Upload BCF to GCE,<b>x</b></b>=Cleanup Backup files,<b>A</b></b>=Upload data to AWS,<br/>" .
    	"<b>X</b>=EXPT=>NCBI <b>S</b>=Orig BAM/CRAM=>NCBI, <b>P</b>=</b>B37=>NCBI<br/>" .
    	"<b>F</b>=FIX";
}

$NOTSET = 0;                			// Not set
$REQUESTED = 1;             			// Task requested
$SUBMITTED = 2;             			// Task submitted to be run
$STARTED   = 3;             			// Task started
$DELIVERED = 19;            			// Data delivered, but not confirmed
$COMPLETED = 20;            			// Task completed successfully
$IGNORETHIS = 80;           			// Task is to be ignored
$CANCELLED = 89;            			// Task cancelled
$FAILEDCHECKSUM = 98;       			// Task failed, because checksum at NCBI bad
$FAILED    = 99;            			// Task failed
$state2str = array(         			// Values here are class for SPAN tag
    $NOTSET => 'notset',
    $REQUESTED => 'requested',
    $SUBMITTED => 'submitted',
    $STARTED => 'started',
    $DELIVERED => 'delivered',
    $COMPLETED => 'completed',
    $CANCELLED => 'cancelled',
    $IGNORETHIS => 'ignorethis',
    $FAILEDCHECKSUM => 'failedchecksum',
    $FAILED => 'failed'
);

$quickcols = array(                     // Map of status column to topmedcmd verb
    'state_arrive'   => 'arrived',
    'state_verify'   => 'verify',
    'state_backup'   => 'backup',
    'state_aws38copy'=> 'awscopy'
);
$quickletter = array(                   	// Map of status column to letter we see
	'state_arrive'   => 'a',
	'state_verify'   => '5',
	'state_backup'   => 'B',
	'state_aws38copy'=> 'A'
);
$separator_actions = array();
if ($LDB['datatype'] == 'genome') {
	$quickcols = array(                     // Map of status column to topmedcmd verb
		'state_arrive'   => 'arrived',
		'state_verify'   => 'verify',
		'state_cram'     => 'cramed',
		'state_backup'   => 'backup',
		'state_qplot'    => 'qplot',
		'state_b37'      => 'mapping37',
		'state_gce38push'=> 'gcepush',
		'state_gce38pull'=> 'gcepull',
		'state_b38'      => 'mapping38',
		'state_gce38bcf' => 'bcf',
		'state_gce38copy'=> 'gcecopy',
		'state_gce38cpbcf'=> 'gcecpbcf',
		'state_gcecleanup'=> 'gcecleanup',
		'state_aws38copy'=> 'awscopy',
		'state_ncbiexpt' => 'sendexpt',
		'state_ncbiorig' => 'sendorig',
		'state_ncbib37'  => 'sendb37',
		'state_fix'      => 'fix'
	);
	$quickletter = array(                   // Map of status column to letter we see
    	'state_arrive'   => 'a',
    	'state_verify'   => '5',
    	'state_cram'     => 'C',
    	'state_backup'   => 'B',
    	'state_qplot'    => 'Q',
    	'state_b37'      => '7',
    	'state_gce38push'=> 's',
    	'state_gce38pull'=> 'r',
    	'state_b38'      => '8',
    	'state_gce38bcf' => 'V',
    	'state_gce38copy'=> 'G',
    	'state_gce38cpbcf' => 'g',
    	'state_gcecleanup' => 'x',
    	'state_aws38copy'=> 'A',
    	'state_ncbiexpt' => 'X',
    	'state_ncbiorig' => 'S',
    	'state_ncbib37'  => 'P',
    	'state_fix'      => 'F'
	);
    $separator_actions = array('Q','7','8','V', 'x', 'A', 'P');
}


$GLOBS['statuscolor'] = "<p><b>Note:</b><br>" .
    "<i>Colors:</i> " .
    "<span class='cancelled'> cancelled </span>&nbsp;" .
    "<span class='completed'> completed </span>&nbsp;" .
    "<span class='delivered'> delivered </span>&nbsp;" .
    "<span class='ignorethis'> ignore </span>&nbsp;" .
    "<span class='failedchecksum'> failed_checksum </span>&nbsp;" .
    "<span class='failed'> failed </span>&nbsp;" .
    "<span class='requested'> requested </span>&nbsp;" .
    "<span class='started'> started </span>&nbsp;" . 
    "<span class='submitted'> submitted </span>&nbsp;" .
    "<br>" .
    "<i>Status:</i> $statusletters<br>" .
    "</p>\n";

$GLOBS['statusruns'] = "<p><b>Note:</b><br>" .
    "<i>Status of data:</i>&nbsp;" .
    "<span class='done'> finished </span>&nbsp;" .
    "<span class='processing'> being processed </span>&nbsp;" .
    "<span class='unknown'> unstarted, etc </span>&nbsp;" .
    "<span class='failed'> at least one failure </span>&nbsp;" .
    "</p>\n";

/*---------------------------------------------------------------
# href = QuickStatus($r, $url)
#   $r = row of data for this BAM
#   $url is a url to wrap around failed states. XX needs to be replaced
#   Generate a short summary of the state for this directory
#   Parameter is the row from the database
#   Returns HTML
---------------------------------------------------------------*/
function QuickStatus($r, $url) {
    global $quickcols, $quickletter, $separator_actions;
    //  Add a small separator to 'group' certain actions
    $h = '';
    $col = '';
    $span='notset';
    $cols = array_keys($quickcols);
    foreach ($cols as $c) {
        $val = DateState($r[$c]);
        //  For those marked as failed, set up link
        $s =  $quickletter[$c];
        if ($val == 'failed' || $val == 'started') {
            $s = str_replace('XX', $quickletter[$c], $url);
        }
        if (in_array($s, $separator_actions)) { $val .= " separator"; }
        $h .= "<span class='$val'>&nbsp;" . $s . "&nbsp;</span>";
    }
    return $h;
}

/*---------------------------------------------------------------
# html = CalcRunStatus($str)
#   Convert status string into shorthand status
#   str should look like: A=done,5=processing,B=unknown etc
#   returns string of html
---------------------------------------------------------------*/
function CalcRunStatus($str) {
    global $quickletter, $separator_actions;
    //return $str;          // To see original state
    $h = '';
    $cols = array_values($quickletter);
    foreach ($cols as $c) {
        $span='notset';
        if (preg_match("/$c=([^,]+),/", $str, $m)) { $span = $m[1]; }
        //  Add a separator class after certain actions
        if (in_array($c, $separator_actions)) { $span .= " separator"; }
        $h .= "<span class='$span'> $c </span>";
    }
    return $h;
}

/*---------------------------------------------------------------
# href = DateState($t)
#   Return state for a particular time. See values at top
#   returns array (state)
---------------------------------------------------------------*/
function DateState($t) {
    global $state2str;
    $state = 'notset';
    if (in_array($t, $state2str)) { $state = $state2str[$t]; }
    if (array_key_exists($t, $state2str)) { $state = $state2str[$t]; }
    return $state;
}

?>
