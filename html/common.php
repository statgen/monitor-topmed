<?php
/*#################################################################
#
# Name:	common.php
#
# Description:
#   Routines common to monitor code
#
# Copyright (C) 2010- University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/

/*---------------------------------------------------------------
# Emsg - Displays message in red
#
# Parameters:
#   msg - message to issue
#   flag - if set return string, else print it
---------------------------------------------------------------*/
function Emsg($msg, $f = 0) {
    $m = "<br/><font color='red'><b>$msg</b></font><br/>\n";
    if ($f) { return $m; }
    print $m;
}

/*---------------------------------------------------------------
# Msg - Displays message in blue
#
# Parameters:
#   msg - message to issue
#   flag - if set return string, else print it
---------------------------------------------------------------*/
function Msg($msg, $f= 0) {
    $m = "<br/><font color='darkblue'><b>$msg</b></font><br/>\n";
    if ($f) { return $m; }
    print $m;
}

/*---------------------------------------------------------------
# Nice_Exit - Displays a message and exits
#
# Parameters:
#   msg - message to issue
#
# Returns:
#   Does not return
---------------------------------------------------------------*/
function Nice_Exit($msg = '') {

    print "<br/><br/><b>$msg</b><br/><br/>\n" .
	"<table align='center' width='80%' border='0'><tr>" .
	"<td align='left'><a href='javascript:close()'>Close this window</a></td>" .
	"<td align='right'><a href='javascript:history.back()'>Return to correct</a></td>" .
	"</tr></table>\n";
    exit;
}

/*---------------------------------------------------------------
# isolate_parms - Given a list of parameter keywords, build
#       a string of PHP code which can be evaluated.
#
# Parameters:
#   aref - array of parameters from POST or GET or _SESSION
#
# Returns:
#   extract input
---------------------------------------------------------------*/
function isolate_parms($aref) {

    $myvars = array();
    while (list ($key, $val) = each($aref)) {
        $v = '';
        if (isset($_POST[$val])) { $v=$_POST[$val]; }
        elseif (isset($_GET[$val])) { $v=$_GET[$val]; }
        elseif (isset($_SESSION[$val])) { $v=$_SESSION[$val]; }
        $myvars[$val]=$v;
        //      Maybe this should be done in some cases?
        //$myvars[$val]=htmlspecialchars($myvars[$val]);        // Keep some characters
        $myvars[$val]=strip_tags($myvars[$val]);   // No HTML tags allowed
    }
    return $myvars;
}

?>
