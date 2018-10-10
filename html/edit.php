<?php
/*#################################################################
#
# Name:	edit.php
#
# Description:
#   Routines common to edit/modify function
#
#   Assumes these functions are available
#       Emsg
#       Nice_Exit
#
# Copyright (C) 2010 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/


/*---------------------------------------------------------------
# Edit - Generate HTML to edit a particular entry
#
# Parameters:
#   table - name of table of data
#   row - array of data to edit from database
#   showcols - array of columns names to be editted and order to show them
#   htmlcomments - hash of column names and HTML comments to for each
#
# Notes:
#   PHP_SELF invoked with fcn=modify&id=primarykeyvalue
#   PHP_SELF form variable names are in format IN[colname]
#
# Returns:
#   string of HTML
---------------------------------------------------------------*/
function Edit($table, $row, $showcols = array(), $htmlcomments = array()) {
    global $iammgr, $peditor;
    $c1 = "<br><font color='darkblue' size='-2'>&nbsp;&nbsp;&nbsp;&nbsp;";
    $c2 = '</font>';

    $DESC = SQL_Desc($table);               // Get array of description information
    if (! isset($DESC['_prikey_'])) { Nice_Exit("No primary key found for '$table'"); }
    $pkey = $DESC['_prikey_'];              // Name of primary key column
    //print "<!-- DESC=\n"; print_r($DESC); print " -->\n";

    //  These are the keys to be shown. Use all cols if nothing provided
    //  By default do not show primary key or any key ending with '_'
    if (! $showcols) {
        $showrols = array();
        while (list($k,$val) = each($row)) {
            if ($k == $pkey) { continue; }
            $x = strlen($k) - 1;
            if (substr($k,$x,1) == '_') { continue; }
            array_push($showcols, $k);
        }
    }
    
    $s = "<table border='0' width='80%' align='left'>\n" .
        "<form action='" . $_SERVER['PHP_SELF'] . "?fcn=modify' method='post'>\n" .
        "<input type='hidden' name='$pkey' value='" . $row[$pkey] . "'>\n";

    foreach ($showcols as $k) {
        $val = $DESC[$k];
        if (substr($k,0,1) == '_') { $ss = substr($k,1); }  // $ss is column name
        else { $ss = $k; }
        if (isset($htmlcomments[$k])) {
            $c = $htmlcomments[$k];
            if ($c) { $c = $c1 . $c . $c2; }    // Add wrapper HTML
        }
        else { $c = ''; }
        $s .= "<tr><td align='right'><b>" . ucfirst($ss) . "</b>&nbsp;</td>";
        //  Calculated fields are just shown
        if (substr($k,0,1) == '_' || (! $iammgr && ! $peditor)) {
            $s .= "<td>" . $row[$k] . "$c</td></tr>\n";
            continue;
        }

        //  Figure out width of text input
        $l = 6;
        if (preg_match('/string\s(\d+)/', $val, $m)) {
            $l = $m[1];
            if ($l > 64) { $l = 64; }           // Don't let it get too long
        }
        if (preg_match('/varchar\s(\d+)/', $val, $m)) {
            $l = $m[1];
            if ($l > 64) { $l = 64; }           // Don't let it get too long
        }
        if (preg_match('/text/', $val)) { $l = 65; }
        if (preg_match('/blob\s(\d+)/', $val, $m))   { $l = $m[1]; }
        if ($l > 64) { 
            $s .= "<td><textarea name='IN[$k]' rows='10' cols='80'>" . $row[$k] . "</textarea>$c</td></tr>\n";
        }
        else {
            $s .= "<td><input type='text' name='IN[$k]' size='$l' value='" . $row[$k] . "'>$c</td></tr>\n";
        }
    }
    $s .= "<input type='submit' value=' Modify Entry '>\n" .
        "<input type='button' name='cancel' value=' Cancel ' onclick='javascript:history.back()'>\n" .
        "<input type='button' name='close' value=' Close ' onclick='javascript:window.close()'>\n" .
        "</form>\n</table>\n<br><br>\n";
    return $s;
}

/*---------------------------------------------------------------
# Modify - Save $_POST[IN] fields in the database
#
# Parameters:
#   table - name of table of data
#   id - id to be updated  (format N or N.M)
#
# Notes:
#   _POST variable names are in format IN[colname]
#
# Returns:
#   Boolean for success
---------------------------------------------------------------*/
function Modify($table, $id) {

    if (! $_POST['IN']) { Emsg('Modify called with no data in $_POST'); return FALSE; }

    $DESC = SQL_Desc($table);           // Get array of description information
    if (! isset($DESC['_prikey_'])) { Nice_Exit("No primary key found for '$table'"); }
    $pkey = $DESC['_prikey_'];              // Name of primary key column

   //  Input fields are an array of arrays $_POST[IN[]]
    $inp = array();
    while (list($k,$val) = each($_POST['IN'])) {
        //  Sep 2018 datetime got picky. Ignore any of those fields
        if ( $DESC[$k] != 'datetime' ) { $inp[$k] = $val; }
    }
    $inp[$pkey] = $id;                      // Be sure primary key is set

    //  Set in defaults if not specified
    //$a = SetDefaults($inp);
    //if ($a) {
    //    reset($a);                          // Broken design
    //    while (list($k,$val) = each($a)) { $inp[$k] = $val; }
    //}
    //if ($_SESSION['id'] && ! VerifyCols($inp)) { return FALSE; }

    //  Call code to update some values as necessary. Can modify $np
    //if ($upd) {
    //    if ($VERBOSE > 1) { print "<!-- Update_INPUT_ROW=\n"; print_r($inp); print "-->\n"; }
    //    $inp = Local_Update($inp);
    //    if ($VERBOSE > 1) { print "<!-- Update_RETURN_ROW=\n"; print_r($inp); print "-->\n"; }
    //}

    //  Modify an existing row
    $sql = "UPDATE $table  SET ";
    while (list($k,$val) = each($inp)) {
        if ($k == $pkey) { continue; }
        $val = stripslashes($val);          // Canonicalize data first
        $val = SQL_Escape($val);
        
        //  July 2016 we suddenly need to quote fields
        if (substr($val,0,1) == "'") { $sql .= "$k=$val,"; }
        else { $sql .= "$k='$val',"; }
	}
    $sql = substr($sql,0,strlen($sql)-1);   // Drop last comma
    $sql .= " WHERE $pkey='$inp[$pkey]'";

    //  Do SQL we have constructed. Be sure we know when it fails
    $result = SQL_Query($sql, 0);
    $e = DB_CheckError($result);
    if ($e) {
        EMsg("Update Failed SQL=$sql<br/><br/>\n$e<br/>");
        return FALSE;
    }
    return TRUE;
}

?>
