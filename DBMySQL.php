<?php
/*#################################################################
#
# Name: DBMySQL.php
#
# Description:
#   Routines for accessing a MySQL database
#   I.e a poor man's Pear::DB cause some ISPs don't install Pear::DB
#   This was converted to use mysqli procedures
#
#   Assumes these external routines:
#   Nice_Exit(msg)
#
#   References these global variables:
#   $dbh - resource set for connection
#   $SQL_ERROR - set here with error message from MySQL
#
# Copyright (C) 2010,2014 Terry Gliedt
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/

/*---------------------------------------------------------------
# DB_Connect - Connect to a database, using persistent connect
#
# Parameters:
#   realm - realm file, containing the connection information
#           See Perl package DBIx::Connector
#
# Returns:
#   nothing, sets global $dbh
---------------------------------------------------------------*/
function DB_Connect($realm) {
    global $dbh;

    $lines = file($realm);
    foreach ($lines as &$l) {
        if (preg_match('/DBD=(\S+)/', $l, $m)) { $type = $m[1]; continue; }
        if (preg_match('/SERVER=host=(\S+)/', $l, $m)) { $host = $m[1]; continue; }
        if (preg_match('/USER=(\S+)/', $l, $m)) { $user = $m[1]; continue; }
        if (preg_match('/PASS=(\S+)/', $l, $m)) { $pass = $m[1]; continue; }
        if (preg_match('/DATABASE=(\S+)/', $l, $m)) { $db = $m[1]; continue; }
    }
    if ($type == '') { Nice_Exit("Failed to read realm file '$realm'"); }
    if ($type != 'mysql') { Nice_Exit("This code does not support anything but MySQL"); }

    $dbh = mysqli_connect($host, $user, $pass, $db);
    if ($dbh && (! mysqli_connect_errno())) { return $dbh; }
    $msg = "Unable to connect to database<BR/>\n" .
        mysqli_connect_errno() . '" ' . mysqli_connect_errno() . "<br/>\n";
    //$msg .= mysqli_info(). "<br/>\n";   // Do not show password usually
    Nice_Exit($msg);
}

/*---------------------------------------------------------------
# SQL_Query - Issue an SQL statement.  Uses $dbh, Sets $SQL_ERROR
#
# Parameters:
#   sql - SQL to execute
#   dieonerror - If set and an error occurs, dies
#
# Returns:
#   result from mysql_query or false and $SQL_ERROR is set
---------------------------------------------------------------*/
function SQL_Query($sql, $dieonerror=1) {
    global $dbh, $SQL_ERROR;
    $SQL_ERROR = '';

    $result = mysqli_query($dbh, $sql);
    if ($result) { return $result; }

    $SQL_ERROR = "SQL prepare/stmt failed. SQL=$sql<br/>\n" .
        'Error=' . mysqli_connect_errno() . '" ' . mysqli_connect_errno() . "<br/>\n";
    //$msg .= mysqli_info(). "<br/>\n";   // Do not show password usually
    if ($dieonerror) { Nice_Exit($SQL_ERROR); }
    return 0;
}

/*---------------------------------------------------------------
# DB_CheckError - Return error message if there was an error
#
# Parameters:
#   dbo - object for DB action
#
# Returns:
#   MySQL error message or null string
---------------------------------------------------------------*/
function DB_CheckError($dbo) {
   global $dbh;
    if (mysqli_error($dbh)) {
        if ($dbo) { return $dbo->getMessage(); }
        return 'Previous SQL command found no results';
    }
    return '';
}

/*---------------------------------------------------------------
# DB_IsError - If an error was detected, generate msg and exit
#
# Parameters:
#   dbo - object for DB action
#   table - name of table
#
# Returns:
#   nothing or does not return
---------------------------------------------------------------*/
function DB_IsError($dbo, $table='unknown') {
   global $dbh;
    if (! mysqli_error($dbh)) { return; }
    $m = "DB failure for <b>'$table'</b> " . $dbo->getMessage() . "<br/>";
    $m .= $dbo->getDebugInfo();     // Not for production mode
    Emsg($m);
    Nice_Exit("<br/><br/>Unable to update database");
}

/*---------------------------------------------------------------
# SQL_Fetch - Get a result row as an associative array
#s
# Parameters:
#   res - SQL query resource
#
# Returns:
#   array
---------------------------------------------------------------*/
function SQL_Fetch($res) {
    return mysqli_fetch_assoc($res);
}

/*---------------------------------------------------------------
# SQL_Query_ID - Get auto generated id usedin the last query
#
# Parameters:
#   dieonerror - If set and an error occurs, dies
#
# Returns:
#   integer or zero
---------------------------------------------------------------*/
function SQL_Query_ID($dieonerror=1) {
    global $dbh, $SQL_ERROR;
    $i = mysqli_insert_id($dbh);
    if ($i) { return $i; }
    $SQL_ERROR = "SQL_Query_ID did not return a value. No AUTO_INCREMENT found<br/>\n";
    if ($dieonerror) { Nice_Exit($SQL_ERROR); exit; }
	return $i;
}

/*---------------------------------------------------------------
# SQL_NumRows - Get number of rows in result
#
# Parameters:
#   res - SQL query resource
#
# Returns:
#   integer
---------------------------------------------------------------*/
function SQL_NumRows($res) {
    return mysqli_num_rows($res);
}

/*---------------------------------------------------------------
# SQL_Escape - Escapes special characters so it is safe in an SQL query
#
# Parameters:
#   s - string
#
# Returns:
#   escaped string
---------------------------------------------------------------*/
function SQL_Escape($s) {
    global $dbh;
    $ss = mysqli_real_escape_string($dbh, $s);
    if ($ss) { return "'" . $ss . "'"; }
    return "'$s'";                  // Helps when we have a syntax error
}

/*---------------------------------------------------------------
# SQL_Desc - Returns an array of columns/datatypes and primary key info
#
# Parameters:
#   table - table to check
#
# Returns:
#   array of $DESC information like this:
#
#   [id] => int 11
#   [_prikey_] => id            Primary column (if any)
#   [_inc_] => id               Auto increment (if any)
#   [dirname] => string 255
#   [status] => string 16
#   [comments] => blob 65535
#   [count] => int 32
#
---------------------------------------------------------------*/
function SQL_Desc($table) {
    $DESC = array();

    $result = SQL_Query("DESCRIBE $table");
    DB_IsError($result);
    $numrows = SQL_NumRows($result);

    for ($i=0; $i<$numrows; $i++) {
        $row = SQL_Fetch($result);
        $c = $row['Field'];
        $t = $row['Type'];
        //  E.g. varchar(32) => varchar 32, int 11  etc
        if (preg_match('/(\S+)\((\d+)\)/', $row['Type'], $m)) {
            $t = $m[1] . ' ' . $m[2];
        }
        $DESC[$c] = $t;
        //  Check for special cases
        if ($row['Key'] == 'PRI') { $DESC['_prikey_'] = $c; continue; }
        if (preg_match('/auto_increment/', $row['Extra'])) { $DESC['_inc_'] = $c; }
	}
    return $DESC;
}

?>
