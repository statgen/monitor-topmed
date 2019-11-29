<?php
/*#################################################################
#
# Name: immgrfcns.php 
#
# Description:
#   Code to support functions for managers
#
# Copyright (C) 2019 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/

/*---------------------------------------------------------------
#   ImmgrFunctions($fcn)
#   Handle functions related to RNA processing. Might not return
---------------------------------------------------------------*/
function ImmgrFunctions($fcn) {
    global $HDR, $HTML, $LDB, $PARMS;
	$JS['SPACER'] = "&nbsp;&nbsp;&nbsp;&nbsp;";
	$JS['CLOSE'] = "<p align='right'><font size='-1'><a href='javascript:window.close()'>Close</a>" . $JS['SPACER'];

	if ($fcn == 'edit') {
		$table = $PARMS['table'];
		$id = $PARMS['id'];
		$t = $LDB[$table];
		$pkey = $table . '_pkey';
		$pkey = $LDB[$pkey];
		$sql = "SELECT * FROM $t WHERE $pkey='$id'";
		$result = SQL_Query($sql);
		$row = SQL_Fetch($result);

		$HTML .= "<h3 align='center'>Edit Database Entry [table=$t entry=$id]</h3>\n" .
			'<center>' . Emsg("Just because you CAN change a field does not mean you should", 1) . "<br/></center>\n" .
			"<div class='indent'>\n";
		$HTML .= Edit($t, $row) . "</div>\n";
		print $HTML;
		print dofooter($HDR['footer']);
		exit;
	}
	if ($fcn == 'modify') {
		$table = $PARMS['table'];
		$id = $PARMS['id'];
		if (Modify($table, $id)) { $HTML .= Msg("Record '$id' modified", 1); }
		else { $HTML .= Emsg("Failed to modify record '$id'", 1); }
		print $HTML .
			"<p align='right'><font size='-1'><a href='javascript:window.close()'>Close</a>&nbsp;&nbsp;&nbsp;</p>\n";
		print dofooter($HDR['footer']);
		exit;
	}
	//	These have to do with job control submission
	//	Unfortunately the way we have our web server configured, the web server
	//	cannot create a file that topmed can read. So this is dead code and incomplete
    if ($fcn == 'setnodes') {
    	if ($PARMS['col'] == 'reset') {		// Delete slrum node list 
    		unlink($LDB['SLURMNODES']);
    		Msg("SLURM Node list deleted");
    	}
    	else {
    		$h = fopen($LDB['SLURMNODES'], 'w');
    		if ($h === false) {
    			fwrite($h, $PARMS['s']);
    			fclose($h);
    			Msg("SLURM Node list replaced");
    		}
    		else { Emsg("Failed to replace SLURM Node list"); }
	   	}
        print $JS['CLOSE'];
        print dofooter($HDR['footer']);
        exit;
    }
    if ($fcn == 'setcancel') {
        $sql = 'UPDATE ' . $LDB['samples'] . " SET " . $PARMS['col'] . "=89 WHERE bamid=" . $PARMS['id'];
        $result = SQL_Query($sql);
        Msg("Cancelled '$col' for BAM '$bamid'  SQL=$sql");
        print $JS['CLOSE'];
        print dofooter($HDR['footer']);
        exit;
    }
    if ($fcn == 'setrequest') {
        $sql = 'UPDATE ' . $LDB['samples'] . " SET " . $PARMS['col'] . "=1 WHERE bamid=" . $PARMS['id'];
        $result = SQL_Query($sql);
        Msg("Requested '$col' for BAM '$bamid'");
        print $JS['CLOSE'];
        print dofooter($HDR['footer']);
        exit;
    }
}

?>
