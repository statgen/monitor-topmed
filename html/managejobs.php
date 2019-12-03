<?php
/*#################################################################
#
# Name: managejobs.php
#
# Description:
#   Code to support SLURM management
#
# Copyright (C) 2019 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/

/*---------------------------------------------------------------
#   JOBFunctions($fcn, $bamid, $samplestate, $op)
#   Handle functions related to RNA processing. Might not return
---------------------------------------------------------------*/
function JOBFunctions($fcn) {
    global $HDR, $LDB, $GLOBS, $PARMS;

	if ($fcn == 'queue') {
		print "<center>" . $GLOBS['showstatus'] . "&nbsp;&nbsp;&nbsp;</center>\n";
		print "<font size='-1'><table border='0' align='center' width='80%'><tr>" .
			"<td align='left'>&nbsp;&nbsp;&nbsp;<a href='javascript:window.location.reload()'>Reload</a></td>" .
			"<td align='right'><a href='javascript:window.close()'>Close</a>&nbsp;&nbsp;&nbsp;</td>" .
			"</tr></table></font><br/>\n";
		$c = $LDB['bindir'] . "/topmedcluster.pl squeue 2>&1";
		print "<pre>\n" . shell_exec($c) . "</pre>\n";
		exit;
	}
	if ($fcn == 'showout') {                // Show output from a SLURM job
		$s='no samplestate provided';
		$samplestate = $PARMS['samplestate'];
		$id = $PARMS['id'];
		if ($samplestate == '5') { $s = 'verify'; }
		if ($samplestate == 'B') { $s = 'backup'; }
		if ($samplestate == 'Q') { $s = 'qplot'; }
		if ($samplestate == 'C') { $s = 'cram'; }
		if ($samplestate == 's') { $s = 'gcepush'; }
		if ($samplestate == 'r') { $s = 'gcepull'; }
		if ($samplestate == 'V') { $s = 'bcf'; }
		if ($samplestate == 'A') { $s = 'awscopy'; }
		if ($samplestate == 'G') { $s = 'gcecopy'; }
		if ($samplestate == 'g') { $s = 'gcecpbcf'; }
		if ($samplestate == 'x') { $s = 'cleanup'; }
		// Hardcoded path cause mario won't play with topmedpath.pl
		$f = '/net/topmed/working/topmed-output/' . $PARMS['id'] . "-$s.out";
		print "<h4 align='center'>SLURM console log for '$s' Sample=$id [" . $PARMS['table'] . "]</h4>\n" .
		"<pre>\nFile='$f'\n";
		if (! file_exists($f)) { print "No console file found for '$id'\n"; }
		else {					// Backticks are deprecated
			$h = popen("/bin/cat $f 2>&1", 'r');
			$read = fgets($h, 65536);
			print $read;
			pclose($h);
		}
		print "</pre>\n";
		print dofooter($HDR['footer']);
		exit;
	}
	if ($fcn == 'permit') {
		$h = HandlePermit($PARMS['op'],$PARMS['datayear'],$PARMS['center'],$PARMS['run'],$PARMS['id']);
		print ControlJobs($h);
		print dofooter($HDR['footer']);
		exit;
	}
	if ($fcn == 'control') {
		print ControlJobs('');
		print dofooter($HDR['footer']);
		exit;
	}
	if ($fcn == 'restart') {
		print RestartJobs('');
		print dofooter($HDR['footer']);
		exit;
	}
	if ($fcn == 'restartjobs') {
		$h = HandleRestartJobs($PARMS['run'],$PARMS['samplestate'],$PARMS['op']);
		if ($h != "") {
			print $h;
			print RestartJobs('');
		}
		print dofooter($HDR['footer']);
		exit;
	}

	if ($fcn == 'logs') {
		$p = getenv('PROJECT');
		print "<center>" . $GLOBS['links'] . " &nbsp;&nbsp;&nbsp;</center>\n<pre>\n";
		$d = "/net/$p/working/$p-output/";
		if (! chdir($d)) { print "Cannot CD to '$d': $!\n"; }
		else {
			$logs=explode("\n", `ls $p*.log`);
			foreach ($logs as $f) {
				if ($f) { print "<b>Showing $f</b>\n" . `tail -6 $f`; }
			}
		}
		print "</pre>\n";
		exit;
	}
}

/*---------------------------------------------------------------
#   html = HandlePermit($op, $center, $run, $id)
#   Update permissions table, return HTML about results
---------------------------------------------------------------*/
function HandlePermit($op, $datayear, $center, $run, $id) {
    global $GLOBS, $LDB;
    $html = "Nothing honey";
    if ($op == '') { return ''; }       // Avoid misleading err msgs

    //  Delete a permission
    if ($op == 'del') {
        $cmd = escapeshellcmd("/usr/cluster/topmed/bin/topmedpermit.pl remove $id");
        $s = `$cmd 2>&1`;
        return "<pre>$s</pre>\n";
        //return "<pre>cmd=$cmd\n$s</pre>\n";
    }
    
    //  Set a permission
    if (! in_array($op, $GLOBS['validfunctions'])) {
        return Emsg("Invalid operation '$op' - How'd you do that?", 1);
    }
    $cmd = escapeshellcmd("/usr/cluster/topmed/bin/topmedpermit.pl add $op $datayear $center $run");
    $s = `$cmd 2>&1`;
    return "<pre>$s</pre>\n";
}

/*---------------------------------------------------------------
#   html = ControlJobs($h)
#   Show details for a particular run
---------------------------------------------------------------*/
function ControlJobs($h) {
    global $LDB, $GLOBS;
    $url = $_SERVER['PHP_SELF'] . "?fcn=permit";    
    $html = '';

    //  Dump table for each permission for now
    $html .= "<h3 align='center'>Control Job Submission</h3>\n";
    $html .= "<center>These controls do not affect jobs already submitted<br>" .
        $GLOBS['showstatus'] . "</center><br/>\n";
    if ($h) { $html .= "<div class='indent'><span class='surprise'>$h</span></div>\n"; }

    //  Prompt for fields which can control permissions  
    $html .= "<h4>Use this to stop future job submissions</h4>\n" .
    	"<form action='$url' method='post'>\n" .
        "<input type='hidden' name='fcn' value='permit'>\n" .
        "<table align='center' width='80%' border='0'>\n" .

        "<tr>" .
        "<td align='right'><b>Data Year</b></td>" .
        "<td>&nbsp;</td>" .
        "<td><input type='text' name='datayear' length='8' value='all'></td>" .
        "<td>&nbsp;</td>" .
        "<td><font color='green'> 1, 2 3 etc.  </font></td>" .
        "</tr>\n" .

        "<tr>" .
        "<td align='right'><b>Center</b></td>" .
        "<td>&nbsp;</td>" .
        "<td><input type='text' name='center' length='8' value='all'></td>" .
        "<td>&nbsp;</td>" .
        "<td><font color='green'> broad, nygc, uw etc.  </font></td>" .
        "</tr>\n" .

        "<tr>" .
        "<td align='right'><b>Name of Run</b></td>" .
        "<td>&nbsp;</td>" .
        "<td><input type='text' name='run' length='16' value='all'></td>" .
        "<td>&nbsp;</td>" .
        "<td><font color='green'> e.g. 2015aug22.weiss.06 etc. </font></td>" .
        "</tr>\n" .

        "<tr>" .
        "<td align='right'><b>Operation</b></td>" .
        "<td>&nbsp;</td>";
    $html .= "<td><select name='op'>";
    foreach ($GLOBS['validfunctions'] as $fcn) { $html .= "<option value='$fcn'>$fcn</option>"; }
    $html .= "</select></td>" .
        "<td><font color='green'>&nbsp; </font></td>" .
        "<td>&nbsp;</td>" .
        "</tr>\n" .
        "<tr>" .
        "<td><input type='submit' value=' Stop Submissions '></td>" .
        "<td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td>" .
        "</tr>\n" .
        "</table>\n</form>\n";
    
    //  This shows what is in effect and allows on to delete an entry
    $hdrcols  = array('operation', 'datayear', 'centername', 'dirname');
    $sql = 'SELECT * FROM ' . $LDB['permissions'];
    $result = SQL_Query($sql);
    $numrows = SQL_NumRows($result);

	//	Unfortunately the way we have our web server configured, the web server
	//	cannot create a file that topmed can read. So this is dead code
	//	Give instructions how to manually change where jobs are queued
	$html .= <<<END
	<br/><h4>Control Where Submitted Jobs Will Run</h4>
	<table align='center' width='80%' border='0'><tr><td>
	<p>You can override which nodes are used when jobs are submitted by creating
	a NINE line file of nodes to be used in <b>/net/topmed/working/topmed-output/slurm.nodes</b>
	This must consist of a simple node list based on the sampleid modulo 9.
	If this is messed up, submitting jobs will fail. 
	This only affects jobs to be submitted, not those already submitted. It does not
	cancel jobs nor change those submitted.
	<font color='red'><b>Remove this when it is no longer needed.</b></font>
	</p>
	<p>You can see jobs queued with this command (as user topmed):<br/>
	<code>  squeue -o '%.18i %.14j %.8u %.2t %.10M %R' -q topmed-qplot</code>
	The last field is the QOS. See all the QOS with <code>sacctmgr -P show qos</code>
	</p>
	<p>You can cancel jobs by using <code>scancel JOBID</code>. The jobid is the first
	field in the 'squeue -q' output. Grep the output to get the jobids you want.
	</p>
	</p><font color='blue'>After canceling jobs, you'll have to reset the state in
	the monitor so they will be resubmitted.</font>. This can be done by using
	<code>topmedcmd.pl set BAMID state_NAME 0</code>. BAMID is part of the SLURM jobname.
	Putting this all together is left as an exercise for the student.
	</p>
	</td></tr></table>
END;

    $html .= "<br/>\n" .
        "<h4>These job submission controls are in effect</h4>\n" .
        "<p>Use '<b>Remove</b>' above to delete a control and allow job submissions to resume.</p>\n" .
        "<table align='center' width='90%' border='1'>\n" .
        "<tr>\n";
    foreach ($hdrcols as $c) {
        $html .= "<th class='heading'>" . ucfirst($c) . "</th>\n";
    }
    $html .= "<th class='heading'>&nbsp;</th></tr>\n" .
        "<tr>";

    for ($i=0; $i<$numrows; $i++) {
        $row = SQL_Fetch($result);
        reset($hdrcols);
        foreach ($hdrcols as $c) {
            $d = $row[$c];
            if ((! isset($d)) || ($d == '') || $d == '0') { $d = 'all'; }
            $html .= "<td align='center'>$d</td>\n";
        }
        $html .= "<td align='center'><a href='$url&amp;op=del&amp;id=$row[id]'>" .
            "<font color='green'>Remove</font></a></td>" .
            "</tr>\n";
    }
    $html .= "</table>\n";
    return $html;
}

/*---------------------------------------------------------------
#   html = RestartJobs($h)
#   Restart failed, running, or queued jobs
---------------------------------------------------------------*/
function RestartJobs($h) {
    global $LDB, $GLOBS, $PARMS;
    $url = $_SERVER['PHP_SELF'] . "?fcn=restart";    
    $html = '';

    //  Dump table for each permission for now
    $html .= "<h3 align='center'>Restart Failed Jobs</h3>\n" .
        "<h4 align='center'>Use this to restart failed jobs for a run. Beware this affects <u>all</u><br>" .
        "samples in the selected state for the specified run.</h4>\n";
            $html .= $GLOBS['showstatus'] . "</center><br/>\n";
    if ($h) { $html .= "<div class='indent'><span class='surprise'>$h</span></div>\n"; }

    //  Prompt for classes of jobs to be restarted
    $html .= "<form action='$url' method='post'>\n" .
        "<input type='hidden' name='fcn' value='restartjobs'>\n" .
        "<table align='left' width='80%' border='0'>\n" .
        "<tr>" .
        "<td align='right'><b>Name or Runid of Run</b></td>" .
        "<td>&nbsp;</td>" .
        "<td><input type='text' name='run' length='16' value='" . $PARMS['run'] . "'></td>" .
        "<td>&nbsp;</td>" .
        "<td><font color='green'> 2015aug22.weiss.06, 244, etc. </font></td>" .
        "</tr>\n" .

        "<tr>" .
        "<td align='right'><b>Current State of Sample</b></td>" .
        "<td>&nbsp;</td>";
    $html .= "<td><select name='samplestate'>" .
        "<option value='99'>Failed</option>" .
        "<option value='2'>Submitted</option>" .
        "<option value='3'>Started</option>" .
        "<option value='19'>Delivered</option>" .
        "<option value='98'>FailedChecksum</option>" .
        "</select></td>" .
        "<td><font color='green'>&nbsp;</font></td>" .
        "<td>&nbsp;</td>" .
        "</tr>\n" .

        "<tr>" .
        "<td align='right'><b>Operation</b></td>" .
        "<td>&nbsp;</td>";
    $html .= "<td><select name='op'>" .
        "<option value='none'>none</option>" .
        "<option value='arrive'>arrive</option>" .
        "<option value='verify'>verify</option>" .
        "<option value='cram'>cram</option>" .
        "<option value='backup'>backup</option>" .
        "<option value='qplot'>qplot</option>" .
        "<option value='gce38push'>gcepush</option>" .
        "<option value='gce38pull'>gcepull</option>" .
        "<option value='gce38post'>post</option>" .
        "<option value='gce38bcf'>bcf</option>" .
        "<option value='gce38copy'>gcecopy</option>" .
        "<option value='gce38cpbcf'>gcecpbcf</option>" .
        "<option value='gcecleanup'>gcecleanup</option>" .
        "<option value='aws38copy'>awscopy</option>" .
        "</select></td>" .
        "<td><font color='green'>&nbsp; </font></td>" .
        "<td>&nbsp;</td>" .
        "</tr>\n" .

        "<tr>" .
        "<td><input type='submit' value=' Restart Operation '></td>" .
        "<td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td>" .
        "</tr>\n" .
        "</table>\n</form>\n";
    
    return $html;
}

/*---------------------------------------------------------------
#   html = HandleRestartJobs($run, $samplestate, $op)
#   Update database for groups of samples in state to new state
---------------------------------------------------------------*/
function HandleRestartJobs($dirname, $samplestate, $op) {
    global $LDB, $GLOBS, $PARMS;
    if ($dirname == 'none' || $op == 'none') {
        return Emsg("Run or Operation was not specified, try again", 1);
    }
	$runspkey=$LDB['runs_pkey'];
	$runstable=$LDB['runs'];
	$samplespkey=$LDB['samples_pkey'];
	$samplestable=$LDB['samples'];
	if ($LDB['datatype'] == 'rnaseq') {
		$runspkey=$LDB['projects_pkey'];
		$runstable=$LDB['projects'];
	}


    if (preg_match('/^\d+$/', $dirname)) { $col = $runspkey; }
    else { $col = 'dirname'; }
    $sql = "SELECT $runspkey FROM $runstable WHERE $col='$dirname'";
    $result = SQL_Query($sql);
    $row = SQL_Fetch($result);
    $runsid = $row[$runspkey];
    if (! $runsid) {
        return Emsg("Run '$dirname' is not known, try again", 1);
    }

    $sql = "UPDATE $samplestable SET ";
    $sql .= "state_$op=0  WHERE $runspkey=$runsid AND state_$op=$samplestate";
    $html = "<h3>Changing State for Samples in '$dirname'</h3>\n" .
        "SQL=$sql<br>\n";
    $result = SQL_Query($sql);
    $changes = SQL_AffectedRows();

    $html .= Emsg("Changed $changes rows, hope that was what you wanted", 1);
    $html .= "<br/><br/>\n";
    return $html;
}

?>
