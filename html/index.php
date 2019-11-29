<?php
/*#################################################################R
#
# Name: index.php
#
# Description:
#   Display information for tracking NHLBI TopMed data
#
# Copyright (C) 2015-2019 Terry Gliedt, University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/
include_once 'local_config.php';
include_once 'header.php';

print doheader($HDR['title'], 1);

print "<h4>Use these links to look at NHLBI TopMed data</h4>" .
	"<ul>\n" .
	"<li><a href='genome.php' target='genome'>Genome Sequence data pages</a></br></li>\n" .
	"<li><a href='rnaseq.php' target='rnaseq'>RNA Sequence data pages</a></br></li>\n" .
	//"<li><a href='origindex.php' target='genome'>Original Genome data pages</a></br></li>\n" .
	"</u>\n";
	
print dofooter($HDR['footer']);
exit;

?>
