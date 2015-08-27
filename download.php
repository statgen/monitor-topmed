<?php
/*#################################################################
#
# Name: download.php
#
# Description:
#   Routines to generate HTML to download a file
#
# Copyright (C) 2014- University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/

/*---------------------------------------------------------------
# Download - Generate header to download a file or string
#
# Parameters:
#   fn - filename to show in header
#   string - string of data to be downloaded
#
# Returns:
#   Does not return if not errors found
#   Else returns an error message
---------------------------------------------------------------*/
function Download($fn='unknown', $string='') {

    if (! $string) { return "Download - no file or string provided"; }

    //  Normal mode, cause browser to download the data
    header("Content-Disposition: attachment; filename=$fn");
    header("Content-length: " . strlen($string));
    header("Content-type: application/download");
    header("Connection: close");
    header("Expires: 0");
    set_time_limit(0);
    print $string;
    exit;
}

/*---------------------------------------------------------------
# Uncompress - Possibly decompress a file
#
# Parameters:
#   infile  - name of input file
#   outfile - name of newly created output file
#
# Returns:
#   boolean if outfile was created
---------------------------------------------------------------*/
function Uncompress($infile, $outfile) {

    $s = `file -L $infile`;         // What kind of file is it

    if (preg_match('/Zip archive/', $s)) {          // ZIP
        $rc = system("unzip -p $infile > $outfile");
        if ($rc == 0) { return TRUE; }
        Nice_Exit("Unable to uncompress ZIP file '$infile'");
    }

    if (preg_match('/gzip compressed/', $s)) {      // GZIP
        $rc = system("gunzip -c $infile > $outfile");
        if ($rc == 0) { return TRUE; }
        Nice_Exit("Unable to uncompress GZIP file '$infile'");
    }

    if (preg_match('/bzip2 compressed/', $s)) {     // BZIP2
        $rc = system("bunzip2 -c $infile > $outfile");
        if ($rc == 0) { return TRUE; }
        Nice_Exit("Unable to uncompress BZIP2 file '$infile'");
    }

    return FALSE;
}

?>
