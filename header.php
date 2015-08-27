<?php
/*#################################################################
#
# Name: header.php
#
# Description:
#   Routines to generate headers and footers
#
# Copyright (C) 2014- University of Michigan
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; See http://www.gnu.org/copyleft/gpl.html
#################################################################*/

/*---------------------------------------------------------------
# doheader - Generate HTML headers
#   Depends on $HDR in local_config.php to set CSS, JS and
#   META tags in the header
#
# Parameters:
#   title - subportion of title
#   banner - optional boolean if banner should be generated
#
# Returns:
#   string of HTML
---------------------------------------------------------------*/
function doheader($title='No Title', $banner=0) {
    global $HDR;

    if ($HDR['title']) { $title = $HDR['title']; }
    header("Content-Type: text/html;charset=iso-8859-1");
    $expires = 60*60*24*1;                  // seconds, minutes, hours, days
    $expires = 30;
    header("Pragma: public");
    header("Cache-Control: maxage=".$expires);
    header('Expires: ' . gmdate('D, d M Y H:i:s', time()+$expires) . ' GMT');
    $m = $HDR['meta'];
    $h = <<<END
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
$m
    <title>$title</title>

END;

    if (isset($HDR['local'])  && $HDR['local']) { $h .= $HDR['local']; }
    $h .= "</head>\n<body>\n<div class='container'>\n";

    if ($banner) {
        $u = $HDR['home'];
        $l = $HDR['logo'];
        $t = $HDR['title'];
        $h .= <<<END
<table class="hdr" width="90%" align="center"><tr>
<td align='center'><a href="$u" class="hdr"><img src="$l" alt="logo" border='0'></a></td>
<td class='hdr2'>$t</td>
</tr></table><br/>

END;
    }
    return $h . "\n\n";
}

/*---------------------------------------------------------------
# dofooter - Generate HTML footer
#
# Parameters:
#   text - HTML centered text in footer
#   class - CSS class for DIV for footer
#
# Returns:
#   string of HTML
---------------------------------------------------------------*/
function dofooter($text='&nbsp;', $class='footer') {

    $h = <<<END
<br/><br/><br/><br/><br/><br/><br/><br/>
<div class='$class'>
<hr size='4' width='80%' noshade='noshade'>
<center>
$text
</center>
</div>
</body></html>

END;
    return $h;
}

?>
