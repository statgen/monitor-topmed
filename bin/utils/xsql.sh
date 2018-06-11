#!/bin/bash
#  xsql.sh sql-command-without-signle-quotes
#
#    Extract columns from a database - generating a file of just columnar output
#
q="'"
eval dosql.pl -batch -cmd ${q}$1${q} | awk '{if(NR>5)print}' | grep -v ' row' | grep -v Bye 
exit $?

