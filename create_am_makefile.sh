#!/bin/bash

#
# This script creates a Makefile.am in the current 
# dir using all the .c files
#


gener ()
{


    files=`ls`

    sub_dir=''
    sources=''
    libs=''

    for f in $files
    do
        if [[ -f $f && `expr $f : '.*\.[c]$'` > 0 ]]
        then
            sources="$sources $f"
        fi
    done


    cur_dir=`pwd`
    cur_dir=`basename $cur_dir`
    lib="lib${cur_dir}_la"

    echo "include \$(top_srcdir)/Make.incl" >> Makefile.am
    echo "AM_CPPFLAGS=\$(PNL_INCLUDES) " >> Makefile.am
    echo "" >> Makefile.am
    echo "noinst_LTLIBRARIES=lib${cur_dir}.la" >> Makefile.am
    echo "${lib}_SOURCES=$sources" >> Makefile.am
    echo "${lib}_LDFLAGS=-avoid-version" >> Makefile.am

}

# dirs=`find . -type d ! -name 'CVS'`
# for d in $dirs
#   do
#   cwd=`pwd`
#   cd $d
#   rm -f Makefile.am
#   gener
#   cd $cwd
# done

rm -f Makefile.am
gener
