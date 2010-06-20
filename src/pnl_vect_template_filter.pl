#!/usr/bin/perl -w


## Jérôme Lelong 2008

## This filter is intended to be used by doxygen to expand the
## template files of PnlVect so that doxygen actually generates the
## whole documentation for the different atomic types.

use strict;

sub dump_file
{
    my ($file) = @_;
    my $FILE;
    if (open($FILE , "<", $file))
    {
        while (<$FILE>) { print $_;}
        close $FILE;
        return;
    }
    die "Unable to open $file";
}


sub parse
{
    my ($INPUT) = @_;
    my $in_comment = 0;
    while (my $line =  <$INPUT>)
    {
        
        if ($line =~ m/^#include[ ]*"(.*)"/ && ! $in_comment)
        {
            my $file = $1;
            dump_file ("include/" . $file) if ($file =~ m/pnl_templates/);
            dump_file ("linalg/" . $file) if ($file =~ m/_source\.c/);
        }
        else
        {
            print $line;
        }
        $in_comment = 1 if ($line =~ m@/\*@);
        $in_comment = 0 if ($line =~ m@\*/@);
    }
}



my $INPUT;
if (@ARGV >= 1)
{
    foreach my $f (@ARGV)
    {
        open($INPUT, "<", $f) or die "unable to open $f";
        parse ($INPUT);
        close $INPUT;
    }
}
else
{
    open ($INPUT, "-");
    parse ($INPUT);
    close $INPUT;
}



