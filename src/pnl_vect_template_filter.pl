#!/usr/bin/perl -w


## Jérôme Lelong 2008

## This filter is intended to be used by doxygen to expand the
## template files of PremiaVect so that doxygen actually generates the
## whole documentation for the different atomic types.

use strict;

my @dirs = ( "../include", "include", ".", "linalg", "sort");

sub dump_file
{
    my ($file) = @_;
    my $FILE;
    foreach my $dir (@dirs)
        {
            if (open($FILE , "<", $dir . "/" . $file))
                {
                    while (<$FILE>) { print $_;}
                    close $FILE;
                    return;
                }
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
            dump_file ($file) if ($file =~ m/pnl_templates/ || $file =~ m/_source\.c/);
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
if (@ARGV == 1)
{
    open($INPUT, "<", $ARGV[0]) or die "unable to open $ARGV[0]";
}
else
{
    open ($INPUT, "-");
}

parse ($INPUT);
close $INPUT;


