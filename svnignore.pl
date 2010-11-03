#!/usr/bin/perl -w

## Jérôme Lelong, July 2009

## This script is intended to go through the Src tree and look for
## directories which do not have a svn:ignore property but which
## contains a Makefile.am. The file Makefile.am is used to detect
## whether we are in a true directory with code.

use strict;
use File::Find;
use File::Basename;
use Getopt::Long;

# for the convenience of &wanted calls, including -eval statements:
use vars qw/*name *dir *prune/;
*name  = *File::Find::name;
*dir   = *File::Find::dir;
*prune = *File::Find::prune;

my $dry_run = 0;
my $help    = 0;
my $verbose = 0;

GetOptions(
    "dry-run" => \$dry_run,
    "help"    => \$help,
    "verbose" => \$verbose
);

##
## Prints the manual page
##
sub print_help {
    print(
"usage : svnignore.pl [--help] [--dry-run] [--verbose]
\t--help : prints this help.
\t--dry-run : do not actually modify any property.
\t--verbose : print the details of the new propeties.
"
         );
}

sub update_ignore {
    if (/Makefile\.am\z/s) {
        my $dir                    = $File::Find::dir;
        my @ignore_list            = `svn propget svn:ignore $dir`;
        my $have_found_makefile    = 0;
        my $have_found_makefile_in = 0;
        my $have_found_tags         = 0;

        foreach (@ignore_list) {
            chomp if (m/^[ ]*\n$/);
            $have_found_makefile = 1
              if ( $have_found_makefile == 0 && m/Makefile[ ]*/ );
            $have_found_makefile_in = 1
              if ( $have_found_makefile_in == 0 && m/Makefile\.in/ );
              $have_found_tags = 1
              if ( $have_found_tags == 0 && m/tags/ );
            last
              if ( $have_found_makefile == 1 && $have_found_makefile_in == 1 
                    && $have_found_tags);
        }

        # svn propset svn:ignore expect a newline separated list.
        push( @ignore_list, "Makefile\n" )    unless ($have_found_makefile);
        push( @ignore_list, "Makefile.in\n" ) unless ($have_found_makefile_in);
        push( @ignore_list, "tags\n" ) unless ($have_found_tags);
        if ( $have_found_makefile == 0 || $have_found_makefile_in == 0
            || $have_found_tags == 0 ) {

            # Joining the list avoids addind useless leading spaces
            my $ignore_list = join( '', @ignore_list );
            if ( $dry_run == 1 ) {
                print "Setting svn:ignore on $dir\n";
            }
            else {
                system("svn propset svn:ignore \"$ignore_list\" $dir");
            }
            print( $ignore_list, "\n" ) if ( $verbose == 1 );
        }
    }
}

if ( $help == 1 ) {
    print_help();
    exit;
}

my %find_opt;

# this option is very important otherwise inside the function wanted,
# we are chdired to the dirname of the current file.
$find_opt{"no_chdir"} = 1;
$find_opt{"wanted"}   = \&update_ignore;
File::Find::find( \%find_opt, "src" );
