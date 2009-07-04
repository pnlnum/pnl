#!/usr/bin/perl -w

## Jérôme Lelong, July 2009
## Steps for adding a model
## 1. Create the directory toto and add it to svn
## 2. Create a template file toto.c and toto.h
## 3. Add a Makefile.am
## 4. At the upper-live, edit the Makefile to take into account
##    the new directory
## 5. Edit configure.in and add the new directory
## 6. Edit premia_obj.c
## 6. Add the documentation
## 7. svn add the different files

use strict;
use Getopt::Long;
use Cwd;

my $model = "";
my $option = "";
my $help    = 0;
my $svn    = 0;

GetOptions(
    "svn"    => \$svn,
    "help"    => \$help,
    "model=s" => \$model,
    "option=s" => \$option
    );

##
## Prints the manual page
##
sub print_help {
    print(
        "usage : svnignore.pl [--help] --model=mod [--option=opt] [--svn]
\t--model=mod : mod is the name of the model to be added.
\t--model=mod : opt is the name of the option to be added in the model mod.
\t--svn : add the created files to the svn repository. This option requires
--model and --option.
\t--help : prints this help.
");
}

##
## Creates the new directory for the model $mod
##
sub create_dir {
    my ($mod) = @_;
    mkdir $mod || die ("Cannot create $mod");
}

##
## Tests if a file alredy exists.
## Returns 1 if yes, 0 otherwise.
##
sub file_exists {
    my ($file) = @_;
    my $cwd = getcwd ();
    if (-f $file) {
        print "File $file already exists in $cwd. Skipping ...\n";
        return 1;
    }
    return 0;
}


##
## Dump a list (without newlines separator) to a file
## Args :
## 1. a file name
## 2. a list
sub dump_to_file{
    my ($filename, @list) = @_;
    my $cwd = getcwd ();
    open (FILE, ">$filename") || die ("Cannot create file $filename in $cwd.");
    print "Creating ${filename} in $cwd\n";
    foreach (@list){
        print FILE ($_, "\n");
    }
    close FILE;
}

##
## Creates $mod.c, $mod.h  in directory $mod.
## Before entering this function, you must alredy be in direcotry Src/mod.
## 
sub create_model_code_templates {
    my ($mod) = @_;
    my @content;
    chdir ($mod);
    # mod.h file
    chdir('..') && return if ( file_exists ("${mod}.h") );
    my $upper_mod = uc($mod);
    @content = ("#ifndef _${upper_mod}_H",
                "#define _${upper_mod}_H",
                "",
                "#include \"optype.h\"",
                "#include \"var.h\"",
                "",
                "#define TYPEMOD $upper_mod",
                "",
                "/* $upper_mod World */ ",
                "typedef struct TYPEMOD {",
                "  VAR T;",
                "  VAR S0;",
                "  VAR Divid;",
                "  VAR R; ",
                "  VAR Sigma;",
                "... Modify the parameters accordingly to the model ...",
                "} TYPEMOD;",
                "",
                "#endif",
                "");
    dump_to_file ("${mod}.h", @content);

    # mod.c file
    chdir('..') && return if ( file_exists ("${mod}.c") );
    @content = ("#include $mod" . ".c",
                "#include \"chk.h\"",
                "#include \"error_msg.h\"",
                "#include \"model.h\"",
                "static int MOD(Init)(Model *model)",
                "{",
                "  TYPEMOD* pt=(TYPEMOD*)(model->TypeModel);",
                "  static int first=1;",
                "",
                "  if (model->init == 0 )",
                "    {",
                "      model->init = 1;",
                "      model->nvar=0;",
                "      pt->T.Vname = \"Current Date\";",
                "      pt->T.Vtype=DATE;",
                "      pt->T.Val.V_DATE=0.;",
                "      pt->T.Viter=ALLOW;",
                "      model->nvar++;",
                "",
                "      pt->S0.Vname = \"Spot\";",
                "      pt->S0.Vtype=PDOUBLE;",
                "      pt->S0.Val.V_PDOUBLE=100.;",
                "      pt->S0.Viter=ALLOW;",
                "      model->nvar++;",
                "",
                "      pt->Divid.Vname = \"Annual Dividend Rate\";",
                "      pt->Divid.Vtype=DOUBLE;",
                "      pt->Divid.Val.V_DOUBLE=0.;",
                "      pt->Divid.Viter=ALLOW;",
                "      model->nvar++;",
                "",
                "      pt->R.Vname = \"Annual Interest Rate\";",
                "      pt->R.Vtype=DOUBLE;",
                "      pt->R.Val.V_DOUBLE=10.;",
                "      pt->R.Viter=ALLOW;",
                "      model->nvar++;",
                "",
                "      pt->Sigma.Vname = \"Volatility\";",
                "      pt->Sigma.Vtype=DOUBLE;",
                "      pt->Sigma.Val.V_DOUBLE=0.2;",
                "      pt->Sigma.Viter=ALLOW;",
                "      model->nvar++;",
                "",
                "      first=0;",
                "    }",
                "",
                "  return OK;",
                "}",
                "",
                "TYPEMOD ${mod};",
                "MAKEMOD(${mod});");
    dump_to_file ("${mod}.c", @content);
    chdir ('..');
}


##
## Creates $mod_opt.c, $mod_opt.h  in directory $mod.
## Before entering this function, you must alredy be in direcotry Src/mod/$mod.
## 
sub create_opt_code_templates {
    my ($model, $opt, $mod_opt) = @_;
    my @content;
    chdir ($mod_opt);
    # mod_opt.h file
    chdir('..') && return if ( file_exists ("${mod_opt}.h") );
    my $upper_mod_opt = uc($mod_opt);
    @content = ("#ifndef _${upper_mod_opt}_H",
                "#define _${upper_mod_opt}_H",
                "",
                "#include \"${model}/${model}.h\"",
                "#include \"${opt}/${opt}.h\"",
                "#include \"pnl_mathtools.h\"",                
                "",
                "#endif",
                "");
    dump_to_file ("${mod_opt}.h", @content);

    # mod_opt.c file
    chdir('..') && return if ( file_exists ("${mod_opt}.c") );
    @content = (
        "#include  \"${mod_opt}.h\"", 
        "",
        "int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)",
        "{",
        "  TYPEOPT* ptOpt=( TYPEOPT*)(Opt->TypeOpt);",
        "  TYPEMOD* ptMod=( TYPEMOD*)(Mod->TypeModel);",
        "  int status=OK;",
        "",
        "  return status;",
        "}",
        "",
        "/* extern PricingMethod MET(...); */",
        "",
        "PricingMethod* MOD_OPT(methods)[]={",
        "NULL",
        "};",
        "",
        "DynamicTest* MOD_OPT(tests)[]={",
        "NULL",
        "};",
        "",
        "Pricing MOD_OPT(pricing)={",
        "ID_MOD_OPT,",
        "MOD_OPT(methods),",
        "MOD_OPT(tests),",
        "MOD_OPT(ChkMix)",
        "};" );
    dump_to_file ("${mod_opt}.c", @content);
    chdir ('..');
}

##
## Edit Makefile.am
##
sub create_makefile{
    my ($mod) = @_;
    chdir ($mod);
    # Makefile.am
    chdir('..') && return if ( file_exists ("Makefile.am") );
    my @content = ("SUBDIRS = ",
                   'include $(top_srcdir)/Make.incl',
                   'AM_CPPFLAGS = $(PNL_INCLUDES) $(PREMIA_MINIMUM_INCLUDES)',
                   "noinst_LTLIBRARIES = lib${mod}.la",
                   "lib${mod}_la_SOURCES = ${mod}.c",
                   "lib${mod}_la_LIBADD = ");
    dump_to_file ("Makefile.am", @content);
    chdir ('..');
}


##
## In Makefiles, it is very common to come accross lines ending by a
## backslash to avoid too long lines. This function searches for the
## last such line and contenate some extra data
## Args :
##     $line : a reference to the current line
##     $MAKE : the filehandle through which the Makefile is read
##     $modified_content : a ref to the array in which we store the modified content
##     $to_add : a string to be concatened to the original data
##
sub jump_to_eol {
    my ($line, $MAKE, $modified_content, $to_add) = @_;
    $_ = $$line;
    if (m/\\$/) {
        push (@$modified_content, $_);
        while (<$MAKE>){
            last if (! m/\\$/);
            push (@$modified_content, $_);
        }
    }
    chomp;
    if (length ($_) + length ($to_add) > 72 && !m/$to_add/){
        push (@$modified_content, "$_ \\\n");
        $_  = "$to_add \n";
    }
    elsif (! m/$to_add/) {
        $_  .= " $to_add \n" ;
    }
    else{
        $_ .= "\n";
    }
    $line = $_;
}

##
## Edit the Makefile.am at the upper level.
## Arg : $mod is the name of the model.
##
sub edit_upper_makefile {
    my ($mod) = @_;
    my $cwd = getcwd ();
    # Edit the upper Makefile.am
    # Read file
    open (MAKE, "<Makefile.am")  || die ("Cannot open file Makefile.am");
    print "Editing Makefile.am in ${cwd} \n";
    # Modify file
    my @modified_content;
    while (<MAKE>){
        if (m/SUBDIRS/){
            jump_to_eol (\$_, *MAKE, \@modified_content, $mod);
        }
        if (m/la_LIBADD/) {
            jump_to_eol (\$_, *MAKE, \@modified_content, "$mod/lib${mod}.la");
        }
        push (@modified_content, $_);
    }
    close (MAKE);
    # Write modification to original file
    open (MAKE, ">Makefile.am")  || die ("Cannot open file Makefile.am");
    print MAKE @modified_content;
    close (MAKE);
}

sub edit_configure {
    my ($mod) = @_;
    my $make = "Src/mod/${mod}/Makefile";
    # Edit the configure.in
    # Read file
    open (CONFIGURE, "<configure.in")  || die ("Cannot open file configure.in");
    print "Editing configure.in\n";
    # Insert $make in the argument list of AC_OUTPUT and make sure the
    # argument list remains sorted.
    my @modified_content;
    while (<CONFIGURE>){
        if (m/AC_OUTPUT/){
            push (@modified_content, $_);
            while (<CONFIGURE>) {
                if (m/\\$/ && $make gt $_) {
                    push (@modified_content, $_);
                }
                else {
                    last;
                }
            }
            $_ = "$make \\\n" .  $_  unless (m/$make/);
        }
        push (@modified_content, $_);
    }
    close (CONFIGURE);
    open (CONFIGURE, ">configure.in")  || die ("Cannot open file configure.in");
    print CONFIGURE @modified_content;
    close (CONFIGURE);    
}

##
## Creates the directory man/tex/mod/ add a template $mod/$mod.h file
##
sub create_mod_doc_templates {
    my ($mod) = @_;
    my @content;
    chdir ($mod);
    # mod_doc.tex file
    chdir('..') && return if ( file_exists ("${mod}_doc.tex") );
    @content = ( "\\documentclass[12pt,a4paper]{article}",
                 "\\input premiamble",
                 "\\input premiadata",
                 "\\begin{document}",
                 "\\input{${mod}_docl.tex}",
                 "",
                 "\\section{Description}",
                 "",
                 "",
                 "\\section{Code Implementation}",
                 "\\verbatiminput{../../../../Src/mod/${mod}/${mod}.h",
                 "\\input premiaend",
                 "\\end{document}",
                 "");
    dump_to_file ("${mod}_doc.tex", @content);
    chdir ('..');
}

##
## Creates the directory man/tex/mod/$mod add a template $mod_opt/$mod_opt.h file
##
sub create_opt_doc_templates {
    my ($mod, $opt, $mod_opt) = @_;
    my @content;
    chdir ($mod_opt);
    # mod_doc.tex file
    chdir('..') && return if ( file_exists ("${mod_opt}_doc.tex") );
    @content = ( "\\documentclass[12pt,a4paper]{article}",
                 "\\input premiamble",
                 "\\input premiadata",
                 "\\begin{document}",
                 "\\input{${mod}_docl.tex}",
                 "",
                 "\\section{Description}",
                 "",
                 "Insert a short description of the model",
                 "",
                 "\\section{Code Implementation}",
                 "\\verbatiminput{../../../../Src/mod/${mod}/${mod}.h",
                 "\\input premiaend",
                 "\\end{document}",
                 "");
    dump_to_file ("${mod_opt}_doc.tex", @content);
    chdir ('..');
}


if ($model eq "") {
    print("You must supply a model name\n");
    print_help ();
    exit ;
}


if ($svn == 1) {
    if ($model eq "" || $option eq "") {
        print("You must supply both model and option names\n");
        print_help ();
        exit ;
    }
    system ("svn add Src/mod/${model}");
    system ("svn add man/tex/mod/${model}");
    system ("svn propset svn:ignore \"Makefile\nMakefile.am\" Src/mod/${model}");
    system ("svn propset svn:ignore \"Makefile\nMakefile.am\" Src/mod/${model}/${model}_${option}");
    exit 0;
}

#Code
chdir ("Src"); chdir ("mod");
create_dir ($model);
create_model_code_templates ($model);
create_makefile ($model);
edit_upper_makefile ($model);
chdir("../..");
edit_configure ($model);
print "\nYou still have to add your model to Src/premia_obj.c
\textern Model " . uc(${model}) . "_model;
\t&" . uc(${model}) . "_model, to the right models array.
\n";

#Doc
chdir ("man/tex/mod");
create_dir ($model);
create_mod_doc_templates ($model);
chdir ('../../..');

if ($option ne "") {
    my $mod_opt = join ('_', ($model, $option));
#Code
    print "\n\n";
    chdir ("Src/mod");
    chdir ($model);
    create_dir ($mod_opt);
    create_opt_code_templates ($model, $option, $mod_opt);
    create_makefile ($mod_opt);
    edit_upper_makefile ($mod_opt);
    chdir ('../../..');    
    edit_configure ("${model}/${mod_opt}");
    print "\nYou still have to add your option to Src/premia_obj.c
\textern Pricing " .  uc(${mod_opt}) . "_pricing;
\t&" . uc(${mod_opt}) . "_pricing, to the right pricings array.
\n";
#Doc
    chdir ("man/tex/mod/${model}");
    create_dir ($mod_opt);
    create_opt_doc_templates ($model, $option, $mod_opt);
    chdir ('../../..');
}

print "\nIf you are happy the created files, rerun with --svn.\n
";
