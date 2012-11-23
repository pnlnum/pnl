dnl Jérôme LELONG Feb 2006


AC_INIT(PNL, 1.3, jerome.lelong@imag.fr)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADERS(src/include/pnl/pnl_config.h:config.h.in)
AC_CONFIG_MACRO_DIR([config])
AM_INIT_AUTOMAKE()

if test "X${CXXFLAGS}" != "X"; then
   user_set_cxxflags=yes
fi
if test "X${CFLAGS}" != "X"; then
   user_set_cflags=yes
fi



AM_PROG_LIBTOOL

## Initialize maintainer mode
AM_MAINTAINER_MODE


# dnl Get system type
# dnl use "config.guess" and "config.sub" from directory "config"
# AC_CANONICAL_HOST


AC_ARG_ENABLE(debug,
    [  --enable-debug          compile with "-g" option, default behaviour is to optimize],
    [enable_debug="yes"], [enable_debug="no"])
AC_ARG_ENABLE(optimize,
    [  --disable-optimize      disable "-0" option off gcc], [enable_optimize="no"], [enable_optimize="yes"])
AC_ARG_ENABLE(verbose,
    [  --disable-silent         tell libtool not to be silent], [], [enable_verbose="no"])
AC_ARG_WITH(internal_blas,
    [  --with-internal-blas    use the internal Blas and Lapack libraries ])

AC_ARG_WITH(external_blas,
    [  --with-external-blas    Do not compile the internal Blas Lapack F2C])

AC_ARG_ENABLE(mpi,
    [  --disable-mpi      do not check for MPI], [enable_mpi="no"], [enable_mpi="yes"])


if test "$enable_verbose" = no; then
   LIBTOOL="$LIBTOOL --silent"
fi

################
## build Paths
################

local_dir=`pwd`
local_home=$local_dir
AC_PREFIX_DEFAULT($local_home)
if test "x$prefix" = "xNONE" ; then
   prefix=$local_home
fi

libdir="$prefix/lib"
includedir="$prefix/include"
LIBPNLDIR=`pwd`
##############################
## define C and C++ compilers
##############################


AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_F77
AC_PROG_MAKE_SET
AC_PROG_INSTALL

RM="rm -f"
##########################
## updating flags
##########################

if test "x$enable_optimize" = xyes  && test "x$enable_debug" = xno ; then
   AC_DEFINE(PNL_RANGE_CHECK_OFF, 1, [turn off range checking by default internally])
fi

# removing -ansi from CFLAGS
CC_OPTIONS="`echo $CFLAGS | sed -e 's/-ansi//g'`"
CXX_OPTIONS=$CXXFLAGS
if test "x$GCC" = xyes ; then
   if test "x$enable_debug" = xyes || test "x$enable_optimize" = xno ; then
      CC_OPTIONS=`echo $CC_OPTIONS | sed -e 's/-O[[0-9]]*//g'`
      CC_OPTIONS=`echo $CC_OPTIONS | sed -e 's/-g[[0-9]]*$//g'`
      CC_OPTIONS=`echo $CC_OPTIONS | sed -e 's/-g[[0-9]]* //g'`
      CC_OPTIONS="${CC_OPTIONS} -Wall -g"
      CXX_OPTIONS=`echo $CXX_OPTIONS | sed -e 's/-O[[0-9]]*//g'`
      CXX_OPTIONS=`echo $CXX_OPTIONS | sed -e 's/-g[[0-9]]*$//g'`
      CXX_OPTIONS=`echo $CXX_OPTIONS | sed -e 's/-g[[0-9]]* //g'`
      CXX_OPTIONS="${CXX_OPTIONS} -Wall -g -fno-inline"
   else
      CC_OPTIONS="${CC_OPTIONS} -Wall"
   fi
   #The -mno-cygwin flag to gcc makes gcc link against standard Microsoft DLLs instead of Cygwin.
   #This is desirable for native Windows programs that don't need a UNIX emulation layer.
   case "$host" in
        *-*-cygwin*)
	 CC_OPTIONS="${CC_OPTIONS} -mno-cygwin"
	 ;;
   esac
   #This option enables compilation of the free or the complete versions of the source code
   #for further explanation cf the directory Unix\FreePremiaScripts\
   if test "$with_premiaversion" != "" ; then
      CC_OPTIONS="${CC_OPTIONS} -DPremiaCurrentVersion=${PREMIAVERSION}"
   fi
fi
AM_CFLAGS=${CC_OPTIONS}
AM_CXXFLAGS=${CXX_OPTIONS}

# if the user has not changed the flags, reset them.
# They are passed to AM_xxxx
if test X$user_set_cxxflags != Xyes; then
   CFLAGS=""
fi
if test X$user_set_cflags != Xyes; then
   CXXFLAGS=""
fi

LDFLAGS="${LDFLAGS} -lm"

AC_C_INLINE
dnl Check for "extern inline", using a modified version of the test
dnl for AC_C_INLINE from acspecific.mt
dnl
PNL_HAVE_INLINE=
PNL_INLINE_DECL=
PNL_INLINE_FUNC=
AC_CACHE_CHECK([for extern inline], ac_cv_c_extern_inline,
                [ac_cv_c_extern_inline=no
                AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[extern $ac_cv_c_inline double foo(double x);
                extern $ac_cv_c_inline double foo(double x) { return x + 1.0 ; } ;
                double foo (double x) { return x + 1.0 ; };]],
                [[  foo(1.0)  ]])],[ac_cv_c_extern_inline="yes"],[])])

if test "$ac_cv_c_extern_inline" != no ; then
    PNL_HAVE_INLINE=1
    PNL_INLINE_DECL="extern"
    PNL_INLINE_FUNC="extern $ac_cv_c_inline"
    AM_CFLAGS="${AM_CFLAGS} -DPNL_HAVE_INLINE"
else
    AC_CACHE_CHECK([for C99 inline], ac_cv_c_99_inline,
                    [ac_cv_c_99_inline=no
                    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[extern inline void* foo() { foo(); return &foo ; };]], 
                    [[  return foo() != 0 ]])],[ac_cv_c_c99_inline="yes"],[])])
    if test "$ac_cv_c_99_inline" == "yes" ; then
        PNL_INLINE_DECL="$ac_cv_c_inline"
        PNL_INLINE_FUNC="$ac_cv_c_inline"
        PNL_HAVE_INLINE=1
    fi
fi
AC_SUBST(PNL_HAVE_INLINE)
AC_DEFINE_UNQUOTED(PNL_HAVE_INLINE, $PNL_HAVE_INLINE, [Define if you have inline])
AC_DEFINE_UNQUOTED(PNL_INLINE_DECL, $PNL_INLINE_DECL, [Keyword for declaring inline functions])
AC_DEFINE_UNQUOTED(PNL_INLINE_FUNC, $PNL_INLINE_FUNC, [Keyword for defining inline functions])




##########################
## check for C functions
##########################
AC_LANG(C)
AC_CHECK_LIB(m,main)
AC_CHECK_FUNCS(tgamma)
AC_CHECK_FUNCS(lgamma)
AC_CHECK_FUNCS(exp10)
AC_CHECK_FUNCS(trunc)
AC_CHECK_FUNCS(isfinite)
AC_CHECK_FUNCS(finite)
AC_CHECK_FUNCS(isnan)
AC_CHECK_FUNCS(isinf)

########################################################
## check for blas
########################################################
PNL_HAVE_BLAS=no
PNL_HAVE_LAPACK=no
BLAS_LIBS=""
LAPACK_LIBS=""
use_internal_blas=0

case "$host" in
    # Darwin/MacOSX
    *-*-darwin*)
    if test "x$with_internal_blas" = "xyes"; then
        use_internal_blas=1
    else
        # When Xcode is installed, Atlas and Lapack can be found in vectlib
        BLAS_LIBS="-lblas"
        LAPACK_LIBS="-llapack"
        ac_cv_f77_mangling="lower, underscore"
    fi
    ;;

    *)
    if test "x$F77" = "x"; then
       with_internal_blas=yes
    fi

    if test "x$with_internal_blas" != "xyes"; then
       builtin(include,config/acx_blas.m4)
       ACX_BLAS([PNL_HAVE_BLAS=yes])

       if test "$PNL_HAVE_BLAS" = yes; then
         AM_CFLAGS="$AM_CFLAGS"
       else
         with_internal_blas=yes
       fi

       builtin(include,config/acx_lapack.m4)
       ACX_LAPACK([PNL_HAVE_LAPACK=yes])

       if test "$PNL_HAVE_LAPACK" = yes; then
         AM_CFLAGS="$AM_CFLAGS"
       else
         with_internal_blas=yes
       fi
    fi

    if test "x$with_internal_blas" = "xyes"; then
        use_internal_blas=1
    fi

    if test "x$BLAS_STATIC_LIBS" != "x" || test "x$LAPACK_STATIC_LIBS" != "x"; then
       STATIC_FLAG="-static"
    fi
    ;;
esac

if test "x$with_external_blas" = "xyes"; then
    use_internal_blas=0
fi
AM_CONDITIONAL([PNL_HAVE_BLASLAPACK], [test $use_internal_blas -eq  0])

#######################
## test for underscores
## note that trailing underscores is already checked
## elsewhere
#######################

FC_TRAILING_UNDERSCORE=no

if test "x$F77" != "x"; then
   AC_MSG_CHECKING([for leading and/or trailing underscores])
cat << EOF > pipof.f
        subroutine pipof
        end
EOF
   eval "FC1=\"$F77\""
   eval "$FC1 -c pipof.f > /dev/null 2>&1"

   output=`$NM pipof.o|grep _pipof 2>&1`
   if test ! -z "$output"; then
     FC_TRAILING_UNDERSCORE=no
   fi

   output=`$NM pipof.o|grep pipof_ 2>&1`
   if test ! -z "$output"; then
     FC_TRAILING_UNDERSCORE=yes
   fi

   $RM pipof.f pipof.o
fi

if test "$FC_TRAILING_UNDERSCORE" = yes; then
   AC_DEFINE(WTU,[],[trailing underscore in fortran symbols])
fi
AC_MSG_RESULT([$FC_LEADING_UNDERSCORE $FC_TRAILING_UNDERSCORE])

#####################
## test for sharpsign
#####################

AC_MSG_CHECKING([for use of sharpsign in CPP])
AC_TRY_COMPILE([#define C2F(name) name##_],
  [C2F(toto)()],
  AC_MSG_RESULT(yes)
  AC_DEFINE(USE_SHARP_SIGN,[],[sharp sign can be used to in macros]),
  AC_MSG_RESULT(no))


#############################################################
## test for dpstrf (Block Cholesky for non definite matrix)
#############################################################

ac_have_dpstrf=no
if test $use_internal_blas -eq 0; then
    libs_save="$LIBS";
    LIBS="$LIBS $BLAS_LIBS $LAPACK_LIBS"
    AC_F77_FUNC(dpstrf)
    AC_MSG_CHECKING([for $dpstrf]) 
    AC_TRY_LINK_FUNC($dpstrf, [ac_have_dpstrf=yes])
    LIBS="$libs_save"
    AC_MSG_RESULT($ac_have_dpstrf)
fi
if test "$ac_have_dpstrf" = "yes" || test $use_internal_blas -eq 1 || 
    test "x$with_external_blas" = "xyes"; then
   AC_DEFINE(PNL_HAVE_DPSTRF, 1, [Define if dpstrf.])
fi


#################
# check for MPI
#################

MPI_LIBS=""
MPI_CFLAGS=""
MPICC=""
MPICXX=""

if test "x$enable_mpi" = xyes; then
    builtin(include,config/ax_mpi.m4)
    AX_MPI([have_mpi=1], [have_mpi=0])

    # We will use CC and not MPICC to compile and link, so we need
    # to tweak MPI_LIBS
    if test $have_mpi -eq 1; then
        MPI_CFLAGS=`$MPICC --showme:compile`
        if test "x " = x"$MPILIBS"; then
            MPI_LIBS=`$MPICC --showme:link`
        else
            MPI_LIBS="$MPILIBS"
        fi
    fi
else
    have_mpi=0
fi
AM_CONDITIONAL([PNL_HAVE_MPI], [test $have_mpi -eq 1])


####################
### checking host
####################
case "$host" in
     *-*-cygwin*)

     AC_DEFINE_UNQUOTED(
			[_WIN32], [1],
			[Running on a Win32 system])

	 AC_DEFINE_UNQUOTED(
			[_CYGWIN], [1],
			[Runnning on Cygwin])

	 AC_DEFINE_UNQUOTED(
			[_WIN32], [1],
			[Running on Cygwin and using Windows Acrobat])
     ;;

    i*86-*-mingw32*)

    AC_DEFINE_UNQUOTED(
			[_WIN32], [1],
			[Running on a Win32 system])
    ;;

    i*86-*-linux-gnu | *-pc-linux-gnu)
	;;

    # Darwin/MacOSX
    *-*-darwin*)
    ;;
esac

AC_SUBST(CC)
AC_SUBST(CXX)
AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(MKDIR_P)
AC_SUBST(BLAS_LIBS)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(BLAS_STATIC_LIBS)
AC_SUBST(LAPACK_STATIC_LIBS)
AC_SUBST(STATIC_FLAG)
AC_SUBST(LIBPNLDIR)


AC_SUBST(MPI_DIR)
AC_SUBST(LIBPNLMPI_LA)
AC_SUBST(MPICC)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_CFLAGS)


AC_OUTPUT(Makefile \
Makeuser.incl \
examples/Makefile \
src/Makefile \
src/fft/Makefile \
src/include/pnl/Makefile \
src/interpol/Makefile \
src/libamos/Makefile \
src/libblas/Makefile \
src/libcephes/Makefile \
src/libf2c/Makefile \
src/libinteg/Makefile \
src/liblapack/Makefile \
src/libminpack/Makefile \
src/librand/Makefile \
src/linalg/Makefile \
src/math/Makefile \
src/mpi/Makefile \
src/objects/Makefile \
src/optim/Makefile \
src/roots/Makefile \
src/sort/Makefile \
src/specfun/Makefile \
)

echo
echo
echo "*******************"
echo "*     Summary     *"
echo "*******************"
echo
echo "C compiler.............:" "$CC"
echo "C options..............:" "$AM_CFLAGS $CFLAGS"
echo "$user_set_cflags"
echo "linker flags...........:" "$LDFLAGS"
echo "installation prefix....:" "$prefix"
if test $have_mpi -eq 1; then
  echo "MPI detected...........:"
  echo "   MPICC...............:" "$MPICC"
  echo "   MPI_CFLAGS..........:" "$MPI_CFLAGS"
  echo "   MPI_LIBS............:" "$MPI_LIBS"
fi
if test $use_internal_blas -eq 1; then
  echo "using internal Blas & Lapack"
else
  echo "Blas lib...............:" "$BLAS_LIBS" "$BLAS_STATIC_LIBS"
  echo "Lapack lib.............:" "$LAPACK_LIBS" "$LAPACK_STATIC_LIBS"
fi
echo
echo

