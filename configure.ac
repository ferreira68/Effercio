#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([Effercio], [1.1.1], [effercio@stjude.org])
AC_CONFIG_SRCDIR([source/])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([1.9])

# Checks for programs.
AC_PROG_CC([mpicc])
AM_PROG_CC_C_O
AC_PROG_INSTALL

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lefence':
#AC_CHECK_LIB([efence], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [sin])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h inttypes.h limits.h stddef.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_TYPE_INT16_T
#AC_TYPE_INT32_T
#AC_TYPE_INT64_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
#AC_TYPE_UINT16_T
#AC_TYPE_UINT32_T
#AC_TYPE_UINT64_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_FORK
AC_CHECK_FUNCS([atexit floor ftruncate getcwd gethostname memset mkdir munmap strchr strdup strerror strrchr strspn strstr])

#static build?
AC_MSG_CHECKING([static build])
AC_ARG_ENABLE([static],
	[AS_HELP_STRING([--enable-static],[Enable building of effercio statically. Default: enabled. To link against shared libraries, use --disable-static])],
	[enable_static="$enableval"],
	[enable_static=yes])
AC_MSG_RESULT($enable_static)
if test x"$enable_static" = x"yes";then
   CFLAGS="$CFLAGS -static "
fi   

#debug
AC_MSG_CHECKING([use debug mode])
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],[Enable debug data logging. Default: no])],
	[enable_debug="$enableval"],
	[enable_debug=no])
AC_MSG_RESULT([$enable_debug])
if test x"$enable_debug" = x"yes"; then
   AC_DEFINE([DEBUG],[],[Debug Mode: Will run verbose and optimization.])
   CFLAGS="$CFLAGS -g -O0 "
fi

# buffer size
AC_MSG_CHECKING([transfer block size])
AC_ARG_WITH([buf-size],
	[AS_HELP_STRING([--with-buff-size=<size],[Set the number of blocks used to transfer data to master node. One block is a kilobyte. Default: 10 blocks.])],
	[with_buff_size="$withval"],
	[with_buff_size=10])
AC_DEFINE_UNQUOTED([MSG_BUFSIZ],[$with_buff_size],[Number of 1k blocks used in MPI transfer.])
AC_MSG_RESULT([$with_buff_size kilobytes])

# Scratch space
AC_MSG_CHECKING([scratch space])
AC_ARG_WITH([scratch],
	[AS_HELP_STRING([--with-scratch=<dir>],[Specify  directory to be used as scratch space. Default: "/scratch_space"])],
	[with_scratch="$withval"],
	[with_scratch=/scratch_space])
[
case "$with_scratch" in
     no) echo Scratch space is required. ; exit -1;;
     *) echo "$with_scratch";;
esac
]

#Autodock
AC_MSG_CHECKING([autodock version])
AC_ARG_WITH([autodock-version],
	[AS_HELP_STRING([--with-autodock-version@<:@version@:>@],[Specify the autodock version. Default: 4.2])],
	[with_autodock_version="$withval"],
	[with_autodock_version=4.2])
[
case "$with_autodock_version" in
     4.2|4.1|4.0|4) echo "$with_autodock_version";;
     *) echo Unsupported version of Autodock. ;exit -1;;
esac
]
AC_DEFINE_UNQUOTED([AUTODOCK_VER],["$with_autodock_version"],[Autodock version number. Used to determine which version of prepare_dpf should be used.])
AC_MSG_CHECKING([autodock path])
AC_ARG_WITH([autodock-path],
	[AS_HELP_STRING([--with-autodock-path=dir],[Specify the path of autodock. Required])],#should be require?? Assume in path perhaps?
	[with_autodock_path="$withval"],
	AC_MSG_ERROR([Autodock path required.]))
AC_MSG_RESULT([$with_autodock_path])
AC_DEFINE_UNQUOTED([AUTODOCK_EXE],["$with_autodock_path"],[Full path of autodock4 executable])

#MGL tools
AC_MSG_CHECKING([MGL Tools path])
AC_ARG_WITH([mgltools],
	[AS_HELP_STRING([--with-mgltools],[Specify the required MGL tools path])],
	[with_mgltools="$withval"],
	AC_MSG_ERROR([MGLTools is required. Provide the path with --with-mgltools]))
AC_MSG_RESULT([$with_mgltools])

AH_TEMPLATE([MGL_TOOLS])
AC_DEFINE_UNQUOTED([MGL_TOOLS],["$with_mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/"],[MGL Tools utility directory])
AH_TEMPLATE([MGL_BIN_DIR])
AC_DEFINE_UNQUOTED([MGL_BIN_DIR],["$with_mgltools/bin"],[MGL Tools bin directory])
AH_TEMPLATE([envMGL_ROOT])
AC_DEFINE_UNQUOTED([envMGL_ROOT],["MGL_ROOT=$with_mgltools"],[])
AH_TEMPLATE([envARCHOSV])
AC_DEFINE_UNQUOTED([envARCHOSV],["MGL_ARCHOSV=x86_64Linux2"],[])
AH_TEMPLATE([envMGL_EXTRALIBS])
AC_DEFINE_UNQUOTED([envMGL_EXTRALIBS],["MGL_EXTRALIBS=$with_mgltools/lib"],[])
AH_TEMPLATE([envMGL_EXTRAINCLUDE])
AC_DEFINE_UNQUOTED([envMGL_EXTRAINCLUDE],["MGL_EXTRAINCLUDE=$with_mgltools/include"],[])
AH_TEMPLATE([envTCL_LIBRARY])
AC_DEFINE_UNQUOTED([envTCL_LIBRARY],["TCL_LIBRARY=$with_mgltools/tcl8.4"],[])
AH_TEMPLATE([envTK_LIBRARY])
AC_DEFINE_UNQUOTED([envTK_LIBRARY],["TK_LIBRARY=$with_mgltools/tk8.4"],[])
AH_TEMPLATE([envLD_LIBRARY])
AC_DEFINE_UNQUOTED([envLD_LIBRARY],["LD_LIBRARY_PATH=$with_mgltools/lib"],[])
AH_TEMPLATE([envPYTHONHOME])
AC_DEFINE_UNQUOTED([envPYTHONHOME],["PYTHONHOME=$with_mgltools"],[])
AH_TEMPLATE([envPYTHONPATH])
AC_DEFINE_UNQUOTED([envPYTHONPATH],["PYTHONPATH=$with_mgltools/MGLToolsPckgs"],[])
AH_TEMPLATE([envPATH])
AC_DEFINE_UNQUOTED([envPATH],["PATH=$PATH:$with_mgltools"],[])

#mopac
AC_MSG_CHECKING([Mopac])
AC_ARG_WITH([mopac],
	[AS_HELP_STRING([--with-mopac=dir],[Specify the Mopac directory. License file must be placed in the Mopac directory. Default: /opt/mopac])],
	[with_mopac="$withval"],
	[with_mopac=/opt/mopac])
AC_MSG_RESULT([$with_mopac])
AH_TEMPLATE([MOPAC_HOME])
AH_TEMPLATE([MOPAC_HOME])
AH_TEMPLATE([envMOPAC_LIC])
AC_DEFINE_UNQUOTED([MOPAC_EXE],["mopac"],[Mopac Executable name])
AC_DEFINE_UNQUOTED([MOPAC_HOME],["$with_mopac/"],[])
AC_DEFINE_UNQUOTED([envMOPAC_LIC],["$with_mopac/"],[])

AC_CONFIG_FILES([Makefile
                 source/Makefile])
AC_OUTPUT
