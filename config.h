/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Full path of autodock4 executable */
#define AUTODOCK_EXE "/opt/AutoDock/4.2.6/gnu/8.2.0/bin/autodock4"

/* Autodock version number. Used to determine which version of prepare_dpf
   should be used. */
#define AUTODOCK_VER "4.2"

/* Debug Mode: Will run verbose and optimization. */
/* #undef DEBUG */

/* Define to 1 if you have the `atexit' function. */
#define HAVE_ATEXIT 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the `floor' function. */
#define HAVE_FLOOR 1

/* Define to 1 if you have the `fork' function. */
#define HAVE_FORK 1

/* Define to 1 if you have the `ftruncate' function. */
#define HAVE_FTRUNCATE 1

/* Define to 1 if you have the `getcwd' function. */
#define HAVE_GETCWD 1

/* Define to 1 if you have the `gethostname' function. */
#define HAVE_GETHOSTNAME 1

/* Define to 1 if you have the `getpagesize' function. */
#define HAVE_GETPAGESIZE 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `memset' function. */
#define HAVE_MEMSET 1

/* Define to 1 if you have the `mkdir' function. */
#define HAVE_MKDIR 1

/* Define to 1 if you have a working `mmap' system call. */
#define HAVE_MMAP 1

/* Define to 1 if you have the `munmap' function. */
#define HAVE_MUNMAP 1

/* Define to 1 if the system has the type `ptrdiff_t'. */
#define HAVE_PTRDIFF_T 1

/* Define to 1 if your system has a GNU libc compatible `realloc' function,
   and to 0 otherwise. */
#define HAVE_REALLOC 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strchr' function. */
#define HAVE_STRCHR 1

/* Define to 1 if you have the `strdup' function. */
#define HAVE_STRDUP 1

/* Define to 1 if you have the `strerror' function. */
#define HAVE_STRERROR 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strrchr' function. */
#define HAVE_STRRCHR 1

/* Define to 1 if you have the `strspn' function. */
#define HAVE_STRSPN 1

/* Define to 1 if you have the `strstr' function. */
#define HAVE_STRSTR 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vfork' function. */
#define HAVE_VFORK 1

/* Define to 1 if you have the <vfork.h> header file. */
/* #undef HAVE_VFORK_H */

/* Define to 1 if `fork' works. */
#define HAVE_WORKING_FORK 1

/* Define to 1 if `vfork' works. */
#define HAVE_WORKING_VFORK 1

/* MGL Tools bin directory */
#define MGL_BIN_DIR "/opt/MGLTools/1.5.7/bin"

/* MGL Tools utility directory */
#define MGL_TOOLS "/opt/MGLTools/1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/"

/* Mopac Executable name */
#define MOPAC_EXE "mopac"

/* */
#define MOPAC_HOME "/opt/mopac/2016/"

/* Number of 1k blocks used in MPI transfer. */
#define MSG_BUFSIZ 10

/* Name of package */
#define PACKAGE "effercio"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "effercio@stjude.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "Effercio"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Effercio 1.1.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "effercio"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.1.1"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "1.1.1"

/* */
#define envARCHOSV "MGL_ARCHOSV=x86_64Linux2"

/* */
#define envLD_LIBRARY "LD_LIBRARY_PATH=/opt/MGLTools/1.5.7/lib"

/* */
#define envMGL_EXTRAINCLUDE "MGL_EXTRAINCLUDE=/opt/MGLTools/1.5.7/include"

/* */
#define envMGL_EXTRALIBS "MGL_EXTRALIBS=/opt/MGLTools/1.5.7/lib"

/* */
#define envMGL_ROOT "MGL_ROOT=/opt/MGLTools/1.5.7"

/* */
#define envMOPAC_LIC "/opt/mopac/2016/"

/* */
#define envPATH "PATH=/opt/openmpi/3.1.6/gnu/8.2.0/bin:/opt/gcc/8.2.0/bin:/usr/lib64/qt-3.3/bin:/usr/local/cuda/bin:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin:/home/aferreira/Software/go/bin:/opt/MGLTools/1.5.7"

/* */
#define envPYTHONHOME "PYTHONHOME=/opt/MGLTools/1.5.7"

/* */
#define envPYTHONPATH "PYTHONPATH=/opt/MGLTools/1.5.7/MGLToolsPckgs"

/* */
#define envTCL_LIBRARY "TCL_LIBRARY=/opt/MGLTools/1.5.7/tcl8.4"

/* */
#define envTK_LIBRARY "TK_LIBRARY=/opt/MGLTools/1.5.7/tk8.4"

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to rpl_realloc if the replacement function should be used. */
/* #undef realloc */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define as `fork' if `vfork' does not work. */
/* #undef vfork */
