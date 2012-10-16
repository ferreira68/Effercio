#include "util.h"

#include <stdarg.h>
#include <errno.h>
#include <stdio.h>

static void reason_exit_vargs(va_list args, const char *format)
{
  vfprintf(stderr,format,args);
  va_end(args);
  if(errno != 0)
  {
      fprintf(stderr,"Reason (%d): %s\n",errno,strerror(errno));
      exit(errno);
  }
  exit(1);
}

static void reason_exit(const char *format, ...)
{
  va_list args;
  va_start(args,format);
  reason_exit_vargs(args,format);
}

void full_assert(int condition,const char *format,...)
{
    va_list args;
    if(condition)
        return;
    va_start(args,format);
    reason_exit_vargs(args,format);
} 

