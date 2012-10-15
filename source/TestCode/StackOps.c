// $Id: StackOps.c 6 2011-11-17 19:55:32Z dcoss $
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define OUTPUT_WIDTH 80;
//#include "defines.h"
#include "../structs.h"
#include "mpi.h"

int main(void)
{
    char *infile = "DebugList";
    jobnode *jobstack = NULL;
    int i,num;
    char nextjob[FILENAME_MAX];

    PushStack(&jobstack,"Test 1");
    PushStack(&jobstack,"Test 2");
    PushStack(&jobstack,"Test 3");
    PushStack(&jobstack,"Test 4");
    PushStack(&jobstack,"Test 5");
    PushStack(&jobstack,"Test 6");
    PushStack(&jobstack,"Test 7");
    PushStack(&jobstack,"Test 8");

    PrintStack(jobstack,"JOBSTACK after creation");
    num = CountStack(jobstack);

    i = 0;
    while(jobstack != NULL) {
        i++;
        printf("i = %d\n",i);
        printf("stack pointer = %p\n",jobstack);
        PopStack(&jobstack,&nextjob);
        printf("top item      = %s\n\n",nextjob);
    }
  
    return EXIT_SUCCESS;

}
