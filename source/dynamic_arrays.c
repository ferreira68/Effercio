// $Id: dynamic_arrays.c 9 2012-03-01 20:01:25Z dcoss $
/*************************************************************************
 * Authors: Antonio M. Ferreira, PhD [1,2]                               *
 *          David Coss, PhD [1]                                          *
 *                                                                       *
 *          (1) High Performance Computing Faclity                       *
 *              Research Informatics, Information Sciences               *
 *                                                                       *
 *          (2) Structural Biology                                       *
 *                                                                       *
 *          St. Jude Children's Research Hospital                        *
 *                                                                       *
 * This file contains routines to manipulate dynamic arrays.  It is very *
 * rudimentary at present and lacks proper error checking, but it gets   *
 * the job done for the moment.                                          *
 *************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memwatch.h"


void add_strtoarray(char ***str_array, int *num_elements,
                    const char *my_string)
{
    // NOTE - This function is used to create NULL TERMINATED arrays of strings.  That
    //        is, the last element of the array must be the NULL pointer.  This is for
    //        creating argument lists for execve and others.  So KEEP TRACK of that
    //        last element!

    (*num_elements)++;
    *str_array = realloc(*str_array,((*num_elements + 1)*sizeof(char **)));
    if (*str_array == NULL) {
        printf("ERROR - Failed to allocate memory for string array\n");
        printf("ERROR - %s\n",strerror(errno));
        exit(errno);
    }

    (*str_array)[(*num_elements) - 1] = malloc((strlen(my_string)+1)*sizeof(char *));
    strcpy((*str_array)[(*num_elements) - 1],my_string);
    (*str_array)[(*num_elements)] = (char *)NULL;

}


void FreeStringArray(char ***str_array, int *num_elements)
{
    int i = *num_elements;
    if(str_array == NULL || *str_array == NULL)
    	return;

    (*str_array)[*num_elements] = strdup("junk");

    for (i=0; i <= *num_elements; i++) {
	free((*str_array)[i]);
	(*str_array)[i] = NULL;
    }
    free(*str_array);
    *str_array = NULL;
    *num_elements = 0;

    return;
}


void add_inttoarray(int **int_array, int *num_elements,
		    const int my_int)
{
    (*num_elements)++;
    *int_array = (int *) realloc(*int_array,
				  (*num_elements)*sizeof(int));
    (*int_array)[(*num_elements) - 1] = my_int;
    
}
