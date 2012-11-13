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
 *                                                                       *
 * Effercio is free software: you can redistribute it and/or modify      *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * Effercio is distributed in the hope that it will be useful,           *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with Effercio. If not, see <http://www.gnu.org/licenses/>.      *
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
