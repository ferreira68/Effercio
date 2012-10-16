// $Id: extract_index.c 9 2012-03-01 20:01:25Z dcoss $
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
 *************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "defines.h"
#include "memwatch.h"


int ExtractIndex(char *input, char key)
{
	char *tempstr = NULL;

    size_t numdigits;
    int retval;

    input = strchr(input,'_');
    if (input == NULL) {
	retval = 1;
    }
    else {
	input = strchr(input,key) + 1;
	numdigits = strspn(input,THE_DIGITS);
	tempstr = (char *) malloc((numdigits + 1) * sizeof(char));
	strncpy(tempstr,input,numdigits+1);
	tempstr[numdigits] = (char) NULL;

	retval = atoi(tempstr);
    }

    free(tempstr);
    tempstr = NULL;

    return retval;
}

