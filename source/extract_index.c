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
	tempstr[numdigits] = 0;

	retval = atoi(tempstr);
    }

    free(tempstr);
    tempstr = NULL;

    return retval;
}

