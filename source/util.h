#ifndef UTIL_H
#define UTIL_H

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
 * This code will take a list of jobs (names) from an input file and run *
 * them in embrassingly parallel mode on a set of processors.            *
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

#include "structs.h"

void full_assert(int condition,const char *format,...);

/**
 * Checks for the location of executables called by the run source files.
 * 
 * @param params JobParameters struct that contains the paths to be checked.
 * @return int 0 if all of the exectuables are accessible, -1 otherwise. 
 */
int check_extern_apps(const JobParameters *params);

/**
 * Checks for the ability to write to directries where data and results are
 * stored.
 *
 * @param params JobParameters struct that contains the paths to be checked.
 * @return int 0 if all of the exectuables are accessible, -1 otherwise. 
 */
int check_directories(const JobParameters *params);

#endif
