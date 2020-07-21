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

#include <errno.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "defines.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "structs.h"
#include "mpi.h"
#include "memwatch.h"


void PrintBanner()
{
    printf("%s\n\n",PLUS_HDR);

    printf("   oooooooooooo  .o88o.  .o88o.                               o8o            \n");
    printf("   `888'     `8  888 `\"  888 `\"                               `\"'            \n");
    printf("    888         o888oo  o888oo   .ooooo.  oooo d8b  .ooooo.  oooo   .ooooo.  \n");
    printf("    888oooo8     888     888    d88' `88b `888\"\"8P d88' `\"Y8 `888  d88' `88b \n");
    printf("    888    \"     888     888    888ooo888  888     888        888  888   888 \n");
    printf("    888       o  888     888    888    .o  888     888   .o8  888  888   888 \n");
    printf("   o888ooooood8 o888o   o888o   `Y8bod8P' d888b    `Y8bod8P' o888o `Y8bod8P' \n");

/*
    printf("             _         _            _            _    __  __ ____ ___ \n");
    printf("            / \\  _   _| |_ ___   __| | ___   ___| | _|  \\/  |  _ \\_ _|\n");
    printf("           / _ \\| | | | __/ _ \\ / _` |/ _ \\ / __| |/ / |\\/| | |_) | | \n");
    printf("          / ___ \\ |_| | || (_) | (_| | (_) | (__|   <| |  | |  __/| | \n");
    printf("         /_/   \\_\\__,_|\\__\\___/ \\__,_|\\___/ \\___|_|\\_\\_|  |_|_|  |___|\n");
*/

    printf("\n%s\n",PLUS_HDR);

    printf("\nAuthors: Antonio M. Ferreira, PhD\n");
    printf("         High Performance Computing Facility\n");
    printf("         Research Informatics/Information Sciences\n");
    printf("         and\n");
    printf("         Department of Structural Biology\n");
    printf("         St. Jude Children's Research Hospital\n");
    printf("e-mail : Antonio.Ferreira@stjude.org\n\n");
    printf("         David Coss, PhD\n");
    printf("         High Performance Computing Facility\n");
    printf("         Research Informatics/Information Sciences\n");
    printf("         St. Jude Children's Research Hospital\n");
    printf("e-mail : David.Coss@stjude.org\n\n\n");
    fflush(stdout);
}


void printf_master(const char *format, ...)
{
    va_list args;
    int rank;

    fflush(NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == MASTER) {
	va_start(args, format);
	vprintf(format,args);
	va_end(args);
    }
    fflush(stdout);
}


void printf_verbose(int verbose, const char *format, ...)
{
    va_list args;

    fflush(NULL);

    if (verbose) {
	va_start(args, format);
	vprintf(format,args);
	va_end(args);
    }
    fflush(stdout);
}


int CopyFile(char source_name[FILENAME_MAX],char destination_name[FILENAME_MAX])
{
    int retval = 0;
    size_t num_read = 0;
    size_t num_write = 0;
    unsigned char *filebuff;

    FILE *source, *destination;

    // Allocate a copy buffer
    filebuff = calloc(FILE_BUFF_SIZE, sizeof(unsigned char));
    if (filebuff == NULL) {
        printf("ERROR: Failed to allocate %lu bytes for file copy buffer!\n",
	       (unsigned long) FILE_BUFF_SIZE);
        printf("ERROR: %s\n",strerror(errno));
        exit(errno);
    }
    // Open the source file
    source = fopen(source_name,"rb");
    if (source == NULL) {
        printf("ERROR: Failed to open source file (%s) for file copy!\n",
	       source_name);
        printf("ERROR: %s\n",strerror(errno));
        exit(errno);
    }
    // Open the destination file
    destination = fopen(destination_name,"w+b");
    if (destination == NULL) {
        printf("ERROR: Failed to open destination file (%s) for file copy!\n",
	       destination_name);
        printf("ERROR: %s\n",strerror(errno));
        exit(errno);
    }
    // Copy the file
    while (!feof(source)) {
        num_read  = fread(filebuff,sizeof(unsigned char),FILE_BUFF_SIZE,source);
	if (num_read == 0) {
            // Need to check feof and ferror here for proper error checking
            printf("ERROR: Read of source file (%s) failed during file copy!\n",
		   destination_name);
	    exit(EXIT_FAILURE);
	}
	num_write = fwrite(filebuff,sizeof(unsigned char),num_read,destination);
	if (num_read != num_write) {
            // Need to check feof and ferror here for proper error checking
            printf("ERROR: Write to destination file (%s) failed during file copy!\n",
		   destination_name);
	    exit(EXIT_FAILURE);
	}
    }
    fflush(destination);
    fsync(fileno(destination));
    fclose(destination);
    fclose(source);
    // Free the transfer buffer and set it to NULL to prevent memory leaks
    free(filebuff);
    filebuff = NULL;
    // We're done
    return(retval);
}


int MoveFile(char from_name[FILENAME_MAX],char to_name[FILENAME_MAX])
{
    int retval = 0;

    retval = CopyFile(from_name,to_name);
    retval = remove(from_name);

    return(retval);
}


#include <math.h>
void GetValueByColumn(const char *string, const char *key, int data_column,
                      double *retval)
{

    int column = 0;
    char *reference;
    char *substring;
    char *endpointer;

    // Initialize the pointers
    reference  = malloc((strlen(string)+1) * sizeof(char));

    if (reference == NULL) {
        printf("ERROR - Call to malloc() failed in GetValueByColumn!!\n");
	exit(EXIT_FAILURE);
    }
    strcpy(reference,string);
    substring = strtok(reference,key);
    while (substring != NULL) {
        column++;
        if (column == data_column) {
            *retval = strtod(substring,&endpointer);
	    break;
	}
	else {
	    substring = strtok(NULL,key);
	}
    }
    // Free up the pointers
    free(reference);

    reference  = NULL;
    substring  = NULL;
    endpointer = NULL;

    return;
}


void GetStringByColumn(const char *string, const char *key, int data_column,
                      char **retstr)
{
    #include <math.h>

    int column = 0;
    char *reference;
    char *substring;

    // Initialize the pointers
    reference  = malloc((strlen(string)+1) * sizeof(char));

    if (reference == NULL) {
        printf("ERROR - Call to malloc() failed in GetValueByColumn!!\n");
	exit(EXIT_FAILURE);
    }
    strcpy(reference,string);
    substring = strtok(reference,key);
    while (substring != NULL) {
        column++;
        if (column == data_column) {
            (*retstr) = realloc((*retstr),(strlen(substring)+1)*sizeof(char));
            memset(*retstr,0, (strlen(substring)+1)*sizeof(char));
            strcpy((*retstr),substring);
	    break;
	}
	else {
	    substring = strtok(NULL,key);
	}
    }
    // Free up the pointers
    free(reference);
    reference  = NULL;
    substring  = NULL;

    return;
}


char *SecondsToHMS(double num_secs)
{
    char *result = calloc(20,sizeof(char));

    int hours = 0;
    int minutes = 0;
    double seconds;

    seconds = num_secs;

    if (seconds > 3600) {
        hours = (int) seconds/3600;
        seconds = seconds - (double) hours*3600;
    }
    if (seconds > 60) {
        minutes = (int) seconds/60;
        seconds = seconds - (double) minutes*60;
    }

    if (hours != 0) {
        sprintf(result,"%d:%02d:%02.0f hours",hours,minutes,seconds);
    }
    else if (minutes != 0) {
        sprintf(result,"%2d:%04.1f minutes",minutes,seconds);
    }
    else {
        sprintf(result,"%.2f seconds",seconds);
    }

    return result;
}


int VerifyDir(const char dirname[FILENAME_MAX], int verbose, const char nodestr[HOST_NAME_MAX])
{
    // This function checks to see if 'dirname' exists.
    // If it does not, it is created.

    int retval = 0;
    struct stat status_buffer;

    // Look for a scratch directory and create one if it doesn't exist
    if (stat(dirname,&status_buffer) == 0) {
	if (verbose) printf("\n%s: Found directory %s\n",nodestr,dirname);
    }
    else {
        if (verbose) printf("\n%s: Directory %s not found!\n",nodestr,dirname);
	mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
	// Need error checking in here!!!
	if (stat(dirname,&status_buffer) == 0) {
	    if (verbose) printf("%s: Successfully created %s\n",nodestr,dirname);
	}
	else {
	    if (verbose) printf("%s: Failed to create %s\n",nodestr,dirname);
	    retval = -1;
	}
    }

    return retval;
}
