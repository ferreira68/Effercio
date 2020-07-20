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

#include "util.h"

#include <stdarg.h>
#include <errno.h>
#include <stdio.h>
#define _GNU_SOURCE             /* See feature_test_macros(7) */
#include <unistd.h>
#include <string.h>


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

/*
static void reason_exit(const char *format, ...)
{
    va_list args;
    va_start(args,format);
    reason_exit_vargs(args,format);
}
*/

void full_assert(int condition,const char *format,...)
{
    va_list args;
    if(condition) return;

    va_start(args,format);
    reason_exit_vargs(args,format);
} 

int check_extern_apps(const JobParameters *params)
{
    int retval;
    if(params == NULL) return -1;

    retval = access(params->mopac_path,X_OK);
    if(retval != 0)
    {
        printf ("%s ERROR: MOPAC executable not found at %s",params->node_tag,params->mopac_path);
        return retval;
    }
  
    retval = access(params->autodock_exe,X_OK);
    if(retval != 0)
    {
        printf ("%s ERROR: Autodock executable not found at %s",params->node_tag,params->autodock_exe);
        return retval;
    }
  
    retval = access(params->mgl_bin_dir,R_OK);
    if(retval != 0)
    {
        printf ("%s ERROR: MGL Tools not found in %s",params->node_tag,params->mgl_bin_dir);
        return retval;
    }
  
    return 0;
}

int check_directories(const JobParameters *params)
{
    int retval;

    // check for no parameters
    if(params == NULL) return -1;
    //
    // Check the receptor .map files directory
    retval = euidaccess(params->receptor_dir ,W_OK);
    if(retval != 0)
    {
        printf("%s ERROR: Directory check failed for %s",params->node_tag,params->receptor_dir);
        return retval;
    }
    //
    // Check the ligands directory
    retval = euidaccess(params->ligand_dir   ,W_OK);
    if(retval != 0)
    {
        printf("%s ERROR: Directory check failed for %s\n",params->node_tag,params->ligand_dir);
        return retval;
    }
    //
    // Check the docking results directory
    retval = euidaccess(params->results_dir  ,W_OK);
    if(retval != 0)
    {
        printf("%s ERROR: Directory check failed for %s\n",params->node_tag,params->results_dir);
        return retval;
    }
    //
    // Check the cluster representatives directory
    retval = euidaccess(params->clusters_dir ,W_OK);
    if(retval != 0)
    {
        printf("%s ERROR: Directory check failed for %s\n",params->node_tag,params->clusters_dir);
        return retval;
    }
    //
    // Check the directory for SEQM optimized structures
    retval = euidaccess(params->optimized_dir,W_OK);
    if(retval != 0)
    {
        printf("%s ERROR: Directory check failed for %s\n",params->node_tag,params->optimized_dir);
        return retval;
    }
    //
    // Check the directory for final results
    retval = euidaccess(params->analysis_dir ,W_OK);
    if(retval != 0)
    {
        printf("%s ERROR: Directory check failed for %s\n",params->node_tag,params->analysis_dir);
        return retval;
    }
    //
    // Check the scratch directory
    retval = euidaccess(params->scratch_dir ,W_OK);
    if(retval != 0)
    {
        printf("%s ERROR: Directory check failed for %s\n",params->node_tag,params->scratch_dir);
        return retval;
    }

    return 0;

}
