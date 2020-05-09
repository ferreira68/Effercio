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

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <unistd.h>
#include "defines.h"
#include "structs.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "memwatch.h"


int RunDock(job_t *job,JobParameters *params, struct STICelement *results)
{
    static int first_job = TRUE;
    //const int FILE_BUFF_SIZE = 2 * 1024 * 1024;
    char filestub[FILENAME_MAX];
    char initdir[FILENAME_MAX];
    char workdir[FILENAME_MAX];
    char receptPDB[FILENAME_MAX];
    char receptStruct[FILENAME_MAX];
    char ligandStruct[FILENAME_MAX];
    char dpf_file[FILENAME_MAX],dlg_file[FILENAME_MAX];
    char tmp_name[FILENAME_MAX];
    char tmp_name2[FILENAME_MAX];
    char cmdstr[ARG_MAX];
    char line[ARG_MAX];
    char cmplx_filename[FILENAME_MAX],cmplxstub[FILENAME_MAX];
    char qmligand_filename[FILENAME_MAX];
    int retval = 0;
    int done = FALSE;
    int i,numclusters;
    long int resultsPos,startmodel[STRUCTS_MAX];
    FILE *dockfile,*complexfile,*receptorfile,*srcfile,*destfile;
    FILE *qmligandfile;
    char prepare_script_name[FILENAME_MAX],dockOpts[ARG_MAX];
    time_t timestamp,autodock_time;
    pid_t pid;
    char **MGL_args = NULL;
    int num_MGL_args = 0;
    char **MGL_env = NULL;
    int num_MGL_env = 0;

    // Parameters pass along by Effercio.c
    char *user = params->username;
    char *nodestr = params->node_tag;
    char *receptDir = params->receptor_dir;
    char *receptName = params->receptor_name;
    char *ligandDir = params->ligand_dir;
    char *resultsDir = params->results_dir;
    char *clustersDir = params->clusters_dir;
    char *scratchDir = params->scratch_dir;
    char **dockingParams = params->dock_params;
    int numParams = params->num_dock_params;
    int verbose = params->verbose;

    // Initialize arrays
    memset(filestub,(char) NULL,FILENAME_MAX);
    memset(initdir,(char) NULL,FILENAME_MAX);
    memset(workdir,(char) NULL,FILENAME_MAX);
    memset(receptPDB, (char) NULL,FILENAME_MAX);
    memset(receptStruct,(char) NULL,FILENAME_MAX);
    memset(ligandStruct,(char) NULL,FILENAME_MAX);
    memset(dpf_file,(char) NULL,FILENAME_MAX);
    memset(dlg_file,(char) NULL,FILENAME_MAX);
    memset(tmp_name,(char) NULL,FILENAME_MAX);
    memset(tmp_name2,(char) NULL,FILENAME_MAX);
    memset(cmdstr,(char) NULL,ARG_MAX);
    memset(line,(char) NULL,ARG_MAX);
    memset(cmplx_filename,(char) NULL,FILENAME_MAX);
    memset(cmplxstub,(char) NULL,FILENAME_MAX);
    memset(qmligand_filename,(char) NULL,FILENAME_MAX);
    memset(prepare_script_name,(char) NULL,FILENAME_MAX);
    memset(dockOpts,(char) NULL,ARG_MAX);

/*
    // Use the following to allow a debugger to attach to the slave process
    #ifdef DEBUG
    int debugWait = TRUE;
    while (debugWait) {
        printf("%s: DEBUG - Process %d is waiting!\n",nodestr,getpid());
        fflush(NULL);
        system("sleep 60");
    }
    #endif
*/

    // Initialize the results data structure
    #ifdef DEBUG
    printf("%s: DEBUG - Initializing variables in RunDock\n",nodestr);
    #endif

    results->S = ExtractIndex(job->name,'S');
    results->T = ExtractIndex(job->name,'T');
    results->I = ExtractIndex(job->name,'I');
    results->C = ExtractIndex(job->name,'C');

    #ifdef DEBUG
    printf("%s: DEBUG - S index = %4d\n",nodestr,results->S);
    printf("%s: DEBUG - T index = %4d\n",nodestr,results->T);
    printf("%s: DEBUG - I index = %4d\n",nodestr,results->I);
    printf("%s: DEBUG - C index = %4d\n",nodestr,results->C);
    #endif

    // Set up some strings for file names
    // TODO: need to do length checking to make sure we don't go longer than
    //       the name arrays
    getcwd(initdir,FILENAME_MAX);

    strcpy(workdir,scratchDir);
    strcat(workdir,"/");
    strcat(workdir,user);

    strcpy(filestub,job->name);

    strcpy(ligandStruct,filestub);
    strcat(ligandStruct,".pdbqt");

    strcpy(receptPDB,receptName);
    strcat(receptPDB,".pdb");

    strcpy(receptStruct,receptName);
    strcat(receptStruct,".pdbqt");

    strcpy(cmplxstub,receptName);
    strcat(cmplxstub,"_");
    strcat(cmplxstub,filestub);

    strcpy(dpf_file,filestub);
    strcat(dpf_file,"_");
    strcat(dpf_file,receptName);
    strcat(dpf_file,".dpf");

    strcpy(dlg_file,filestub);
    strcat(dlg_file,"_");
    strcat(dlg_file,receptName);
    strcat(dlg_file,".dlg");

    // Look for a scratch directory and create one if it doesn't exist
    if (VerifyDir(workdir,verbose,nodestr) == 0) {
	chdir(workdir);
    }
    else {
	printf("Cannot change to directory %s",workdir);
	exit(-1);
    }

    // For now let's just do everything with C calls
    // TODO: Let MPI handle all of these file transfers
    
    // Copy the .map files to the working directory
    for (i=1;i<OUTPUT_WIDTH;i++) printf_verbose(verbose,"-"); printf_verbose(verbose,"\n");
    printf_verbose(verbose,"%s: Ligand     = %s\n",nodestr,job->name);
    printf_verbose(verbose,"%s: Receptor   = %s\n",nodestr,receptName);
    printf_verbose(verbose,"%s: Source Dir = %s\n",nodestr,receptDir);


    // Copy the ligand file to the local scratch directory
    sprintf(tmp_name,"%s/%s",ligandDir,ligandStruct);
    CopyFile(tmp_name,ligandStruct);

    // Decide which script to use
    if (AUTODOCK_VER == "4.0") {
	strcpy(prepare_script_name,"prepare_dpf4.py");
    }
    else if (AUTODOCK_VER == "4.1") {
	strcpy(prepare_script_name,"prepare_dpf41.py");
    }
    else if (AUTODOCK_VER == "4.2") {
	strcpy(prepare_script_name,"prepare_dpf42.py");
    }

    #ifdef DEBUG
    printf("%s: DEBUG - Preparation script = %s/%s\n",nodestr,MGL_TOOLS,prepare_script_name);
    #endif

    // Set the argument list for the prepare script
    // The first argument should be the file name
    add_strtoarray(&MGL_args,&num_MGL_args,prepare_script_name);
    
    // Set the ligand file
    add_strtoarray(&MGL_args,&num_MGL_args,"-l");
    add_strtoarray(&MGL_args,&num_MGL_args,ligandStruct);

    // Set the receptor file
    add_strtoarray(&MGL_args,&num_MGL_args,"-r");
    add_strtoarray(&MGL_args,&num_MGL_args,receptStruct);
    
    // Set the parameters to be used for docking
    if (numParams > 0) {
	for(i=0;i<numParams;i++){
	    add_strtoarray(&MGL_args,&num_MGL_args,"-p");
	    add_strtoarray(&MGL_args,&num_MGL_args,dockingParams[i]);
	}
    }

    // Dump the argument list for DEBUG purposes
    #ifdef DEBUG
    printf("%s: DEBUG - MGL_args is:\n",nodestr);
    for (i=0;i<num_MGL_args;i++) {
        printf("%s: DEBUG - MGL_args[%d] = %s\n",nodestr,i,MGL_args[i]);
    }
    #endif

    // Set up the environment for the prepare script
    // NOTE: These are obtained by running the get_environment script in the source directory
    add_strtoarray(&MGL_env,&num_MGL_env,envMGL_ROOT);
    add_strtoarray(&MGL_env,&num_MGL_env,envARCHOSV);
    add_strtoarray(&MGL_env,&num_MGL_env,envMGL_EXTRALIBS);
    add_strtoarray(&MGL_env,&num_MGL_env,envMGL_EXTRAINCLUDE);
    add_strtoarray(&MGL_env,&num_MGL_env,envTCL_LIBRARY);
    add_strtoarray(&MGL_env,&num_MGL_env,envTK_LIBRARY);
    add_strtoarray(&MGL_env,&num_MGL_env,envLD_LIBRARY);
    add_strtoarray(&MGL_env,&num_MGL_env,envPYTHONHOME);
    add_strtoarray(&MGL_env,&num_MGL_env,envPYTHONPATH);
    add_strtoarray(&MGL_env,&num_MGL_env,envPATH);
    //sprintf(cmdstr,"PWD=%s",getenv("PWD"));
    sprintf(cmdstr,"PWD=%s",workdir);
    add_strtoarray(&MGL_env,&num_MGL_env,cmdstr);

    // Dump the environment list for DEBUG purposes
    #ifdef DEBUG
    printf("%s: DEBUG - MGL_env is:\n",nodestr);
    for (i=0;i<num_MGL_env;i++) {
      printf("%s: DEBUG - MGL_env[%d] = %s\n",nodestr,i,MGL_env[i]);
    }
    #endif

    // Run the script to create the .dpf file here
    sprintf(cmdstr,"%s/%s",MGL_TOOLS,prepare_script_name);
    if (verbose) {
	printf("%s: Calling %s to create docking parameter file\n",nodestr,prepare_script_name);
	fflush(stdout);
    }
    #ifdef DEBUG
    printf("%s: DEBUG - cmdstr = %s\n",nodestr,cmdstr);
    fflush(stdout);
    #endif
    pid = vfork();
    if (pid == -1) {
	// fork error code here
	printf("%s: ERROR - Failed to initialize child process via vfork()\n",nodestr);
	printf("%s: ERROR - %s\n",nodestr,strerror(errno));
	retval = -1;
	return(retval);
    }
    if (pid == 0) {
      int prepare_retval = execve(cmdstr,MGL_args,MGL_env);
      if(prepare_retval)
	{
	  printf("%s: ERROR - parameter file creation ended with the following error code: %d",prepare_retval);
	  return -1;
	}
    }
    else {
	wait(NULL);
    }
    if (verbose) {
	printf("%s: %s has finished\n",nodestr,prepare_script_name);
	fflush(stdout);
    }


    // Run Autodock here
    timestamp = time(NULL);
    destfile = fopen(dlg_file, "w+");
    if (destfile == NULL) {
	printf("%s: ERROR - Failed to open docking results file (%s)\n",nodestr,dlg_file);
	printf("%s: ERROR - %s\n",nodestr,strerror(errno));
	retval = -1;
	return(retval);
    }
    fprintf(destfile,"%s: began execution on %s",nodestr,asctime(localtime(&timestamp)));
    fsync(fileno(destfile));
    fclose(destfile);
    printf_verbose(verbose,"%s: Running AutoDock on %s\n",nodestr,dpf_file);
    sprintf(cmdstr,"%s -p %s >> %s",AUTODOCK_EXE,dpf_file,dlg_file);
    #ifdef DEBUG
    printf("%s: DEBUG - cmdstr = %s \n",nodestr,cmdstr);
    #endif
    retval = system(cmdstr);
    printf_verbose(verbose,"%s: AutoDock has finished with return value %d\n",nodestr,retval);
    
    autodock_time = time(NULL);

    #ifndef DEBUG
    // NOTE: For debugging purposes, the following two remove calls will
    //       be commented out, leaving the files in the scratch directory
    //
    // Remove the ligand structure file
    if(retval == 0)
    remove(ligandStruct);

    // Remove the docking parameter file (it's in the log file now)
    if(retval == 0)
      remove(dpf_file);
    #endif

    // Extract low-energy structures here
    dockfile = fopen(dlg_file,"r");
    if (dockfile == NULL) {
	printf("%s: ERROR - Failed to open docking results file (%s)!\n",nodestr,dlg_file);
	return (-1);
    }
    // Find out how many clusters we have
    rewind(dockfile);
    for (i=0; i < STRUCTS_MAX ; i++) startmodel[i] = 0;
    numclusters = 0;
    resultsPos = 0;
    while (fgets(line,ARG_MAX,dockfile) != NULL) {
	if (strstr(line,DOCK_RESULTS_HDR) != NULL) {
	    resultsPos = ftell(dockfile);
	}
	if (strstr(line,CLUSTER_MARKER) != NULL) {
	    startmodel[numclusters] = ftell(dockfile);
	    numclusters++;
	}
    }
    if (resultsPos == 0) {
	printf("%s: No results found in file (%s)!!\n",nodestr,dlg_file);
	return(DEFAULT_ERR_CODE);
    }
    #ifdef DEBUG
    printf("%s: DEBUG - Found %d clusters in %s\n",nodestr,numclusters,dlg_file);
    printf("%s: DEBUG - at locations ",nodestr);
    for (i=0 ; i < numclusters ; i++) printf("%ld ",startmodel[i]);
    printf("\n");
    fflush(stdout);
    #endif

    // Generate the cluster representatives
    // Open the receptor file
    receptorfile = fopen(receptStruct,"r");
    // receptorfile = fopen(receptPDB,"r");
    if (receptorfile == NULL) {
	printf("%s: ERROR - Failed to open receptor structure file (%s)!\n",
	       nodestr,receptPDB);
	return (-1);
    }
    
    for (i=0; i < numclusters ; i++) {
	double complex_charge = ZERO;
	double atom_charge = ZERO;
	double receptor_charge = ZERO;
	double ligand_charge = ZERO;
	struct DOCKresult *dockdata = NULL;
	//struct QMresult *qmdata;
	//int cluster_size;
	done = FALSE;
	struct ClusterRep *clusterdata = NULL;
	struct ClusterRep *temp = NULL;

	// Initialize the pointers in clusterdata
	if(results->reps == NULL)
	  results->reps = InitClusterRep(results->reps);
	clusterdata = results->reps;
	while (clusterdata->next != NULL) clusterdata = clusterdata->next;
	clusterdata->index = i + 1;

	dockdata = &clusterdata->docked;
	dockdata->time = difftime(autodock_time,timestamp);

	//sprintf(cmplx_filename,"%s/%s_%03d.pdbqt",clustersDir,cmplxstub,(i+1));
	sprintf(cmplx_filename,"%s/%s_%03d.pdb",clustersDir,cmplxstub,(i+1));
	printf_verbose(verbose,"%s: Placing docked complex in %s\n",nodestr,cmplx_filename);
	complexfile = fopen(cmplx_filename,"w+");
	if (complexfile == NULL) {
	    printf("%s: ERROR - Failed to open complex structure file (%s)!\n",
		   nodestr,cmplx_filename);
	    return(-1);
	}
	rewind(complexfile);

	//sprintf(qmligand_filename,"%s/%s_%03d.pdbqt",clustersDir,filestub,(i+1));
	sprintf(qmligand_filename,"%s/%s_%03d.pdb",clustersDir,filestub,(i+1));
	printf_verbose(verbose,"%s: Placing docked ligand structure in %s\n",nodestr,qmligand_filename);
	qmligandfile = fopen(qmligand_filename,"w+");
	if (qmligandfile == NULL) {
	    printf("%s: ERROR - Failed to open ligand structure file (%s)!\n",
		   nodestr,qmligand_filename);
	    return(-1);
	}
	rewind(qmligandfile);

	// Copy the receptor information to the head of the docked complex file
	rewind(receptorfile);
	while (fgets(line,ARG_MAX,receptorfile) != NULL) {
	    if (strstr(line,"ATOM") == line ||
		strstr(line,"HETATM") == line ||
		strstr(line,"TER") == line ||
		strstr(line,"REMARK") == line ||
		strstr(line,"AUTHOR") == line) {
		fputs(line,complexfile);
	    }
	    if (strstr(line,"ATOM") == line || strstr(line,"HETATOM") == line) {
		atom_charge = ZERO;
		GetValueByColumn(line," ",12,&atom_charge);
		receptor_charge = receptor_charge + atom_charge;
	    }
	}
	fflush(complexfile);
	// Now append the docked structure to the docked complex file
	fseek(dockfile,startmodel[i],SEEK_SET);
	ligand_charge = ZERO;
	char* lineptr;
	char  tempstr[4];
	while (!done && (fgets(line,ARG_MAX,dockfile) != NULL)) {
	    if (strstr(line,MODEL_END) == line) done = TRUE;
	    else if ((strstr(line,"ATOM") == line) || (strstr(line,"TER") == line)) {
		if (strstr(line,"ATOM") == line) {
		    // Set some PDB fields explicitly
		    lineptr = line;
		    // Set the residue name for the ligand
		    lineptr = lineptr + 17;
		    strcpy(tempstr,"LIG");
		    strncpy(lineptr,tempstr,3);
		    // Set the residue number for the ligand
		    lineptr = lineptr + 5;
		    strcpy(tempstr,"   1");
		    strncpy(lineptr,tempstr,4);
		    // Keep the charge but drop the trailing information.
		    // MOPAC should be able to deal with the charge
		    lineptr = lineptr + 55;
		    strcpy(tempstr," \n\0");
		    strncpy(lineptr,tempstr,3);
		}
                // NOTE: Autodock 4.2 has a bug in its output routines.  All "Br" atoms
                //       are written as simply "r" in the .dlg file.  This needs to be
                //       fixed for MOPAC.
	        if (strstr(line,"ATOM") == line) {
                    lineptr = line + 13;
                    if (strstr(lineptr,"r") == lineptr) {
                        lineptr--;
                        strncpy(lineptr,"B",1);
                    }
                }
                // NOTE: Autodock 4.2 has a bug in its output routines.  All "Cl" atoms
                //       are written as simply "r" in the .dlg file.  This needs to be
                //       fixed for MOPAC.
	        if (strstr(line,"ATOM") == line) {
                    lineptr = line + 13;
                    if (strstr(lineptr,"l") == lineptr) {
                        lineptr--;
                        strncpy(lineptr,"C",1);
                    }
                }
		fputs(line,complexfile);
		fputs(line,qmligandfile);
	    }

	    // Extract the results into the dockjob data structure
	    // Get the number of poses in this cluster
	    if (strstr(line,CLUSTER_SIZE_MARKER) != NULL) {
		// Here atom_charge is used as a temporary variable
		GetValueByColumn(line," ",9,&atom_charge);
		clusterdata->size  = (int) atom_charge;
	    }
	    // Get the charge
	    /*
	    if (strstr(line,"ATOM") == line) {
		atom_charge = ZERO;
		GetValueByColumn(line," ",12,&atom_charge);
		ligand_charge = ligand_charge + atom_charge;
	    }
            */
            // Let's let AutoDock calculate this for us
            if (strstr(line,DOCK_LIGAND_CHARGE) != NULL) {
                GetValueByColumn(line," ",6,&ligand_charge);
            }
	    // Get the free energy of binding
	    if (strstr(line,G_binding_MARKER) != NULL) 
		GetValueByColumn(line," ",8,&dockdata->G_binding);
	    // Get the Ki
	    if (strstr(line,Ki_EST_MARKER) != NULL) {
		GetValueByColumn(line," ",7,&dockdata->Ki_DOCK);
		GetStringByColumn(line," ",8,&dockdata->Ki_unit);
	    }
	    // Get the interaction energy
	    if (strstr(line,E_INTER_MARKER) != NULL)
		GetValueByColumn(line," ",7,&dockdata->E_inter);
	    // Get the nonbonded interaction energy
	    if (strstr(line,E_NONBOND_MARKER) != NULL)
		GetValueByColumn(line," ",9,&dockdata->E_nonbond);
	    // Get the electrostatic energy
	    if (strstr(line,E_ELECTROSTAT_MARKER) != NULL)
		GetValueByColumn(line," ",5,&dockdata->E_electrostat);
	    // Get the internal energy
	    if (strstr(line,E_INTERNAL_MARKER) != NULL)
		GetValueByColumn(line," ",8,&dockdata->E_internal);
	    // Get the torsional free energy
	    if (strstr(line,G_TORS_MARKER) != NULL)
		GetValueByColumn(line," ",7,&dockdata->G_tors);
	    // Get the unbound energy
	    if (strstr(line,E_UNBOUND_MARKER) != NULL)
		GetValueByColumn(line," ",7,&dockdata->E_unbound);
	    // Get the RMSD deviation
	    if (strstr(line,DOCK_RMSD_MARKER) != NULL)
		GetValueByColumn(line," ",7,&dockdata->rmsd_ref);


	    memset(line,(char) NULL,ARG_MAX);
	}

	fsync(fileno(complexfile));
	if (fclose(complexfile) == -1) {
	    printf("%s: ERROR - Failed to close complex file (%s)\n",nodestr,cmplx_filename);
	    printf("%s: ERROR - %s\n",nodestr,strerror(errno));
	    retval = -1;
	    return(retval);
	}

	fsync(fileno(qmligandfile));
	if (fclose(qmligandfile) == -1) {
	    printf("%s: ERROR - Failed to close ligand file (%s)\n",nodestr,qmligand_filename);
	    printf("%s: ERROR - %s\n",nodestr,strerror(errno));
	    retval = -1;
	    return(retval);
	}

	const double ROUNDING_FACTOR = 10.0;
	receptor_charge = ((double) floor(receptor_charge*ROUNDING_FACTOR+0.5))/ROUNDING_FACTOR;
	ligand_charge = ((double) floor(ligand_charge*ROUNDING_FACTOR+0.5))/ROUNDING_FACTOR;
	complex_charge = receptor_charge + ligand_charge;
	results->charge = ligand_charge;
	results->total_charge = complex_charge;
        #ifdef DEBUG
       	printf("%s: DEBUG - CHARGES receptor = %6.3f  ligand = %6.3f  total = %6.3f\n",
	       nodestr,receptor_charge,ligand_charge,complex_charge);
        #endif

	printf("\n");

	if ((i+1) < numclusters) {
	    temp = InitClusterRep(temp);
	    clusterdata->next = temp;
	}
    }

    fsync(fileno(dockfile));
    if (fclose(dockfile) == -1) {
	printf("%s: ERROR - Failed to close docking results file (%s)\n",nodestr,dlg_file);
	printf("%s: ERROR - %s\n",nodestr,strerror(errno));
	retval = -1;
	return(retval);
    }
    fsync(fileno(receptorfile));
    if (fclose(receptorfile) == -1) {
	printf("%s: ERROR - Failed to close receptor file (%s)\n",nodestr,receptPDB);
	printf("%s: ERROR - %s\n",nodestr,strerror(errno));
	retval = -1;
	return(retval);
    }


    // Move results to the right place
    // - Store the docking log file
    sprintf(tmp_name,"%s/%s",resultsDir,dlg_file);
    MoveFile(dlg_file,tmp_name);
    
    if (verbose) {
	for (i=1;i<OUTPUT_WIDTH;i++) printf_verbose(verbose,"-"); 
	printf_verbose(verbose,"\n");
    }

    first_job = FALSE;
    chdir(initdir);

    // Use fsync in an attempt to flush the file handles for deleted files
    srcfile = fopen(workdir,"rw");
    fsync(fileno(srcfile));
    fclose(srcfile);
    destfile = fopen(resultsDir,"rw");
    fsync(fileno(destfile));
    fclose(destfile);

    // We don't need MGL_args and MGL_env any more
    for (i=0;i<num_MGL_args;i++) {
	free(MGL_args[i]);
	MGL_args[i] = NULL;
    }
    free(MGL_args);
    MGL_args = NULL;

    for (i=0;i<num_MGL_env;i++) {
	free(MGL_env[i]);
	MGL_env[i] = NULL;
    }
    free(MGL_env);
    MGL_env = NULL;

    return(retval);

}
