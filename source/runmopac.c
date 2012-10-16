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
 * This routine takes the name of a .pdbqt file (without the file        *
 * extension), generates an appropriate .dat file for a MOPAC2009 job,   *
 * runs the job and then extracts the results into a STICelement, which  *
 * was provided by the calling routine.                                  *
 *************************************************************************/

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
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

int MOPACException(int retval,struct STICelement *MOPAC_STIC, const char *format, ...)
{
  va_list args;
  int S,T,I,C;
  struct ClusterRep *reps;
  va_start(args, format);
  vprintf(format,args);
  va_end(args);
  fflush(stdout);
  
  if(MOPAC_STIC != NULL)
    {
      S = MOPAC_STIC->S;T = MOPAC_STIC->T;I = MOPAC_STIC->I;C = MOPAC_STIC->C;
      reps = MOPAC_STIC->reps;
      SetupSTIC(MOPAC_STIC);
      MOPAC_STIC->S = S;MOPAC_STIC->T = T;MOPAC_STIC->I = I;MOPAC_STIC->C = C;
      MOPAC_STIC->reps = reps;
    }

  return retval;
}

void MOPAC_add_keywords(char *deck,char **params, int num_params)
{
  int param_counter = 0;
  if(params == NULL || num_params == 0 || deck == NULL)
    return;
  for(;param_counter < num_params;param_counter++)
    if(params[param_counter] != NULL)
      {
	strcat(deck," ");
	strcat(deck,params[param_counter]);
	strcat(deck," ");
      }
  
  strcat(deck,"+\n");
}

int RunMOPAC(job_t *job,
	     JobParameters *parameters,
	     struct STICelement *MOPAC_STIC)
{
  struct QMresult *MOPACresults;
  double *STIC_charge;
  int retval = 0;
  char initdir[FILENAME_MAX];
  char workdir[FILENAME_MAX];
  char filestub[FILENAME_MAX];
  char pdbqtfile_name[FILENAME_MAX];
  char temp_name[FILENAME_MAX];
  char src_name[FILENAME_MAX];
  char dest_name[FILENAME_MAX];
  char testdeck_name[FILENAME_MAX];
  char testout_name[FILENAME_MAX];
  char inputdeck_name[FILENAME_MAX];
  char outputdeck_name[FILENAME_MAX];
  char line[ARG_MAX];
  char cmdstr[ARG_MAX];
  double atom_charge, total_charge, MOPAC_charge;
  FILE *pdbqtfile;
  FILE *testdeck;
  FILE *testout;
  FILE *inputdeck;
  FILE *outputdeck;
  char **MOPAC_args = NULL;
  int num_MOPAC_args = 0;
  char **MOPAC_env = NULL;
  int num_MOPAC_env = 0;
  int charge_OK = FALSE;
  int should_restart = FALSE;// Indicates whether or not a restart is taking place.
  int i;
  pid_t pid;

  // Parameters pass from Effercio.c
  char *structDir;
  char *resultsDir = parameters->optimized_dir;// Null checked in Effercio.c
  char *method = parameters->qm_method;
  int useMOZYME = parameters->use_mozyme;
  int doFreeEnergy = parameters->UseFreeEnergy;
  char *receptor = parameters->receptor_name;
  char *user = parameters->username;
  char *nodestr = parameters->node_tag;
  char *scratch_dir = parameters->scratch_dir;
  int verbose = parameters->verbose;
  double time_left = parameters->total_time - time(NULL) + parameters->start_time;
  time_left -= 180;// Provided a buffer, in which restart files may be sent to the master node.

  char MOPAC_header[512];
  char MOPAC_footer[512];

  int num_underscores = 0;
  int IncludesReceptor = FALSE;
  char *find;

  for (i=0; i < strlen(job->name); i++) if (job->name[i] == '_') num_underscores++;

  find = job->name;
  if (strstr(job->name,receptor) == job->name) {
    IncludesReceptor = TRUE;
    find = &job->name[strlen(receptor)];
    if(strlen(find) == 0)
      return MOPACException(MOPAC_ERROR,"%s: ERROR - No ligand name and/or STIC indices in mopac job named %s.\n",nodestr,job->name);
    if(find[0] == '_')
      find = &find[1];
  }
  //MOPAC_STIC->S = ExtractIndex(job,'S');
  //MOPAC_STIC->T = ExtractIndex(job,'T');
  //MOPAC_STIC->I = ExtractIndex(job,'I');
  //MOPAC_STIC->C = ExtractIndex(job,'C');
  MOPAC_STIC->S = ExtractIndex(find,'S');
  MOPAC_STIC->T = ExtractIndex(find,'T');
  MOPAC_STIC->I = ExtractIndex(find,'I');
  MOPAC_STIC->C = ExtractIndex(find,'C');

#ifdef DEBUG
  printf("DEBUG - %s contains %d '_' characters\n",job->name,num_underscores);
  printf("DEBUG - receptor = %s\n",receptor);
  printf("DEBUG - S = %d\tT = %d\tI = %d\tC = %d\n",MOPAC_STIC->S,MOPAC_STIC->T,
	 MOPAC_STIC->I,MOPAC_STIC->C);
#endif

  if(job->type == PRESCREENING)
    structDir = parameters->ligand_dir;
  else
    structDir = parameters->clusters_dir;

  if(MOPAC_STIC->reps == NULL)
    MOPAC_STIC->reps = InitClusterRep(MOPAC_STIC->reps);
  MOPACresults = &MOPAC_STIC->reps->optimized;
  STIC_charge = &MOPAC_STIC->charge;
   

  // Initialize arrays
  memset(initdir,(char) NULL,FILENAME_MAX);
  memset(workdir,(char) NULL,FILENAME_MAX);
  memset(filestub,(char) NULL,FILENAME_MAX);
  memset(temp_name, (char) NULL,FILENAME_MAX);
  memset(src_name, (char) NULL,FILENAME_MAX);
  memset(dest_name, (char) NULL,FILENAME_MAX);
  memset(pdbqtfile_name,(char) NULL,FILENAME_MAX);
  memset(testdeck_name,(char) NULL,FILENAME_MAX);
  memset(testout_name,(char) NULL,FILENAME_MAX);
  memset(inputdeck_name,(char) NULL,FILENAME_MAX);
  memset(outputdeck_name,(char) NULL,FILENAME_MAX);
  memset(line,(char) NULL,ARG_MAX);
  memset(cmdstr,(char) NULL,ARG_MAX);

  // Set up the filenames
  getcwd(initdir,FILENAME_MAX);

  strcpy(workdir,scratch_dir);
  strcat(workdir,"/");
  strcat(workdir,user);

  strcpy(filestub,job->name);
    
  strcpy(pdbqtfile_name,filestub);
  strcat(pdbqtfile_name,".pdbqt");

  strcpy(testdeck_name,filestub);
  strcat(testdeck_name,"_test.dat");
    
  strcpy(testout_name,filestub);
  strcat(testout_name,"_test.out");

  strcpy(inputdeck_name,filestub);
  strcat(inputdeck_name,".");
  strcat(inputdeck_name,method);
  strcat(inputdeck_name,".dat");

  strcpy(outputdeck_name,filestub);
  strcat(outputdeck_name,".");
  strcat(outputdeck_name,method);
  strcat(outputdeck_name,".out");

  // Write out the starting values of everything for DEBUG purposes
#ifdef DEBUG
  printf("%s: DEBUG - initdir = %s\n",nodestr,initdir);
  printf("%s: DEBUG - workdir = %s\n",nodestr,workdir);
  printf("%s: DEBUG - filestub = %s\n",nodestr,filestub);
  printf("%s: DEBUG - pdbqtfile_name = %s\n",nodestr,pdbqtfile_name);
  printf("%s: DEBUG - inputdeck_name = %s\n",nodestr,inputdeck_name);
  printf("%s: DEBUG - user = %s\n",nodestr,user);
  printf("%s: DEBUG - STIC_charge = %.2f\n",nodestr,*STIC_charge);
  printf("%s: DEBUG - method = %s\n",nodestr,method);
  printf("%s: DEBUG - useMOZYME = %d\n",nodestr,useMOZYME);
  printf("%s: DEBUG - doFreeEnergy = %d\n",nodestr,doFreeEnergy);
#endif

    
  // Check for working directory
  if (VerifyDir(workdir,verbose,nodestr) == 0) {
    chdir(workdir);
  }
  else {
    printf("Cannot change to directory %s",workdir);
    exit(MOPAC_ERROR);
  }

  // Bring the .pdbqt file to the local scratch directory
  sprintf(temp_name,"%s/%s",structDir,pdbqtfile_name);
  if (access(temp_name,F_OK) == -1) {
    memset(pdbqtfile_name,(char) NULL,FILENAME_MAX);
    strcpy(pdbqtfile_name,filestub);
    strcat(pdbqtfile_name,".pdb");
    sprintf(temp_name,"%s/%s",structDir,pdbqtfile_name);
  }
  CopyFile(temp_name,pdbqtfile_name);
#ifdef DEBUG
  printf("DEBUG - %s: Copied %s\n",nodestr,temp_name);
  printf("DEBUG - %s:     to %s/%s\n",nodestr,workdir,pdbqtfile_name);
#endif

  // Open the .pdbqt file.  We'll need this for a bunch of things later
  pdbqtfile = fopen(pdbqtfile_name,"r");
  if (pdbqtfile == NULL) {
    printf("%s: ERROR - Failed to open MOPAC input deck for writing (%s)\n",nodestr,pdbqtfile_name);
    return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
  }

  // Get the charge from the .pdbqt file if none is supplied
#ifdef DEBUG
  printf("%s: DEBUG - supplied total_charge is %.1f\n",nodestr,*STIC_charge);
#endif
  if (*STIC_charge == -9876543.21) {
    // Need to calculate the charge for the MOPAC input deck
    total_charge = ZERO;
    rewind(pdbqtfile);
    while (fgets(line,ARG_MAX,pdbqtfile) != NULL) {
      if (strstr(line,"ATOM") == line) {
	atom_charge = ZERO;
	GetValueByColumn(line," ",12,&atom_charge);
	total_charge = total_charge + atom_charge;
      }
    }
    rewind(pdbqtfile);
    *STIC_charge = total_charge;
  }
  else {
    // Use the supplied charge
    total_charge = *STIC_charge;
  }

  if(strlen(job->mopac_res) != 0)
    {
      struct stat res_status;

      if(stat(job->mopac_res,&res_status) == -1)
    	{
	  fprintf(stderr,"%s: RunMopac could not stat the restart file. Therefore, the MOPAC run will start over.\n",nodestr);
	  fprintf(stderr,"%s: Reason: %s",nodestr,strerror(errno));
    	}
      else
    	{
	  char restart_filename[FILENAME_MAX];
	  sprintf(restart_filename,"%s.%s.res",filestub,method);
	  CopyFile(job->mopac_res,restart_filename);
	  should_restart = TRUE;
    	}
    }

#ifdef DEBUG
  printf("%s: DEBUG - Got molecular charge of %.2f\n",nodestr,total_charge);
#endif

  // Do some error-correcting here in case the empirical charge models got it wrong
  // In particular CORINA seems to generate structures with serious problems
  while (!charge_OK) {
#ifdef DEBUG
    printf("%s: DEBUG - total charge in test loop is %.3f\n",nodestr,total_charge);
    printf("%s: DEBUG - MOPAC charge in test loop is %.3f\n",nodestr,MOPAC_charge);
    printf("%s: DEBUG - Charge difference = %.4f\t\tTolerance = %.4f\n",nodestr,
	   fabs(MOPAC_charge - total_charge),CHARGE_TOL);
#endif
    testdeck = fopen(testdeck_name,"w+");
    if (testdeck == NULL) {
      printf("%s: ERROR - Failed to open MOPAC input deck for writing (%s)\n",nodestr,testdeck_name);
      return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
    }
    memset(MOPAC_header,(char) NULL, 512);
    // Added user specified header keywords
    MOPAC_add_keywords(MOPAC_header,parameters->mopac_header_params, parameters->num_mopac_header_params);
    sprintf(line," %s ",method);
    strcat(MOPAC_header,line);
    strcat(MOPAC_header,"LET CHARGES ");
    if (useMOZYME) {
      strcat(MOPAC_header,"MOZYME ");
    }
    strcat(MOPAC_header,MOPAC_KEYWDS_SOLV);
    sprintf(line," CHARGE=%.0f\n %s: Test run to validate charges\n\n",total_charge,job->name);
    strcat(MOPAC_header,line);
    memset(line,(char) NULL,ARG_MAX);
    fputs(MOPAC_header,testdeck);


    rewind(pdbqtfile);
    while (fgets(line,ARG_MAX,pdbqtfile) != NULL) {
      const int LINE_OFFSET = 76;
      if (strstr(line,"ATOM") == line || strstr(line,"HETATM") == line ||
	  strstr(line,"TER") == line) {
	int line_length = 0;
	line_length = strlen(line);
	char *lineptr;
	lineptr = line;
	lineptr += LINE_OFFSET;
	strncpy(lineptr,"\n",1);
	int line_counter;
	for (line_counter=LINE_OFFSET+1;line_counter<line_length;line_counter++) {
	  lineptr++;
	  lineptr[0] = (char) NULL;
	}
	fputs(line,testdeck);
      }
    }
    fsync(fileno(testdeck));
    if (fclose(testdeck) == -1) {
      printf("%s: ERROR- Failed to close MOPAC input deck (%s)\n",nodestr,testdeck_name);
      return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
    }
    // Create argv and environmental variables for MOPAC
    add_strtoarray(&MOPAC_args,&num_MOPAC_args,MOPAC_EXE);

    // Set the name of the input deck
    add_strtoarray(&MOPAC_args,&num_MOPAC_args,testdeck_name);

    // Dump the argument list for debug purposes
#ifdef DEBUG
    printf("%s: DEBUG - MOPAC_args is:\n",nodestr);
    for (i=0;i<num_MOPAC_args;i++) {
      printf("%s: DEBUG - MOPAC_args[%d] = %s\n",nodestr,i,MOPAC_args[i]);
    }
#endif

    if (strlen(envMOPAC_LIC) > 0)
      {
	sprintf(cmdstr,"MOPAC_LICENSE=%s",envMOPAC_LIC);
	add_strtoarray(&MOPAC_env,&num_MOPAC_env,cmdstr);
      }
    add_strtoarray(&MOPAC_env,&num_MOPAC_env,envPATH);
    sprintf(cmdstr,"PWD=%s",workdir);
    add_strtoarray(&MOPAC_env,&num_MOPAC_env,cmdstr);

    // Dump the environment list for DEBUG purposes
#ifdef DEBUG
    printf("%s: DEBUG - MOPAC_env is:\n",nodestr);
    for (i=0;i<num_MOPAC_env;i++) {
      printf("%s: DEBUG - MOPAC_env[%d] = %s\n",nodestr,i,MOPAC_env[i]);
    }
#endif

    // Run MOPAC
    sprintf(cmdstr,"%s/%s",MOPAC_HOME,MOPAC_EXE);
    if (verbose) {
      printf("%s: Charge validation step for %s (using %s)\n",nodestr,filestub,MOPAC_EXE);
      fflush(stdout);
    }
#ifdef DEBUG
    printf("%s: DEBUG - cmdstr = %s\n",nodestr,cmdstr);
#endif
    pid = vfork();
    if (pid == -1) {
      // fork error code here
      printf("%s: ERROR - Failed to initialize child process via vfork()\n",nodestr);
      return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
    }
    if (pid == 0) {
      if(execve(cmdstr,MOPAC_args,MOPAC_env))
	{
	  return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - execve failed for test deck.\nReason given:\n%s\n",nodestr,strerror(errno));
	}
    }
    else {
      int testdeck_retval;
      wait(&testdeck_retval);
      if(testdeck_retval != 0)
	return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - execve failed for test deck. See log for details.\n",nodestr);
	  
    }
    if (verbose) {
      printf("%s: %s has finished\n",nodestr,MOPAC_EXE);
      fflush(stdout);
    }

    testout = fopen(testout_name,"r");
    if (testout == NULL) {
      printf("%s: ERROR - Failed to open MOPAC input deck for reading (%s)\n",nodestr,testout_name);
      return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
    }
	
    rewind(testout);
    while(fgets(line,ARG_MAX,testout) != NULL) {
      if (strstr(line,MOPAC_GEOMETRY_ERROR) != NULL) {
	printf("%s: ERROR - The calculation for %s was aborted.\n",nodestr,filestub);
	printf("%s: ERROR - MOPAC has discovered a problem with the input geometry!\n",nodestr);
	return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - Repair the geometry for %s and run the calculation again.\n",nodestr,filestub);
      }
      if (strstr(line,MOPAC_CALC_CHARGE) != NULL) {
	GetValueByColumn(line," ",5,&MOPAC_charge);
#ifdef DEBUG
	printf("%s: DEBUG - total charge after run in test loop is %.3f\n",nodestr,total_charge);
	printf("%s: DEBUG - MOPAC charge after run in test loop is %.3f\n",nodestr,MOPAC_charge);
	printf("%s: DEBUG - Charge difference = %.4f\t\tTolerance = %.4f\n",nodestr,
	       fabs(MOPAC_charge - total_charge),CHARGE_TOL);
#endif
	if (fabs(MOPAC_charge - total_charge) < CHARGE_TOL ) {
	  charge_OK = TRUE;
	  break;
	}
	else {
	  printf("%s: WARNING - Problem with job %s\n",nodestr,job->name);
	  printf("%s: WARNING - The total charge for your system is incorrect!!\n",nodestr);
	  printf("%s: WARNING - You specified a charge of %.1f and MOPAC thinks this should be %.1f\n",
		 nodestr,total_charge,MOPAC_charge);
	  printf("%s: WARNING - %s will use the MOPAC charge\n",nodestr,PROG_NAME);
	  printf("%s: WARNING - Please check your results and input files carefully\n",nodestr);
	  total_charge = MOPAC_charge;
	  break;
	}
      }
    }

    fsync(fileno(testout));
    if (fclose(testout) == -1) {
      printf("%s: ERROR- Failed to close MOPAC input deck (%s)\n",nodestr,testout_name);
      return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
    }
	
    FreeStringArray(&MOPAC_env,&num_MOPAC_env);
    FreeStringArray(&MOPAC_args,&num_MOPAC_args);

#ifndef DEBUG
    remove(testdeck_name);
    remove(testout_name);
    sprintf(temp_name,"%s_test.arc",filestub);
    remove(temp_name);
#endif
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Now that we've verified that the charge is correct ... do the real work

  inputdeck = fopen(inputdeck_name,"w+");
  if (inputdeck == NULL) {
    printf("%s: ERROR - Failed to open MOPAC input deck for writing (%s)\n",nodestr,inputdeck_name);
    return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
  }
    
  // Set up the header lines for the MOPAC input deck
  memset(MOPAC_header,(char) 0, 512);
  memset(MOPAC_footer,(char) 0, 512);
  MOPAC_add_keywords(MOPAC_header,parameters->mopac_header_params,parameters->num_mopac_header_params);
  MOPAC_add_keywords(MOPAC_footer,parameters->mopac_footer_params,parameters->num_mopac_footer_params);
    
  sprintf(line," %s ",method);
  strcat(MOPAC_header,line);
  if (verbose) printf("%s: DEBUG - Time remaining for MOPAC calculation is %.0f seconds\n",nodestr,time_left);
#ifdef DEBUG
  time_left = (double) MOPAC_DEBUG_TIME;
  printf("%s: DEBUG - Using remaining time of %.0f instead!!!!\n",nodestr,time_left);
#endif
  sprintf(line,"T=%.0f DUMP=4320 ",time_left);
  strcat(MOPAC_header,line);
  //strcat(MOPAC_header,"XYZ GRAPH ");
  strcat(MOPAC_header,"XYZ ");
  strcat(MOPAC_header, MOPAC_KEYWDS_SOLV);
  //strcat(MOPAC_header," MULLIK ");
  sprintf(line," %s ",MOPAC_KEYWDS_CHRG);
  strcat(MOPAC_header,line);

  if(should_restart)
    {
      printf_verbose(parameters->verbose,"%s: Restarting %s\n",nodestr,job->name);
      strcat(MOPAC_header," RESTART ");
    }

  if (useMOZYME) {
    strcat(MOPAC_header,"LBFGS MOZYME PDBOUT +\n ");
  }
  else {
    //strcat(MOPAC_header,"BFGS SYBYL +\n ");
    strcat(MOPAC_header,"SYBYL +\n ");
  }

  sprintf(line,"CHARGE=%.0f",total_charge);
  strcat(MOPAC_header,line);
  //strcat(MOPAC_header," PRECISE LET DDMIN=0.0 GNORM=0.15");
  strcat(MOPAC_header," RELSCF=0.1 LET DDMIN=0.0 GNORM=0.15");

#ifdef DEBUG
  strcat(MOPAC_header," NOOPT");
#endif
  // If this is only a ligand in a docked pose, don't optimize
  if(job->type == QMLIG) {
    strcat(MOPAC_header," 1SCF");
  }

  // Add the comment line for the job
  strcat(MOPAC_header,"\n");
  strcat(MOPAC_header,job->name);
  strcat(MOPAC_header,": Geometry optimization step\n\n");

  sprintf(line," %s ",method);
  strcat(MOPAC_footer,line);
  sprintf(line,"T=%.0f ",time_left*FREE_ENERGY_TIME_RATIO);
  if(time_left > 2700)// If there are only 45 minutes left (or less) dump the restart file more frequently
    sprintf(line,"DUMP=4320 ");
  else
    sprintf(line,"DUMP=240 ");
  strcat(MOPAC_footer,line);
  if (useMOZYME) {
    strcat(MOPAC_footer,"MOZYME ");
  }
  //strcat(MOPAC_footer,"PRECISE THERMO LET ");
  //strcat(MOPAC_footer,"PRECISE THERMO ");
  strcat(MOPAC_footer,"RELSCF=0.1 THERMO ");
  strcat(MOPAC_footer,MOPAC_KEYWDS_SOLV);
  strcat(MOPAC_footer," +\n OLDGEO GEO-OK ");
  sprintf(line,"CHARGE=%.0f",total_charge);
  strcat(MOPAC_footer,line);
#ifdef DEBUG
  strcat(MOPAC_footer," LET");
#endif
  strcat(MOPAC_footer,"\n ");
  strcat(MOPAC_footer,job->name);
  strcat(MOPAC_footer,": Free energy calculation\n");

  if (!doFreeEnergy) MOPAC_footer[0] = (char) NULL;

#ifdef DEBUG
  printf("%s: DEBUG - MOPAC_header\n%s\n%s\n%s\n",nodestr,DASH_HDR,MOPAC_header,DASH_HDR);
  printf("\n%s: DEBUG - MOPAC_footer\n%s\n%s\n%s\n",nodestr,DASH_HDR,MOPAC_footer,DASH_HDR);
#endif
    
  // Write the MOPAC header lines
  fputs(MOPAC_header,inputdeck);

  // Write the geometry
  rewind(pdbqtfile);
  while (fgets(line,ARG_MAX,pdbqtfile) != NULL) {
    const int LINE_OFFSET = 76;
    if (strstr(line,"ATOM") == line || strstr(line,"HETATM") == line ||
	strstr(line,"TER") == line) {
      int line_length = 0;
      line_length = strlen(line);
      char *lineptr;
      lineptr = line;
      lineptr += LINE_OFFSET;
      strncpy(lineptr,"\n",1);
      int line_counter;
      for (line_counter=LINE_OFFSET+1;line_counter<line_length;line_counter++) {
	lineptr++;
	lineptr[0] = (char) NULL;
      }
      fputs(line,inputdeck);
    }
  }

  // Write the MOPAC footer lines
  
  fputs("\n",inputdeck);
  fputs(MOPAC_footer,inputdeck);

  fsync(fileno(pdbqtfile));
  if (fclose(pdbqtfile) == -1) {
    printf("%s: ERROR- Failed to close MOPAC input deck (%s)\n",nodestr,pdbqtfile_name);
    return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
  }

  fsync(fileno(inputdeck));
  if (fclose(inputdeck) == -1) {
    printf("%s: ERROR- Failed to close MOPAC input deck (%s)\n",nodestr,inputdeck_name);
    return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
  }

  // Create argv and environmental variables for MOPAC
  // The first string in the argument list should be the executable name
  add_strtoarray(&MOPAC_args,&num_MOPAC_args,MOPAC_EXE);

  // Set the name of the input deck
  add_strtoarray(&MOPAC_args,&num_MOPAC_args,inputdeck_name);

  // Dump the argument list for debug purposes
#ifdef DEBUG
  printf("%s: DEBUG - MOPAC_args is:\n",nodestr);
  for (i=0;i<num_MOPAC_args;i++) {
    printf("%s: DEBUG - MOPAC_args[%d] = %s\n",nodestr,i,MOPAC_args[i]);
  }
#endif

  if (strlen(envMOPAC_LIC) > 0) {
    sprintf(cmdstr,"MOPAC_LICENSE=%s",envMOPAC_LIC);
    add_strtoarray(&MOPAC_env,&num_MOPAC_env,cmdstr);
  }
  add_strtoarray(&MOPAC_env,&num_MOPAC_env,envPATH);
  sprintf(cmdstr,"PWD=%s",getenv("PWD"));
  add_strtoarray(&MOPAC_env,&num_MOPAC_env,cmdstr);

  // Dump the environment list for DEBUG purposes
#ifdef DEBUG
  printf("%s: DEBUG - MOPAC_env is:\n",nodestr);
  for (i=0;i<num_MOPAC_env;i++) {
    printf("%s: DEBUG - MOPAC_env[%d] = %s\n",nodestr,i,MOPAC_env[i]);
  }
#endif

  // Run MOPAC
  sprintf(cmdstr,"%s/%s",MOPAC_HOME,MOPAC_EXE);
  if (verbose) {
    printf("%s: Semi-empirical QM step for %s (using %s)\n",nodestr,filestub,MOPAC_EXE);
    fflush(stdout);
  }
#ifdef DEBUG
  printf("%s: DEBUG - cmdstr = %s\n",nodestr,cmdstr);
#endif
  pid = vfork();
  if (pid == -1) {
    // fork error code here
    printf("%s: ERROR - Failed to initialize child process via vfork()\n",nodestr);
    return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
  }
  if (pid == 0) {
    execve(cmdstr,MOPAC_args,MOPAC_env);
  }
  else {
    wait(NULL);
  }
  if (verbose) {
    printf("%s: %s has finished\n",nodestr,MOPAC_EXE);
    fflush(stdout);
  }
    
  // Parse output
  memset(job->mopac_res,0,FILENAME_MAX);// this will ensure that IFF a restart file is to be used, it is indicated by the fact that job->mopac_res has a positive string length.
  outputdeck = fopen(outputdeck_name,"r");
  if (outputdeck == NULL) {
    printf("%s: ERROR - Failed to open MOPAC output deck for reading (%s)\n",nodestr,outputdeck_name);
    return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
  }
  rewind(outputdeck);

  memset(line,(char) NULL,ARG_MAX);
  int num_runs = 1;
  int run_count = 0;
  double temp_val = ZERO;
  if (doFreeEnergy) num_runs = 2;
  while(fgets(line,ARG_MAX,outputdeck) != NULL) {
    if (strstr(line,MOPAC_EXCESS_CYCLES) != NULL) {
      printf("%s: ERROR - MOPAC calculation failed from too many optimization steps!\n",nodestr);
      printf("%s: ERROR - results for %s are invalid\n",nodestr,filestub);
      retval = MOPAC_ERROR;
    }
    if (strstr(line,MOPAC_LAMBDA_EXCESS) != NULL) {
      printf("%s: ERROR - MOPAC calculation failed from problems in the geometry line search!\n",nodestr);
      printf("%s: ERROR - results for %s are invalid\n",nodestr,filestub);
      retval = MOPAC_ERROR;
    }
    if (strstr(line,MOPAC_GEOMETRY_ERROR) != NULL) {
      printf("%s: ERROR - The calculation for %s was aborted.\n",nodestr,filestub);
      printf("%s: ERROR - MOPAC has discovered a problem with the input geometry!\n",nodestr);
      printf("%s: ERROR - Repair the geometry for %s and run the calculation again.\n",nodestr,filestub);
      retval = MOPAC_ERROR;
    }
    if (strstr(line,MOPAC_GEOM_UNREC) != NULL) {
      printf("%s: ERROR - The calculation for %s failed.\n",nodestr,filestub);
      return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - MOPAC was unable to recognize the input geometry!\n",nodestr);
    }
    if (strstr(line,MOPAC_ENTHALPY) != NULL)
      GetValueByColumn(line," ",6,&MOPACresults->Hf);
    if (strstr(line,MOPAC_VDW_AREA) != NULL)
      GetValueByColumn(line," ",6,&MOPACresults->vdW_A);
    if (strstr(line,MOPAC_E_DIELEC) != NULL)
      GetValueByColumn(line," ",4,&MOPACresults->E_dielec);
    if (strstr(line,MOPAC_COSMO_AREA) != NULL)
      GetValueByColumn(line," ",4,&MOPACresults->COSMO_A);
    if (strstr(line,MOPAC_COSMO_VOL) != NULL)
      GetValueByColumn(line," ",4,&MOPACresults->COSMO_V);
    if (strstr(line,MOPAC_numSCF) != NULL) {
      GetValueByColumn(line," ",4,&temp_val);
      MOPACresults->num_SCFs = (int) temp_val;
    }
    if (strstr(line,MOPAC_DIPOLE) != NULL) {
      for (i=0; i<3; i++) {
	if (fgets(line,ARG_MAX,outputdeck) == NULL) {
	  printf("%s: ERROR - Problem reading dipole section of MOPAC outputdeck\n",nodestr);
	  return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,outputdeck_name);
	}
      }
      GetValueByColumn(line," ",2,&MOPACresults->mu_x);
      GetValueByColumn(line," ",3,&MOPACresults->mu_y);
      GetValueByColumn(line," ",4,&MOPACresults->mu_z);
      GetValueByColumn(line," ",5,&MOPACresults->mu_total);
    }
    if (strstr(line,MOPAC_TIME) != NULL && run_count <= num_runs) {
      GetValueByColumn(line," ",4,&temp_val);
      if (MOPACresults->time >= ZERO)
	MOPACresults->time += temp_val;
      else
	MOPACresults->time = temp_val;
    }
    if (strstr(line,MOPAC_SUCCESSFUL) != NULL) {
      run_count++;
      if (run_count == 1) {
	int numchars;
	numchars = strlen(method) + 1;
	MOPACresults->method = realloc(MOPACresults->method,numchars * sizeof(char));
	memset(MOPACresults->method,(char) NULL,numchars);
	strcpy(MOPACresults->method,method);
      }
    }
    if(strstr(line,MOPAC_NO_TIME_LEFT) != NULL)
      {
	sprintf(temp_name,"%s.%s.res",filestub,method);
	sprintf(job->mopac_res,"%s/%s",resultsDir,temp_name);
	printf_verbose(parameters->verbose,"%s: Saving restart file %s as %s\n",nodestr,temp_name,job->mopac_res);

	CopyFile(temp_name,job->mopac_res);
      }
    if (strstr(line,MOPAC_ZPE_MARKER) != NULL)
      GetValueByColumn(line," ",4,&MOPACresults->ZPE);
    if (strstr(line,MOPAC_THERMO_MARKER) != NULL) {
      // This in here to be able to add functionality for using a user-specified temperature
      // or temperatures in the future.
      double TargetTemp = 298.00E0;
      int foundTemp = FALSE;
      double TestTemp = ZERO;
      while(fgets(line,ARG_MAX,outputdeck) != NULL) {
	GetValueByColumn(line," ",1,&TestTemp);
	if (TestTemp == TargetTemp) {
	  foundTemp = TRUE;
	  for (i=0; i<4; i++) fgets(line,ARG_MAX,outputdeck);
	  GetValueByColumn(line," ",2,&temp_val); // This is the enthalpy at the TargetTemp
	  GetValueByColumn(line," ",4,&MOPACresults->Cp);
	  MOPACresults->Cp = MOPACresults->Cp/1.0E03;
	  GetValueByColumn(line," ",5,&MOPACresults->S);
	  MOPACresults->S = MOPACresults->S/1.0E3;
	  // Calculate the free energy
	  MOPACresults->G = temp_val - TargetTemp * MOPACresults->S;
	  break;
	}
      }
      if (!foundTemp) {
	return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - Unable to find data for temperature of %.2f\n",nodestr,TargetTemp);
      }
    }
  }

  find = job->name;
  if (num_underscores > 1) {
    find = strrchr(job->name,'_');
    find++;
    MOPAC_STIC->reps->index = atoi(find);
#ifdef DEBUG
    printf("DEBUG - find = %s\t\tCluster index = %d\n",find,MOPAC_STIC->reps->index);
#endif
    if (!IncludesReceptor) {
      struct QMresult *qm_temp;
      qm_temp = &MOPAC_STIC->reps->optimized;
      if (doFreeEnergy) qm_temp->G_ligand = qm_temp->G;
      else qm_temp->G_ligand = qm_temp->Hf;

      // Now set everyone else back to their initial values
      qm_temp->Hf = DOUBLE_INIT;
      qm_temp->S = DOUBLE_INIT;
      qm_temp->Cp = DOUBLE_INIT;
      qm_temp->G = DOUBLE_INIT;
      qm_temp->ZPE = DOUBLE_INIT;
      qm_temp->E_dielec = DOUBLE_INIT;
      qm_temp->mu_x = DOUBLE_INIT;
      qm_temp->mu_y = DOUBLE_INIT;
      qm_temp->mu_z = DOUBLE_INIT;
      qm_temp->mu_total = DOUBLE_INIT;
      qm_temp->time = DOUBLE_INIT;
      qm_temp->num_SCFs = INT_INIT;
      qm_temp->COSMO_A = DOUBLE_INIT;
      qm_temp->COSMO_V = DOUBLE_INIT;
      qm_temp->vdW_A = DOUBLE_INIT;
      qm_temp->G_prot = DOUBLE_INIT;
      qm_temp->G_binding = DOUBLE_INIT;
      qm_temp->Ki_QM = DOUBLE_INIT;
    }
  }

  // Close the file ... we're done with it
  if (fclose(outputdeck) == -1) {
    printf("%s: ERROR- Failed to close MOPAC output deck (%s)\n",nodestr,outputdeck_name);
    return MOPACException(MOPAC_ERROR,MOPAC_STIC,"%s: ERROR - %s\n",nodestr,strerror(errno));
  }

  // Move results files to the appropriate place
  if (strcmp(resultsDir,scratch_dir) != 0) {
    struct stat file_status;
    int stat_retval = 0;
    // move the MOPAC output file
    sprintf(dest_name,"%s/%s.%s.out",resultsDir,filestub,method);
    CopyFile(outputdeck_name,dest_name);
    remove(outputdeck_name);
    // move the structure file
    if (useMOZYME) {
      sprintf(src_name,"%s.%s.pdb",filestub,method);
      sprintf(dest_name,"%s/%s.%s.pdb",resultsDir,filestub,method);
    }
    else {
      sprintf(src_name,"%s.%s.syb",filestub,method);
      sprintf(dest_name,"%s/%s.%s.syb",resultsDir,filestub,method);
    }
    stat_retval = stat(src_name,&file_status);
    if(stat_retval != 0 && errno == ENOENT && strlen(job->mopac_res) != 0)
      printf_verbose(parameters->verbose,"%s: Could not find structure file, but a restart has been scheduled.\n",nodestr);
    else
      {
	CopyFile(src_name,dest_name);
	remove(src_name);
      }
  }
    
  // Clean up
#ifndef DEBUG
  if (retval == 0) {
    remove(inputdeck_name);
    sprintf(temp_name,"%s.%s.arc",filestub,method);
    remove(temp_name);
    sprintf(temp_name,"%s.%s.den",filestub,method);
    remove(temp_name);
    sprintf(temp_name,"%s.%s.gpt",filestub,method);
    remove(temp_name);
    sprintf(temp_name,"%s.%s.res",filestub,method);
    remove(temp_name);
    sprintf(temp_name,"%s.%s.temp",filestub,method);
    remove(temp_name);
  }
#endif

  // Change back to the original working directory
  chdir(initdir);

  FreeStringArray(&MOPAC_env,&num_MOPAC_env);
  FreeStringArray(&MOPAC_args,&num_MOPAC_args);

  return retval;
}
