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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <pwd.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <dirent.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include "defines.h"
#include "structs.h"
#include "deque.h"
#include "RBTree.h"
#include "mpi.h"
#include "effercio_db.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "memwatch.h"

deque* CreateJobList(char inpfile[FILENAME_MAX], JobParameters *params);
char *SecondsToHMS(double num_secs);
job_t* InitJob(job_t *newjob);
int VerifyDir(char dirname[FILENAME_MAX], int verbose, char nodestr[HOST_NAME_MAX]);
void printf_verbose(int verbose, const char *format, ...);


/**
 * Pops a job from the job tree and sends it to the destination slave node.
 * The job, along with the destination node rank, is stored
 * in the busy_list stack.
 *
 * Returns pointer to the RBTree node for the job send to the slave node, destination.
 */
deque_node* issue_order(deque* jobs, deque* busy_list,int destination, int rank, JobParameters *params)
{
	job_t *job_info;
	char strJobType[512];
	double float_buffer[FLOAT_TRANSFER_COUNT], time_left;

	if(jobs == NULL)
		return NULL;
	if (destination != rank) {
		deque_node* curr_job = pop_front(jobs);

		// If there are only ten minutes (or less) left in the calculation,
		// push_back curr_job into the busy list. Then this will
		// build up a list of jobs that should be start after the next restart
		// and minimize the number of jobs that will be started without reporting
		// back to the master.
		time_left = params->total_time - time(NULL) + params->start_time;
		if(time_left < 600)
		{
            if (params->verbose)
			    printf("Only have %f seconds remaining, delaying %s\n",time_left,((job_t*)curr_job->data)->name);
			push_back(busy_list,curr_job);
			return NULL;
		}

		if(curr_job == NULL)
		{
			printf_master("WARNING - issue_order requested, but there are no jobs.");
			return NULL;
		}
		else if(curr_job->type != DEQUE_JOB)
		{
			printf_master("ERROR - Requested a job type from job tree, but the tree is not a job type tree. Type: %d",curr_job->type);
			exit(RBTREE_TYPE_ERROR);
		}
		job_info = (job_t*)curr_job->data;

        if (params->verbose) {
		    JobString(strJobType,job_info);
		    printf_master("Sending %s to Process %d for %s run\n",job_info->name,destination,strJobType);
        }

		float_buffer[0] = job_info->input_data.Hf;
		float_buffer[1] = job_info->input_data.G;
		MPI_Send(job_info->name,FILENAME_MAX,MPI_CHAR,destination,job_info->type,MPI_COMM_WORLD);
		MPI_Send(&curr_job,sizeof(RBTree*),MPI_BYTE,destination,job_info->type,MPI_COMM_WORLD);
		MPI_Send(float_buffer,FLOAT_TRANSFER_COUNT,MPI_DOUBLE,destination,job_info->type,MPI_COMM_WORLD);
		MPI_Send(job_info->mopac_res,strlen(job_info->mopac_res)+1,MPI_CHAR,destination,job_info->type,MPI_COMM_WORLD);
		job_info->host_rank = destination;

		//Put curr_job into the busy_list
		push_back(busy_list,curr_job);
        #ifdef DEBUG
		    printf_master("DEBUG - Added %s (0x%x) to busy list for %d\n",
                          ((job_t*)curr_job->data)->name,curr_job,((job_t*)curr_job->data)->host_rank);
		    fflush(stdout);
        #endif
		return curr_job;
	}
	return NULL;
}

/**
 * Adds the slave node, cpu_rank, to the wait queue, free_cpus.
 *
 * Returns a pointer to the cpunode of the suspended cpu (cpu_rank)
 */
cpunode* suspend_cpu(cpunode** free_cpus, int cpu_rank)
{
	return PushCPUNode(free_cpus, cpu_rank);
}

/**
 * Based on the job type, the function needed to run the job is
 * called. Output is stored in return_buffer. The return value
 * of the job's run function is returned.
 */
int RunJob(job_t *job,struct STICelement *stic, JobParameters *params)
{
  int retval = 0;
#ifdef FAKEDATA
	printf("%s: Ran %s\n",params->node_tag,job->name);
	stic->S = ExtractIndex(job->name,'S');
	stic->T = ExtractIndex(job->name,'T');
	stic->I = ExtractIndex(job->name,'I');
	stic->C = ExtractIndex(job->name,'C');
	stic->G = rand()/-1e9;
	stic->Hf = rand()/-1e9;
	stic->reps = InitClusterRep(stic->reps);
	if(job->type == DOCK)
	{
		stic->reps->docked.G_binding = rand()/-1e3;
		stic->reps->docked.Ki_DOCK = rand()/-1e3;
	}
	else if(job->type == QM)
		stic->reps->optimized.G = rand()/-1e3;
	else if(job->type == QMLIG)
		stic->reps->optimized.G_ligand = rand()/-1e3;
	strcpy(stic->reps->docked.Ki_unit,"mM");
	return 0;
#endif
	switch(job->type)
	{
	case PRESCREENING:
	  {
	    params->use_mozyme = FALSE;
		stic->charge = DOUBLE_INIT;
		retval = RunMOPAC(job,params,stic);
		stic->Hf = stic->reps->optimized.Hf;
		if (params->UseFreeEnergy) 
		  stic->G = stic->reps->optimized.G;
		stic->reps->optimized.Hf        = DOUBLE_INIT;
		stic->reps->optimized.S         = DOUBLE_INIT;
		stic->reps->optimized.Cp        = DOUBLE_INIT;
		stic->reps->optimized.G         = DOUBLE_INIT;
		stic->reps->optimized.ZPE       = DOUBLE_INIT;
		stic->reps->optimized.E_dielec  = DOUBLE_INIT;
		stic->reps->optimized.mu_x      = DOUBLE_INIT;
		stic->reps->optimized.mu_y      = DOUBLE_INIT;
		stic->reps->optimized.mu_z      = DOUBLE_INIT;
		stic->reps->optimized.mu_total  = DOUBLE_INIT;
		stic->reps->optimized.time      = DOUBLE_INIT;
		stic->reps->optimized.num_SCFs  = INT_INIT;
		stic->reps->optimized.COSMO_A   = DOUBLE_INIT;
		stic->reps->optimized.COSMO_V   = DOUBLE_INIT;
		stic->reps->optimized.vdW_A     = DOUBLE_INIT;
		stic->reps->optimized.G_prot    = DOUBLE_INIT;
		stic->reps->optimized.G_ligand  = DOUBLE_INIT;
		stic->reps->optimized.G_binding = DOUBLE_INIT;
		stic->reps->optimized.Ki_QM     = DOUBLE_INIT;
		break;
	  }
	case QM:
		params->use_mozyme = TRUE;
	case QMLIG:
		// stic->charge = DOUBLE_INIT;
		retval = RunMOPAC(job,params,stic);
		break;
	case DOCK:
		return RunDock(job,params,stic);
	default:
                printf("(%s): ERROR - Unknown job type: %d\n",params->node_tag,job->type);
		break;
	}
	return retval;
}


/**
 * Performs functions on the job which should happen after
 * that particular type of job finishes. The job node is then removed
 * from the job Tree. Any state stored in finished_job using calloc is
 * freed. HOWEVER, this function does not free finished_job itself.
 * It may be used, and eventually should be freed, after this function
 * is called.
 *
 * Returns finished_job
 */
deque_node* finish_job(deque *jobs, deque *busy_list,
		deque_node *finished_job,
		       const char *output, int return_val, JobParameters *params)
{

	job_t *job_info;
	int job_type;
	double G_prot = DOUBLE_INIT;

	if(finished_job == NULL || finished_job->type != DEQUE_JOB)
		return NULL;

	job_info = (job_t*)finished_job->data;
	job_type = job_info->type;

	if(busy_list != NULL)
		remove_node(busy_list,finished_job);

	if(job_info->type == PRESCREENING)
	{
		job_info->host_rank = -1;
		job_info->type = DOCK;
		if(return_val == 0)
		  push_back(jobs,finished_job);
		else
		  {
		    FreeJob(job_info);
		    free(finished_job);finished_job = NULL;
		  }
	}
	else if(job_info->type == DOCK)
	{
		if(params->doMOPAC && return_val == 0)
		  {
		    int repcounter = 0;
		    struct ClusterRep *curr_rep = NULL;
		    job_info->host_rank = -1;
		    job_info->type = QM;
		    curr_rep = job_info->input_data.reps;
		    while(curr_rep != NULL)
		    {
		    	int struct_count = 0;
		    	repcounter++;
		    	for (struct_count = 0; struct_count < 2; struct_count++) {
		    		deque_node *new_qm_job = NULL;
		    		job_t *new_job = NULL;
		    		new_qm_job = InitDequeNode(new_qm_job);
		    		new_qm_job->type = DEQUE_JOB;
		    		new_qm_job->data = InitJob((job_t*)new_qm_job->data);
		    		new_job = (job_t*)new_qm_job->data;
		    		//new_job->name = (char*)calloc(FILENAME_MAX,sizeof(char));
		    		if (struct_count == 0) {
		    			sprintf(new_job->name,"%s_%03d",job_info->name,repcounter);
		    			new_job->type = QMLIG;
		    		}
		    		else {
		    			sprintf(new_job->name,"%s_%s_%03d",params->receptor_name,job_info->name,repcounter);
		    			new_job->type = QM;
		    		}
		    		new_job->host_rank = -1;
		    		new_job->input_data = job_info->input_data;
		    		new_job->input_data.reps = NULL;
		    		push_back(jobs,new_qm_job);
		    	}
		    	curr_rep = curr_rep->next;
		    }
		  }
		FreeJob(job_info);
		free(finished_job);finished_job = NULL;
	}
	else
	{
            if(return_val != 0 || strlen(job_info->mopac_res) == 0) {
                FreeJob(job_info);
                free(finished_job);finished_job = NULL;
            }
            else
	        push_back(jobs,finished_job);
	}

	return finished_job;
}

int SendMapFiles(MPI_Comm transfer_comm, int rank, int root_group_rank, JobParameters *params, const char *transfer_dir)
{
	DIR *dirhandle;
	struct dirent *file;
	char tmp_name[FILENAME_MAX],d_name[FILENAME_MAX];
	FILE *map_file = NULL;
	char *massive_buffer;
        char filestub[FILENAME_MAX];
	int num_read = 0;
	enum {loop_terminate = 0,read_write,new_file};
	int continue_signal = new_file;//2 = new file, 1 = keep reading/writing. 0 terminate.
	int is_master = (rank == MASTER);
	int filename_size;
	time_t start_time, stop_time;
	long int bytes_transfered = 0;

    #ifdef QUICK_DEBUG
	    return 0;
    #endif

    #ifdef DEBUG
	    printf("DEBUG - %s: Entered SendMapFiles (rank = %d)\n",params->node_tag,rank);
    #endif 
	// Let's try calloc here
	// massive_buffer = (char*)malloc(FILE_BUFF_SIZE);
	massive_buffer = (char*)calloc(1,FILE_BUFF_SIZE);

	start_time = stop_time = time(NULL);

        sprintf(filestub,"%s.",params->receptor_name);

	if(is_master)
	{
		continue_signal = new_file;
		dirhandle = opendir(transfer_dir);
		if(dirhandle == NULL)
		  {
		    fprintf(stderr,"Master could not open the receptor directory: %s\n",transfer_dir);
		    fprintf(stderr,"Reason: %s\n",strerror(errno));
		    exit(FILESYS_ERROR);
		  }
		start_time = time(NULL);
	}
	else
	  {
	    if (VerifyDir(transfer_dir,1,params->node_tag) != 0)
	      {
		fprintf(stderr,"%s: Could not create directory needed for receptor transfer %s\n",params->node_tag,transfer_dir);
		exit(FILESYS_ERROR);
	      }
	  }
	MPI_Barrier(transfer_comm);//all start the loop together

	while(continue_signal != loop_terminate)//loop terminates when readdir has no more files.
	{
		if(is_master)
		{
			struct stat file_status;
			file=readdir(dirhandle);
			if(file == NULL)
			{
				continue_signal = loop_terminate;
			}
			else
			{
				// if(strlen(file->d_name) == 0 || file->d_name[0] == '.')
				if(strncmp(filestub,file->d_name,strlen(filestub)) != 0)
					continue_signal = new_file;
				else
				{
					sprintf(tmp_name,"%s/%s",transfer_dir,file->d_name);
					if(stat(tmp_name,&file_status) || S_ISDIR(file_status.st_mode))
						continue_signal = new_file;
					else
					{
						strcpy(d_name,file->d_name);
						filename_size = strlen(d_name) + 1;// + 1 for the null character
						continue_signal = read_write;
						bytes_transfered += file_status.st_size;
						printf_verbose(params->verbose,"Sending to transfer nodes: %s\n",d_name);
					}
				}
			}
		}
		MPI_Bcast(&continue_signal,1,MPI_INT,root_group_rank,transfer_comm);
		if(continue_signal == loop_terminate)
			break;
		else if(continue_signal == new_file)
			continue;

		MPI_Bcast(&filename_size,1,MPI_INT,root_group_rank,transfer_comm);
		MPI_Bcast(&d_name,filename_size,MPI_CHAR,root_group_rank,transfer_comm);
		if(!is_master)
			sprintf(tmp_name,"%s/%s",transfer_dir,d_name);
		if(map_file != NULL)
			fclose(map_file);

		if(is_master)
		{
			map_file = fopen(tmp_name,"r");
		}
		else
			map_file = fopen(tmp_name,"w");

		if(map_file == NULL)
		  {
		    fprintf(stderr,"%s: (rank %d) could not open file: %s\n",params->node_tag,rank,tmp_name);
		    fprintf(stderr,"%s: Reason: %s\n",params->node_tag,strerror(errno));
		    exit(FILESYS_ERROR);
		  }

		MPI_Bcast(&continue_signal,1,MPI_INT,root_group_rank,transfer_comm);

		while(continue_signal == read_write)//Loop terminates when master (via getmapsegment) says there is nothing more to read/write (i.e. num_read = 0)
		{
			if(is_master)
			{
			  num_read = fread(massive_buffer,1,FILE_BUFF_SIZE,map_file);
				if(num_read == 0)
					continue_signal = new_file;
			}

			MPI_Bcast(&continue_signal,1,MPI_INT,root_group_rank,transfer_comm);
			if(continue_signal != read_write)
				break;
			MPI_Bcast(&num_read,1,MPI_INT,root_group_rank,transfer_comm);
			MPI_Bcast(massive_buffer,num_read,MPI_CHAR,root_group_rank,transfer_comm);
			if(!is_master)
			{
				fwrite(massive_buffer,1,num_read,map_file);
			}
			else
				if(feof(map_file))
					continue_signal = new_file;

			MPI_Bcast(&continue_signal,1,MPI_INT,root_group_rank,transfer_comm);

		}
	}

	if(map_file != NULL)
		fclose(map_file);

	if(is_master)
	  {
	    stop_time = time(NULL);
	    // printf_verbose(params->verbose,"Total transfered: %fMb\n",bytes_transfered/1e6);
	    // printf_verbose(params->verbose,"Total time in seconds for transfer: %f\n",difftime(stop_time,start_time));
	    printf_verbose(params->verbose,"Transferred %.2f Mb in %.2f seconds",bytes_transfered/1e6,difftime(stop_time,start_time));
	    printf_verbose(params->verbose," (Average transfer rate: %.2f Mb/s)\n",bytes_transfered/(1e6*difftime(stop_time,start_time)));
	  }

	free(massive_buffer);
	return 0;
}

void FPackBufferSTICTree(FILE *stic_file,RBTree *stics,const char *cmpd_name)
{
	void *buffer;
	size_t size;
	if(stics == NULL)
		return;

	FPackBufferSTICTree(stic_file,stics->left,cmpd_name);

	PackBufferSTIC(&buffer,&size,cmpd_name,0,(struct STICelement*)stics->data);
	fwrite(buffer,1,size,stic_file);
	fflush(stic_file);
	free(buffer);

	FPackBufferSTICTree(stic_file,stics->right,cmpd_name);
}

void FSaveCompoundList(FILE *compound_file, FILE *stic_file, RBTree *CompoundList)
{
	tpl_node *tn;
	struct compoundData compound;
	if(CompoundList == NULL)
		return;
	FSaveCompoundList(compound_file, stic_file, CompoundList->left);

	compound = *((struct compoundData*)((CompoundTree*)CompoundList->data)->data);
	tn = tpl_map(COMPOUND_TPL_MAP,&compound);
	tpl_pack(tn,1);
	tpl_dump(tn,TPL_FD,fileno(compound_file));
	tpl_free(tn);

	FPackBufferSTICTree(stic_file,((CompoundTree*)CompoundList->data)->stics,compound.ID);

	FSaveCompoundList(compound_file,stic_file,CompoundList->right);
}

void SaveCompoundTree(RBTree *CompoundList)
{
	FILE *compound_file, *stic_file;
	compound_file = fopen(COMPOUND_RESTART_FILE,"w");
	stic_file = fopen(STIC_RESTART_FILE,"w");

	FSaveCompoundList(compound_file,stic_file,CompoundList);
	fclose(compound_file);
	fclose(stic_file);
}

void SaveState(deque *jobs, deque *busy_list, RBTree *CompoundList, JobParameters *params)
{
	job_t job;
	deque_node *curr_node;
	tpl_node *tn;

	if(!params->restart_job)
		return;

	tn = tpl_map(JOBLIST_TPL_MAP,&job);
	curr_node = jobs->head;
	while(curr_node != NULL)
	{
		job = *((job_t*)curr_node->data);
		tpl_pack(tn,1);
		curr_node = curr_node->right;
	}
	tpl_dump(tn,TPL_FILE,JOBS_RESTART_FILE);
	tpl_free(tn);

	tn = tpl_map(JOBLIST_TPL_MAP,&job);
	curr_node = busy_list->head;
	while(curr_node != NULL)
	{
		job = *((job_t*)curr_node->data);
		tpl_pack(tn,1);
		curr_node = curr_node->right;
	}
	tpl_dump(tn,TPL_FILE,BUSY_RESTART_FILE);
	tpl_free(tn);

	SaveCompoundTree(CompoundList);
}

int RestoreCompoundTree(RBTree **CompoundList)
{
	FILE *compound_tpl, *stic_tpl;
	int stic_idx = 0;
	tpl_node *tn_cmpd;
	struct compoundData new_compound;
	struct STICelement *new_stic = NULL;

	if(CompoundList == NULL)
		return TRUE;

	compound_tpl = fopen(COMPOUND_RESTART_FILE,"r");
	stic_tpl = fopen(STIC_RESTART_FILE,"r");

	tn_cmpd = tpl_map(COMPOUND_TPL_MAP,&new_compound);

	while(tn_cmpd != NULL && tpl_load(tn_cmpd,TPL_FD,fileno(compound_tpl)) == 0)
	{
		if(!tpl_unpack(tn_cmpd,1))
			continue;
		RBTree *compound_node = NULL;
		// CompoundTree *new_cmpd_tree = (CompoundTree*)malloc(sizeof(CompoundTree));
		CompoundTree *new_cmpd_tree = (CompoundTree*)calloc(1,sizeof(CompoundTree));
		new_cmpd_tree->data = InitCompound(new_cmpd_tree->data);
		new_cmpd_tree->stics = NULL;
		stic_idx = 0;
		for(;stic_idx < new_compound.num_STIC;stic_idx++)
		{
			int retval;
			RBTree *new_stic_node = NULL;
			new_stic = InitSTIC(new_stic);
			retval = UnpackFileSTIC(stic_tpl,new_compound.ID, 0, new_stic);
			if(retval)
			{
				fprintf(stderr,"RestoreCompoundTree: Could not unpack stic. Expected: %d. Unpaced: %d\n",new_compound.num_STIC,stic_idx);
				exit(-1);
			}

			new_stic_node = InitRBTree(new_stic_node,STIC);
			new_stic_node->data = new_stic;
			insertRBT(&new_cmpd_tree->stics,new_stic_node);
		}
		*new_cmpd_tree->data = new_compound;
		compound_node = InitRBTree(compound_node,COMPOUND);
		compound_node->data = new_cmpd_tree;
		insertRBT(CompoundList,compound_node);
	}
	return TRUE;
}

int RestoreState(cpunode **free_cpus,
		deque *jobs, deque *busy_list,
		RBTree **CompoundList,
		JobParameters *params)
{
	tpl_node *tn;
	job_t loadedjob;
	int retval;

	if(access(JOBS_RESTART_FILE,R_OK) != 0)
		return FALSE;
	if(access(BUSY_RESTART_FILE,R_OK) != 0)
			return FALSE;
	if(access(COMPOUND_RESTART_FILE,R_OK) != 0)
			return FALSE;
	if(access(STIC_RESTART_FILE,R_OK) != 0)
			return FALSE;

	//load save jobs list
	clear_deque(jobs);
	clear_deque(busy_list);
	tn = tpl_map(JOBLIST_TPL_MAP,&loadedjob);
	tpl_load(tn,TPL_FILE,JOBS_RESTART_FILE);
	while(tpl_unpack(tn,1))
	{
		deque_node *new_job_node = NULL;
		job_t *newjob = NULL;
		newjob = InitJob(newjob);
		new_job_node = InitDequeNode(new_job_node);
		new_job_node->type = DEQUE_JOB;
		new_job_node->data = newjob;
		*newjob = loadedjob;
		push_back(jobs,new_job_node);
	}
	tpl_free(tn);

	//move busy list to the front of the job stack
	tn = tpl_map(JOBLIST_TPL_MAP,&loadedjob);
	tpl_load(tn,TPL_FILE,BUSY_RESTART_FILE);
	while(tpl_unpack(tn,1))
	{
		deque_node *new_job_node = NULL;
		job_t *newjob = NULL;
		newjob = InitJob(newjob);

		new_job_node = InitDequeNode(new_job_node);
		new_job_node->type = DEQUE_JOB;
		new_job_node->data = newjob;
		*newjob = loadedjob;
		push_front(jobs,new_job_node);
	}
	tpl_free(tn);

	retval = RestoreCompoundTree(CompoundList);
	if(!retval)
		printf("WARNING - RestoreState could not restore the compound tree.");

    #ifdef DEBUG
	    printf_master("DEBUG - Compound Tree after restart:\n");
	    FPrintCompoundTree(stdout, *CompoundList);
    #endif

	return retval;

}

void SetupEffercio( cpunode **free_cpus,
		deque *jobs, deque *busy_list,
		RBTree **CompoundList,
		int rank, int node_count,
		JobParameters *params)
{
	int destination = 0;
	int restored = FALSE;

	if(params->restart_job)
		restored = RestoreState(free_cpus,jobs, busy_list,CompoundList,params);


	// Kick off the initial round of jobs, one to each slave
	// Send out as many jobs as there are slave nodes, unless
	// there are fewer jobs than nodes. In which case, some nodes
	// get nothing (hopefully later they will) and are pushed
	// into a wait queue (free_cpus).
	for (; destination < node_count; destination++) {
		//Loop iterates through all of the slaves.
		//If we run out of jobs to do, leave this distribution loop
		if(rank == destination)
			continue;
		if(jobs->head == NULL)
			suspend_cpu(free_cpus,destination);
		else
			issue_order(jobs,busy_list,destination,rank,params);
	}
}

void DefaultParameters(JobParameters *params)
{
	if(params == NULL)
		return;
	params->dock_params = calloc(1,sizeof(char *));
	params->mopac_footer_params = calloc(1,sizeof(char*));
	params->mopac_header_params = calloc(1,sizeof(char*));
	params->num_dock_params = 0;
	params->num_mopac_footer_params = 0;
	params->num_mopac_header_params = 0;
	memset(params->receptor_dir,0,FILENAME_MAX);
	memset(params->receptor_name,0,FILENAME_MAX);
	memset(params->ligand_dir,0,FILENAME_MAX);
	memset(params->results_dir,0,FILENAME_MAX);
	memset(params->optimized_dir,0,FILENAME_MAX);
	memset(params->clusters_dir,0,FILENAME_MAX);
	memset(params->analysis_dir,0,FILENAME_MAX);
	memset(params->username,0,LOGIN_NAME_MAX);
	memset(params->node_tag,0,HOST_NAME_MAX);
	params->qm_method = "PM6"; // can use PM6-D3H4 as well
	params->doMOPAC = FALSE;
	params->restart_job = FALSE;
	params->prescreen = FALSE;
	params->input_units = INPUT_KCAL;
	params->start_time = time(NULL);
	params->total_time = DAY_IN_SEC * 365;
	params->target_G = params->target_Hf = DOUBLE_INIT;
}

int main(int argc, char **argv)
{
	//int rank, i, size, nints, len, num_jobs;
	int rank, i, size, num_jobs;

	char infilename[FILENAME_MAX]; // file containing the compound names
	char restart_filename[FILENAME_MAX];
	//FILE *jobfile;
	char hostname[HOST_NAME_MAX];
	char sys_command[ARG_MAX];
	char nextjob[FILENAME_MAX];
	char home_dir[FILENAME_MAX];

	deque *jobs = NULL;//Queue of jobs to be run
	deque *busy_list = NULL;//List of jobs currently running

	//For efficiency's sake, keep track of next job to be run.
	//By using this and issue_order(...), Effercio will treat
	//the tree like it is a queue for faster job spawning.
	cpunode *free_cpus = NULL;//List of CPUs sitting idle
	int master_node = FALSE;
	size_t BUFFER_SIZE = (size_t) (MSG_BUFSIZ * MEM_BLOCK_SIZE);

	char *timestring;
	RBTree *CompoundList = NULL;

	// MPI data
	MPI_Group transfer_group;
	MPI_Comm transfer_comm;


	JobParameters params;

	memset(infilename,0,FILENAME_MAX);
	memset(restart_filename,0,FILENAME_MAX);
	memset(hostname,0,HOST_NAME_MAX);
	memset(sys_command,0,ARG_MAX);
	memset(nextjob,0,FILENAME_MAX);
	memset(home_dir,0,FILENAME_MAX);


	// Initialize Parameters that are passed along to jobs
	DefaultParameters(&params);




	const char *opts_short = ":hvRGl:m:d:c:a:qo:p:s:t:k:";
	enum {FLAGS_MOPAC_HEADER = '{',FLAGS_MOPAC_FOOTER,FLAGS_SCRATCH_DIR='w'};
	const struct option opts_long[] = {
			{"help"        ,0,NULL,'h'}, // print usage information
			{"verbose"     ,0,NULL,'v'}, // turn on verbose mode
			{"ligands"     ,1,NULL,'l'}, // directory containing the ligand .pdbqt files
			{"maps"        ,1,NULL,'m'}, // directory containing the receptor .pdbqt and .map files
			{"docks"       ,1,NULL,'d'}, // directory to store the docking results
			{"cluster-reps",1,NULL,'c'}, // directory to store the cluster representatives
			{"analysis"    ,1,NULL,'a'}, // directory to store analysis reports
			{"screen"      ,0,NULL,'s'}, // flag to prescreen the ligands with MOPAC
			{"quantum"     ,0,NULL,'q'}, // flag to turn on MOPAC refinement
			{"optimized"   ,1,NULL,'o'}, // directory to store the MOPAC geometry optimizations
			{"param"       ,1,NULL,'p'}, // parameter string for AutoDock
			{"time"        ,1,NULL,'t'}, // total time for job in seconds
			{"Gibbs"       ,0,NULL,'G'}, // flag to turn on use of free energies rather than enthalpies
			{"restart"     ,0,NULL,'R'}, // restart from an incomplete previous job
			{"mopac"       ,1,NULL,'k'}, // Keyword(s) to be added to the mopac header *and* footer
			{"mopac-header",1,NULL,FLAGS_MOPAC_HEADER}, // Keyword(s) to be added to the mopac header
			{"mopac-footer",1,NULL,FLAGS_MOPAC_FOOTER}, // Keyword(s) to be added to the mopac footer
			{"scratch-dir" ,1,NULL,FLAGS_SCRATCH_DIR},// Scratch Directory
			{NULL, 0, NULL,0}
	};

	char opt_switch = (char) 0;

	MPI_Status status;

    #ifdef DEBUG
	    srand((unsigned)(time(0)));
    #endif


	MPI_Init(&argc,&argv);
	printf("Finished MPI_Init\n");
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (gethostname(hostname,HOST_NAME_MAX) == -1) {
		printf("ERROR - Failed to get hostname\n");
		printf("ERROR - %s\n",strerror(errno));
		exit(EXIT_FAILURE);
	}
	if (getcwd(home_dir,FILENAME_MAX) == NULL) {
		printf("ERROR - Failed to get current directory\n");
		printf("ERROR - %s\n",strerror(errno));
		exit(EXIT_FAILURE);
	}

	// Check the length of the strings and truncate if necessary
	if ((strlen(hostname) + 15) <= HOST_NAME_MAX) {
		sprintf(params.node_tag,"Process %4d (%s)",rank,hostname);
	}
	else {
		int templength;
		templength = strlen(hostname) + 15 - HOST_NAME_MAX;
		hostname[templength] = 0;
        if (params.verbose) {
		    sprintf(params.node_tag,"Process %4d running on %s with pid of %d",rank,hostname,getpid());
        }
        else {
		    sprintf(params.node_tag,"Process %4d (%s)",rank,hostname);
        }
	}

	if (rank == MASTER) master_node = TRUE;

	// Send start_time to everyone
	MPI_Bcast(&params.start_time,sizeof(time_t),MPI_BYTE,MASTER,MPI_COMM_WORLD);

///////////////////////////////////////////////////////////////////////////////
//
// This is a "trap" for debugging with DDD/gdb.  This block can be commented
// out if you simply wish to debug via the debugging output
// 
#ifdef GDBDEBUG
	if (master_node) {
		int WaitDebug = 1;
		printf("MASTER PROCESS IS WAITING: pid = %d\n",getpid());
		while (WaitDebug) ;
	}
#endif
//
///////////////////////////////////////////////////////////////////////////////
	if (master_node) {
		PrintBanner();
		printf("%s started on: %s\n",PROG_NAME,asctime(localtime(&params.start_time)));
		fflush(stdout);
	}

	// Send username to all nodes
	if(master_node)
		strcpy(params.username,getenv("USER"));
	MPI_Bcast(&params.username,LOGIN_NAME_MAX,MPI_CHAR,MASTER,MPI_COMM_WORLD);

	// Set the default locations for output
	params.useDefaultScratch = 1;
	strcpy(params.scratch_dir,"/scratch_dir");
	strcpy(params.receptor_dir ,home_dir);
	strcpy(params.ligand_dir   ,home_dir);
	strcpy(params.results_dir  ,home_dir);
	strcpy(params.clusters_dir ,home_dir);
	strcpy(params.optimized_dir,home_dir);
	strcpy(params.analysis_dir ,home_dir);

	strcat(params.receptor_dir ,"/receptors");
	strcat(params.ligand_dir   ,"/ligands");
	strcat(params.results_dir  ,"/dockings");
	strcat(params.clusters_dir ,"/optprep");
	strcat(params.optimized_dir,"/optimized");
	strcat(params.analysis_dir ,"/analysis");

	strcat(restart_filename,PROG_NAME);
	strcat(restart_filename,".restart");
	strcat(params.mgl_bin_dir,MGL_BIN_DIR);
	strcat(params.autodock_exe,AUTODOCK_EXE);
	sprintf(params.mopac_path,"%s%s",MOPAC_HOME,MOPAC_EXE);

	// Parse the command line and do some mild error checking
	int opt_count = 0; // number of options processed
	int curr_flag = 0;
	while((opt_switch = getopt_long(argc,argv,opts_short,opts_long,&curr_flag)) != -1){
	  if (optarg == NULL) opt_count = opt_count + 1; 
	  else opt_count = opt_count + 2; 

		switch(opt_switch) {
		case 'h':
            #ifdef DEBUG
			printf_master("DEBUG - Found -h argument\n");
            #endif
			PrintHelp();
			exit(0);
		case 'v':
            #ifdef DEBUG
			printf_master("DEBUG - Found -v argument\n");
            #endif
			params.verbose = TRUE;
			printf_master("WARNING: ********************************\n");
			printf_master("WARNING: *** Running in verbose mode! ***\n");
			printf_master("WARNING: *** Expect LOTS of output    ***\n");
			printf_master("WARNING: ********************************\n\n");
			break;
		case 'l':
			strcpy(params.ligand_dir,optarg);
            #ifdef DEBUG
			printf_master("DEBUG - Found -l argument\n");
			printf_master("DEBUG - ligand_dir = %s\n",params.ligand_dir);
            #endif
			break;
		case 'm':
			strcpy(params.receptor_dir,optarg);
            #ifdef DEBUG
			printf_master("DEBUG - Found -m argument\n");
			printf_master("DEBUG - receptor_dir = %s\n",params.receptor_dir);
            #endif
			break;
		case 'd':
			strcpy(params.results_dir,optarg);
            #ifdef DEBUG
			printf_master("DEBUG - Found -d argument\n");
			printf_master("DEBUG - results_dir = %s\n",params.results_dir);
            #endif
			break;
		case 'c':
			strcpy(params.clusters_dir,optarg);
            #ifdef DEBUG
			printf_master("DEBUG - Found -c argument\n");
			printf_master("DEBUG - clusters_dir = %s\n",params.clusters_dir);
            #endif
			break;
		case 'a':
			strcpy(params.analysis_dir,optarg);
            #ifdef DEBUG
			printf_master("DEBUG - Found -a argument\n");
			printf_master("DEBUG - analysis_dir = %s\n",params.analysis_dir);
            #endif
			break;
		case 'q':
            #ifdef DEBUG
			printf_master("DEBUG - Found -q argument\n");
            #endif
			params.doMOPAC = TRUE;
			printf_master("WARNING: ********************************************************\n");
			printf_master("WARNING: *** Performing two-stage docking with semi-empirical ***\n");
			printf_master("WARNING: *** refinement of the docked structures.  This WILL  ***\n");
			printf_master("WARNING: *** require a GREAT DEAL of computational time.  Be  ***\n");
			printf_master("WARNING: *** certain you know what you're doing!!  Such runs  ***\n");
			printf_master("WARNING: *** typically require about 1 CPU-day per ligand.    ***\n");
			printf_master("WARNING: ********************************************************\n");
			printf_master("\n");
			break;
		case 'o':
			strcpy(params.optimized_dir,optarg);
            #ifdef DEBUG
			printf_master("DEBUG - Found -o argument\n");
			printf_master("DEBUG - optimized_dir = %s\n",params.optimized_dir);
            #endif
			break;
		case 'p':
			add_strtoarray(&params.dock_params,&params.num_dock_params,optarg);
            #ifdef DEBUG
			printf_master("DEBUG - Found -p argument\n");
			printf_master("DEBUG - docking parameters:\n");
			for (i=0;i<params.num_dock_params;i++)
				printf_master("DEBUG - %s\n",params.dock_params[i]);
            #endif
			break;
		case 't':
			// Need to work on this to allow Minutes, Hours, or Days
			// Right now we deal only in seconds
			params.total_time = (time_t) strtod(optarg,NULL);
			break;
		case 's':
            #ifdef DEBUG
           printf_master("DEBUG - Found -s argument\n");
           #endif
           params.prescreen = TRUE;
           printf_master("Input geometries will be prescreened using semiempirical quanutm mechanics.\n");
           break;
			break;
		case 'G':
            #ifdef DEBUG
			printf_master("DEBUG - Found -G argument\n");
            #endif
			params.UseFreeEnergy = TRUE;
			printf_master("WARNING: ***********************************************************\n");
			printf_master("WARNING: *** Maxwell-Boltzmann averaging of the docking results  ***\n");
			printf_master("WARNING: *** will use calculated Gibbs free energies.  This WILL ***\n");
			printf_master("WARNING: *** require a GREAT DEAL of computational time.  Be     ***\n");
			printf_master("WARNING: *** certain you know what you're doing!! Computing the  ***\n");
			printf_master("WARNING: *** free energy of a protein-ligand complex typically   ***\n");
			printf_master("WARNING: *** requires about 1 CPU-day per ligand.                ***\n");
			printf_master("WARNING: ***********************************************************\n");
			printf_master("\n");
			break;
		case 'R':
			params.restart_job = TRUE;
            #ifdef DEBUG
			printf_master("DEBUG - Found -R argument\n");
            #endif
			printf_master("This is a restart of a previous run.  Using '%s' as the restart file.\n",
					restart_filename);
			break;
		case 'k':
            #ifdef DEBUG
            printf_master("DEBUG - Adding the following to mopac headers and footers:\n%s\n",optarg);
            #endif
			add_strtoarray(&params.mopac_header_params,&params.num_mopac_header_params,optarg);
			add_strtoarray(&params.mopac_footer_params,&params.num_mopac_footer_params,optarg);
			break;
		case FLAGS_MOPAC_HEADER:
            #ifdef DEBUG
            printf_master("DEBUG - Adding the following to mopac headers:\n%s\n",optarg);
            #endif
			add_strtoarray(&params.mopac_header_params,&params.num_mopac_header_params,optarg);
			break;
		case FLAGS_MOPAC_FOOTER:
            #ifdef DEBUG
            printf_master("DEBUG - Adding the following to mopac footers:\n%s\n",optarg);
            #endif

			add_strtoarray(&params.mopac_footer_params,&params.num_mopac_footer_params,optarg);
			break;
		// case FLAGS_SCRATCH_DIR:
		case 'w':
            #ifdef DEBUG
            printf_master("DEBUG - Setting scratch directory to %s\n\n",optarg);
            #endif
            strcpy(params.scratch_dir,optarg);
            params.useDefaultScratch = 0;
            break;
        #ifdef DEBUG
		case '?':
            printf_master("DEBUG - Unrecognized option: ",argv);
            break;
        #endif
		default:
			break;
		}
	}

	if (argc - optind < 2) {
        printf_master("ERROR: %s requires at least two arguments!!! There were only %d, %d of which have been parsed\n\n",
                      PROG_NAME, argc, optind);
		PrintHelp();
		MPI_Abort(MPI_COMM_WORLD,-1);
	}
	strcpy(params.receptor_name,argv[optind++]);
	strcpy(infilename   ,argv[optind]);

    // Override certain switches for debugging purposes
    #ifdef DEBUG
    // DEBUG implies verbose output
	params.verbose = TRUE;
    //
    // The following two flags will turn on the SEQM parts to check
    // if they are working properly
	// params.prescreen = TRUE;
	// params.doMOPAC = TRUE;
    #endif

    // If we're using MOPAC, then we need to prescreen to get the ligand energies
    // for the Boltzmann averaging
    if (params.doMOPAC) params.prescreen = TRUE;

	// Print some useful info about the run
	printf_master("Total time alloted    : %d seconds\n",params.total_time);
	printf_master("Using receptor name   : %s\n",params.receptor_name);
	printf_master("Ligand names from     : %s\n",infilename);
	printf_master("Receptor files (dir)  : %s\n",params.receptor_dir);
	printf_master("Ligand files (dir)    : %s\n",params.ligand_dir);
	printf_master("Docking results (dir) : %s\n",params.results_dir);
	printf_master("Cluster structs (dir) : %s\n",params.clusters_dir);
	printf_master("MOPAC results (dir)   : %s\n",params.optimized_dir);
	printf_master("Analysis results (dir): %s\n",params.analysis_dir);
	printf_master("MGLTools binary (dir) : %s\n",params.mgl_bin_dir);
	printf_master("Autodock executable   : %s\n",params.autodock_exe);
	printf_master("MOPAC executable      : %s\n",params.mopac_path);
	printf_master("Restart file          : %s\n",restart_filename);

#ifdef DEBUG
	// Free all the strings in dock_params and repopulate it with some parameters
	// more suitable for debugging purposes
	FreeStringArray(&params.dock_params,&params.num_dock_params);
	params.dock_params = realloc(params.dock_params,sizeof(char**));
	add_strtoarray(&params.dock_params,&params.num_dock_params,"ga_run=2");
	add_strtoarray(&params.dock_params,&params.num_dock_params,"ga_num_evals=10000");
	printf_master("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	printf_master("!! Ignoring user-supplied docking parameters for DEBUG run !!\n");
	printf_master("!!!!!!!!!! Using the following parameters instead !!!!!!!!!!!\n");
	printf_master("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
#endif

	for (i=0;i<params.num_dock_params;i++) {
		printf_master("\nAutoDock Parameters\n");
		printf_master("%s\n",DASH_HDR);
		printf_master("    %s\n",params.dock_params[i]);
	}
	printf_master("\n");

	printf_master("Running with %d processors\n\n",size);

	// Check access to executables
	if(check_extern_apps(&params) != 0)
	  {
	    printf("%s ERROR - Cannot access external executables.\n",params.node_tag);
	    if(errno)
	      MPI_Abort(MPI_COMM_WORLD,errno);
	    else
	      MPI_Abort(MPI_COMM_WORLD,1);
	  }
	if(check_directories(&params) != 0)
	  {
	    printf("%s ERROR - Cannot access input/output directories.\n",params.node_tag);
	    if(errno)
	      MPI_Abort(MPI_COMM_WORLD,errno);
	    else
	      MPI_Abort(MPI_COMM_WORLD,1);
	  }
	
	//
	// MAIN LOOP
	//
	int *unique_ranks = calloc(1,sizeof(int));
	int num_unique = 0;
	//int unique_count = num_unique;
	if (master_node) {
		char **unique_hosts = calloc(1,sizeof(char **));
		char tempname[HOST_NAME_MAX],outstring[OUTPUT_WIDTH];
		int outwidth = 0;

		memset(tempname,(char) 0,HOST_NAME_MAX);
		memset(outstring,(char) 0,OUTPUT_WIDTH);


		fflush(NULL);

		if (params.verbose) {
		  printf("Master process is running on %s (pid: %lu)\n",hostname,getpid());
			printf("Slave processes are running on:\n");
			for (i=1;i<OUTPUT_WIDTH;i++) printf("-"); printf("\n");
		}
		fflush(stdout);
		outwidth = 0;
		for (i=1; i < size; i++) {
			int is_unique = TRUE;

			MPI_Recv (tempname,HOST_NAME_MAX, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG,
					MPI_COMM_WORLD,&status);
			if (params.verbose) {
				if ((outwidth + strlen(tempname)) <= OUTPUT_WIDTH) {
					printf("%s",tempname);
					outwidth = outwidth + strlen(tempname);
					if ((outwidth + 1) <= OUTPUT_WIDTH) {
						printf(" ");
						fflush(stdout);
						outwidth = outwidth + 1;
					}
					else {
						printf("\n");
						fflush(stdout);
						outwidth = 0;
					}
				}
				else {
					printf("\n%s ",tempname);
					fflush(stdout);
					outwidth = strlen(tempname);
				}
			}
			is_unique = TRUE;
			if (i > 1) {
				// Here's where we check for unique hostnames
				int j = 0;
				while (j < num_unique) {
					if (strcmp(unique_hosts[j],tempname) == 0) {
						is_unique = FALSE;
						break;
					}
					j++;
				}
			}
			if (is_unique) {
				add_strtoarray(&unique_hosts,&num_unique,tempname);
				// Decrement num_unique before sending to add_inttoarray as
				// the function will increment it again
				num_unique--;
				add_inttoarray(&unique_ranks,&num_unique,status.MPI_SOURCE);
			}
		}
		if (params.verbose) {
			printf("\n");
			for (i=1;i<OUTPUT_WIDTH;i++) printf("-"); printf("\n");
			printf("Unique Hosts are:\n");
			for (i=0;i<num_unique;i++) {
				printf("%s (rank = %d)\n",unique_hosts[i],unique_ranks[i]);
			}
			printf("\n");
		}


		fflush(stdout);
		// Free unique_hosts
		int temp_num_unique = num_unique;
		FreeStringArray(&unique_hosts,&temp_num_unique);
		add_inttoarray(&unique_ranks,&num_unique,MASTER);

	}//END MASTER SECTION
	else {
		MPI_Send(&hostname,HOST_NAME_MAX,MPI_CHAR, MASTER, rank, MPI_COMM_WORLD);
	}

	// Use only the unique nodes for transferring files
	// Send num_unique to everyone
	MPI_Bcast(&num_unique,1,MPI_INT,MASTER,MPI_COMM_WORLD);
	if(!master_node)
		unique_ranks = realloc(unique_ranks,num_unique*sizeof(int));

	MPI_Bcast(unique_ranks,num_unique,MPI_INT,MASTER,MPI_COMM_WORLD);
	//MPI_Group everyone;
	//MPI_Comm_group(MPI_COMM_WORLD, &everyone);
	//MPI_Group_incl(everyone,num_unique,unique_ranks,&transfer_group);
	//MPI_Comm_create(MPI_COMM_WORLD,transfer_group,&transfer_comm);

	for (i=0;i<num_unique;i++) {
	  int j;
	  for (j=0;j<size;j++) {
            if (j == unique_ranks[i]) {
	      if (rank == j) params.transfer_node = TRUE;
            }
	  }
	}
	MPI_Comm_split(MPI_COMM_WORLD,params.transfer_node,rank,&transfer_comm);
	int root_group_rank;
	if(master_node)
	  MPI_Comm_rank(transfer_comm,&root_group_rank);
	MPI_Bcast(&root_group_rank,1,MPI_INT,MASTER,MPI_COMM_WORLD);

	//Send Receptor map files to all of the transfer nodes (those in the transfer_group)
	if(master_node)
	{
	  printf("Master is calling SendMapFiles (%d, %d)\n",rank,root_group_rank);
	  SendMapFiles(transfer_comm, rank, root_group_rank,&params,params.receptor_dir);
	}
	else
	{
		char workdir[FILENAME_MAX];
		int group_rank;
		strcpy(workdir,params.scratch_dir);
		if (params.useDefaultScratch)
		{
			strcat(workdir,"/");
			strcat(workdir,params.username);
		}
		printf("%s: DEBUG - scratch directory set to %s\n",params.node_tag,workdir);
		if(params.transfer_node)
		  {
		    printf("%s is calling SendMapFiles (%d, %d)\n",params.node_tag,rank,root_group_rank);
		    SendMapFiles(transfer_comm, rank, root_group_rank, &params,workdir);
		  }
	}
	MPI_Barrier(MPI_COMM_WORLD);




	// Free unique_ranks
	free(unique_ranks);
	unique_ranks = NULL;

	// Send qm_method to everyone
	MPI_Bcast(&params.qm_method,strlen(params.qm_method)+1,MPI_CHAR,MASTER,MPI_COMM_WORLD);

	// Let everyone sync up before we start
	MPI_Barrier(MPI_COMM_WORLD);
	fflush(stdout);

	///NEW QUEUE LOOP
	//In this new queue system, jobs waiting to be run are stored
	//in "jobs", which is a Red-Black Tree. Slave nodes waiting for
	//a job are stored in free_cpus.
	//Once a slave is free to work, it is popped from the free_cpu stack
	//and master sends the jobname and pointer to the job node for that
	//job (RBTree*) which is then stored in busy_list. Once
	//the slave finishes the work, it sends back the return buffer
	//and the job node pointer, which is then used by the master to
	//efficiently remove the finished job from the queue.
	if (master_node) {
	  int destination, sender, job_type, finished_dock_count;
	  int unpack_stic_status = 0;
	  char *temp_mopac_res;
	  finished_dock_count = 0;
	  temp_mopac_res = (char*)calloc(FILENAME_MAX,sizeof(char));
		jobs = CreateJobList(infilename,&params);
		busy_list = InitDeque(busy_list);
		SetupEffercio(&free_cpus,jobs,busy_list,&CompoundList,rank,size,&params);
		num_jobs = CountJobs(jobs) + CountJobs(busy_list);
		// printf_master("Will run %d jobs\n",num_jobs);
		printf_master("Will process %d ligands\n",num_jobs);


		// Wait for the results to come back and kick off another job
		// if necessary
		while(busy_list->head != NULL) {
			int return_val, have_incoming_msg;
			char *result_name;
			void *receive_buffer;
			// char *Compound_ID = malloc(sizeof(char *));
			char *Compound_ID = calloc(1,sizeof(char *));
			struct STICelement *results = NULL;
			size_t buffer_size;
			deque_node *received_node = NULL;
			job_t *received_job = NULL;

			//Save state
            #ifdef DEBUG
			printf_master("DEBUG - Saving state\n");fflush(stdout);
            #endif
			SaveState(jobs,busy_list,CompoundList,&params);
            #ifdef DEBUG
			printf_master("DEBUG - State saved\n");fflush(stdout);
            #endif

			results = InitSTIC(results);

			// Check to see if a message has arrived from a slave.
			// If so, go ahead an receive it.
			// If not, check to see if time has elapsed and shut everything down
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &have_incoming_msg,&status);
			while (!have_incoming_msg)
			{
				// struct timespec pause_interval;
				// pause_interval.tv_sec = 45;
				if(params.total_time - difftime(time(NULL),params.start_time) < 120)
				{
					// Time's up!
					printf_verbose(params.verbose,"Time remaining: %f\nShutting down.",
						       params.total_time - difftime(time(NULL),params.start_time));
				 	// Save the current state and shutdown 
					SaveState(jobs,busy_list,CompoundList,&params);
					MPI_Abort(MPI_COMM_WORLD,0);
				}

				// nanosleep(&pause_interval,NULL);
				MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &have_incoming_msg,&status);
			}


			// Get the results from the slave node
			// Data sent: Stic buffer size, Stic Buffer, job tree node pointer
			MPI_Recv(&buffer_size,1,MPI_INT,status.MPI_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD,&status);

			// Allocate the receive buffer
			// receive_buffer = malloc(buffer_size);
			receive_buffer = calloc(1,buffer_size);
			MPI_Recv(receive_buffer,buffer_size,MPI_BYTE,status.MPI_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Recv(&received_node,sizeof(deque_node*),MPI_BYTE,status.MPI_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            #ifdef DEBUG
			    printf("DEBUG - Receiving mopac_res\n");
            #endif
			memset(temp_mopac_res,0,FILENAME_MAX*sizeof(char));
			MPI_Recv(temp_mopac_res,FILENAME_MAX,MPI_CHAR,status.MPI_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			return_val = status.MPI_TAG;
			
            #ifdef DEBUG
			    printf("DEBUG - Received mopac_res for %s: \"%s\"\n",((job_t*)received_node->data)->name,temp_mopac_res);
            #endif
			if(strlen(temp_mopac_res) == 0)
			{
				free(((job_t*)received_node->data)->mopac_res);
				((job_t*)received_node->data)->mopac_res = (char*)calloc(FILENAME_MAX,sizeof(char));
			}
			else
			{
				free(((job_t*)received_node->data)->mopac_res);
				((job_t*)received_node->data)->mopac_res = strdup(temp_mopac_res);
			}
            if (params.verbose) {
			    printf("Filename copied for %s: \"%s\"\n",((job_t*)received_node->data)->name,
                       ((job_t*)received_node->data)->mopac_res);
            }

			suspend_cpu(&free_cpus,status.MPI_SOURCE);

			// Do something with the results here
			received_job = (job_t*)received_node->data;
			result_name = strdup(received_job->name);
			unpack_stic_status = UnpackBufferSTIC(receive_buffer, buffer_size,&result_name, &return_val, results);
			if(unpack_stic_status != 0)
			{
				// If unpacking fails, re-initialize the stic and warn the
				// user in the log file.
				int S,T,I,C;
				struct ClusterRep *tmp_reps;
				S = results->S;T = results->T;I = results->I;C = results->C;
				tmp_reps = results->reps;
				SetupSTIC(results);
				results->S = S;results->T = T;results->I = I;results->C = C;
				results->reps = tmp_reps;
				printf("ERROR - unpacking of stic buffer failed for %s\n",result_name);
			}

            #ifdef DEBUG
			    printf_master("DEBUG - STIC results received from slave are:\n");
			    FPrintSTIC(stdout,results);
                fflush(stdout);
            #endif
			received_job->input_data = *results;

			sender = status.MPI_SOURCE;
			if(strncmp(result_name,params.receptor_name,strlen(params.receptor_name)) == 0)
			  {
			    // TODO: Decide what to do if there is no ligand name. Then, add that check here.
			    size_t beg = strlen(params.receptor_name);
			    if(result_name[beg] == '_')
			      beg++;
			    GetStringByColumn(&result_name[beg],"_",1,&Compound_ID);
			  }
			else
			  GetStringByColumn(result_name,"_",1,&Compound_ID);			
            #ifdef DEBUG
			    printf("STIC for %s:\n",Compound_ID);
			    FPrintSTIC(stdout, results);
			    fflush(stdout);
            #endif

			// If all dock jobs have been completed, dump the CompoundList to a
			// file. If the file cannot be opened, dump it to standard to and
			// make a warning.
			//
			// Also, DOCK creates cluster reps. Therefore we want to initializes
			// their target energies from the input deck.
			if(received_node->data && received_node->type == DEQUE_JOB)
			  job_type = ((job_t*)received_node->data)->type;
			else
			  job_type = NUM_JOB_TYPES;

			if(job_type == DOCK)
			{
				struct ClusterRep *curr_rep = results->reps;
				while(curr_rep != NULL)
				{
					if(params.target_G != DOUBLE_INIT)
						curr_rep->optimized.G_prot = params.target_G;
					else
						curr_rep->optimized.G_prot = params.target_Hf;
					curr_rep = curr_rep->next;
				}
				finished_dock_count++;
			} 
            #ifdef DEBUG
                printf("DEBUG - temp_mopac_res prior to merge is \"%s\"\n",temp_mopac_res);
            #endif
			if(strlen(temp_mopac_res) == 0)// Only write stics when results finish, i.e. no restarts
			  MergeSTICData(Compound_ID,results,&CompoundList);

            if (params.verbose)
            {
			    printf_master("%s\n",AT_HDR);
			    printf_master("\nCOMPOUND LIST AFTER MERGE OF %s\n",result_name);
			    FPrintCompoundTree(stdout, CompoundList);
			    printf_master("%s\n",AT_HDR);
			    fflush(stdout);
            }

			finish_job(jobs,busy_list,received_node,receive_buffer,return_val,&params);

			if(job_type == DOCK && finished_dock_count == num_jobs)
			{
				char *dock_result_name = "Effercio.dockings";
				FILE *dock_results_file = fopen(dock_result_name,"w");
                #ifdef DEBUG
				printf("DEBUG - dumping dock data\n");
                #endif

				// Write the compound tree now that effercio has finished. If the
				// file cannot be opened, write it to standard out so that a least
				// the data will not be lost.
				if(dock_results_file == NULL)
				{
					printf("WARNING - Could not open %s. Reason: %s",dock_result_name, strerror(errno));
					printf("WARNING - Writing final results to standard output instead");
					dock_results_file = stdout;
				}
		                //BoltzmannAvgCompoundTree(CompoundList, params.analysis_dir,params.UseFreeEnergy);
				FPrintCompoundTree(dock_results_file, CompoundList);

				if(dock_results_file != stdout)
					fclose(dock_results_file);
			}

			//Are there jobs to be done and CPUs to do them?
			//If so, send out work.
			while (jobs->head != NULL && HaveFreeCPUs(free_cpus)){
				cpunode *free_node = NULL;
				free_node = PopCPUNode(&free_cpus);
				if(free_node != NULL)
				{
					issue_order(jobs,busy_list,free_node->rank,rank,&params);
					free(free_node);
				}
			}
			FreeSTIC(results);
			free(Compound_ID);
			free(result_name);
			free(receive_buffer);

		}

		printf_master("Finished with Job Queue. Sending quit signal to nodes.\n");
		//Send terminate queue loop command to node.
		char unnecessary[FILENAME_MAX];unnecessary[0] = 'a';unnecessary[1] = 0;
		for (destination = 0; destination < size; destination++)
		  if(destination != rank)
		    MPI_Send(unnecessary,FILENAME_MAX,MPI_CHAR,destination,QUIT_SIGNAL,MPI_COMM_WORLD);

		//Free Memory
        #ifdef DEBUG
            printf_master("Freeing memory pointers\n");
        #endif
		if (jobs != NULL) free(jobs);
		if (busy_list != NULL) free(busy_list);
		FreeCPUNode(free_cpus);
		free(temp_mopac_res);
	}
	else { // Slave node
		char jobname[FILENAME_MAX];
		//char mopac_res[FILENAME_MAX];
		char *return_buffer;
		deque_node *current_job = NULL;
		double float_buff[FLOAT_TRANSFER_COUNT];
		job_t newjob;
		newjob.name = jobname;
		//newjob.mopac_res = mopac_res;
		//strcpy(mopac_res,"");
                newjob.mopac_res = calloc(FILENAME_MAX,sizeof(char));

		return_buffer = calloc(BUFFER_SIZE,sizeof(char));

		//Get initial job. If the MPI_TAG is QUIT_SIGNAL, do no work.
		MPI_Recv(jobname,FILENAME_MAX,MPI_CHAR,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		if(status.MPI_TAG != QUIT_SIGNAL)
		{
			MPI_Recv(&current_job,sizeof(deque_node*),MPI_BYTE,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Recv(float_buff,FLOAT_TRANSFER_COUNT,MPI_DOUBLE,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			MPI_Recv(newjob.mopac_res,FILENAME_MAX,MPI_CHAR,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		}


		//Continuously do work until the master node issues the Quit Signal
		while(status.MPI_TAG != QUIT_SIGNAL) {
			size_t num_bytes = 0;
			int retval = 0;
			struct STICelement *STIC = NULL;

			STIC = InitSTIC(STIC);

			SetupSTIC(newjob.input_data);
			newjob.input_data.Hf = STIC->Hf = float_buff[0];
			newjob.input_data.G = STIC->G = float_buff[1];
			newjob.type = status.MPI_TAG;

            #ifdef DEBUG
			    char strJobType[512];
			    JobString(strJobType,&newjob);
			    printf("%s: DEBUG - Got %s from master node for %s job\n",params.node_tag,jobname,strJobType);
			    printf("%s: DEBUG - status.MPI_TAG = %d\n",params.node_tag,status.MPI_TAG);
			    printf("%s: DEBUG - charge on STIC before = %.1f\n",params.node_tag,STIC->charge);
            #endif

			// Run the calculation here
			retval = RunJob(&newjob,STIC,&params);
			if(retval < 0)
			  retval *= -1;
            #ifdef DEBUG
			    printf("%s: DEBUG - charge on STIC after = %.1f\n",params.node_tag,STIC->charge);
			    printf("%s: DEBUG - %s results returned:\n",params.node_tag,strJobType);
			    FPrintSTIC(stdout, STIC);
            #endif

			PackBufferSTIC(&return_buffer,&num_bytes,jobname,retval,STIC);

			MPI_Send(&num_bytes,1,MPI_INT,MASTER,retval,MPI_COMM_WORLD);
			MPI_Send(return_buffer,num_bytes,MPI_BYTE,MASTER,retval,MPI_COMM_WORLD);
			MPI_Send(&current_job,sizeof(RBTree*),MPI_BYTE,MASTER,retval,MPI_COMM_WORLD);
			//MPI_Send(newjob.mopac_res,sizeof(char)*(strlen(newjob.mopac_res)+1),MPI_CHAR,MASTER,retval,MPI_COMM_WORLD);
			MPI_Send(newjob.mopac_res,FILENAME_MAX,MPI_CHAR,MASTER,retval,MPI_COMM_WORLD);



			//Receive either next job or order to shutdown.
			MPI_Recv(jobname,FILENAME_MAX,MPI_CHAR,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			if(status.MPI_TAG != QUIT_SIGNAL)
			{
				MPI_Recv(&current_job,sizeof(RBTree*),MPI_BYTE,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				MPI_Recv(float_buff,FLOAT_TRANSFER_COUNT,MPI_DOUBLE,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				newjob.mopac_res[0] = (char)0;
				MPI_Recv(newjob.mopac_res,FILENAME_MAX,MPI_CHAR,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			}

			FreeSTIC(STIC);
			memset(return_buffer,0,num_bytes);

		}

        free(newjob.mopac_res);
		free(return_buffer);

        #ifdef DEBUG
		    printf("%s: DEBUG - No more jobs for me to do\n",params.node_tag);
        #endif

	}
	//END NEW QUEUE LOOP

    //Save state
    if (params.verbose) {
        printf_master("DEBUG - Saving state\n");fflush(stdout);
    }
    SaveState(jobs,busy_list,CompoundList,&params);
    #ifdef DEBUG
        printf_master("DEBUG - State saved\n");fflush(stdout);
    #endif

	if (master_node)
	{
		char *full_result_name = "Effercio.results";
		char *sql_result_name = "Effercio.sql";
		FILE *full_results_file = fopen(full_result_name,"w");
		FILE *sql_results_file = NULL;

		// Write the compound tree now that effercio has finished. If the
		// file cannot be opened, write it to standard out so that a least
		// the data will not be lost.
		if(full_results_file == NULL)
		{
			printf("WARNING - Could not open %s. Reason: %s\n",full_result_name, strerror(errno));
			printf("WARNING - Writing final results to standard output instead\n");
			full_results_file = stdout;
		}
		else
		{
			printf("Writing final results to %s\n",full_result_name);
		}

		// Average results and write to the results file
		printf("Performing Boltzmann averaging ... ");
		BoltzmannAvgCompoundTree(CompoundList, params.analysis_dir,params.UseFreeEnergy);
		printf("done.\n");
		FPrintCompoundTree(full_results_file, CompoundList);
		
		// Write SQL statements for the results
		sql_results_file = fopen(sql_result_name,"w");
		if(sql_results_file == NULL)
		  {
		    printf("WARNING - Could not open %s. Reason: %s\n",sql_result_name, strerror(errno));
		    printf("WARNING - Writing SQL statements to standard output instead\n");
		    sql_results_file = stdout;
		  }
		else
		{
			printf("Writing SQL statements to %s\n",sql_result_name);
		}
		FPrintSQL(sql_results_file,CompoundList,&params);

		FreeRBTree(CompoundList);
		if(full_results_file != stdout)
		  fclose(full_results_file);
		if(sql_results_file != stdout)
		  fclose(sql_results_file);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (master_node) {
		time_t end_time = time(NULL);
		printf("\n%s finished on: %s\n",PROG_NAME,asctime(localtime(&end_time)));
		timestring = SecondsToHMS(difftime(end_time,params.start_time));
		printf("Total time = %s\n\n",timestring);
		free(timestring);
	}

	// Should probably have a function to free params rather than this piecemeal approach
	FreeStringArray(&params.dock_params,&params.num_dock_params);
	FreeStringArray(params.mopac_footer_params,&params.num_mopac_footer_params);
	FreeStringArray(params.mopac_header_params,&params.num_mopac_header_params);

	MPI_Finalize();

	return 0;
}
