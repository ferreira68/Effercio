#include "state.h"
#include "defines.h"
#include <unistd.h>

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
			retval = UnpackFileSTIC(stic_tpl,&new_compound.ID, 0, new_stic);
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


