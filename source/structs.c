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
 * This sets up a simple stack with push and pop operations to contain   *
 * the list of drug names to be processed.  It lacks significant error   *
 * checking, but should get the basic job done.                          *
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

#include <limits.h>
#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defines.h"
#include "structs.h"
#include "io_utils.h"
#include "deque.h"
#include "RBTree.h"
#include "mpi.h"
#include "memwatch.h"

// Defined in analysis.c
extern void BoltzmannAvgSTIC(struct STICelement *stic, double *Z_list_DOCK, double *Z_list_QM, double *Z_DOCK, double *Z_QM, double *dock,double *qm);
extern void AdjustUnits(double *val, char *units);

/**
 * This is a stack for cpu nodes, represented by their MPI rank.
 *
 * Returns a pointer to the newly created cpu node for the
 * node corresponding to rank.
 */
cpunode* PushCPUNode(cpunode **stack,int rank)
{
	cpunode* newnode = (cpunode *)calloc(1,sizeof(cpunode));
	newnode->rank = rank;
	if(*stack == NULL)
	{
		newnode->link = NULL;
		*stack = newnode;
	}
	else
	{
		newnode->link = *stack;
		*stack = newnode;
	}
	return newnode;
}

/**
 * Returns the first (i.e. top) value on the stack and removes
 * that node from the stack.
 */
cpunode* PopCPUNode(cpunode **stack)
{
	cpunode *popped_value;
	if(*stack == NULL)
		return NULL;
	popped_value = *stack;
	*stack = popped_value->link;
	popped_value->link = NULL;

	return popped_value;
}

/**
 * Indicates whether or not there are free cpus to which
 * work may be sent.
 */
int HaveFreeCPUs(const cpunode *node_list)
{
	return (node_list != NULL);
}

/**
 * Frees the memory allocated to the cpunode and all of the
 * cpunodes to which it links.
 */
void FreeCPUNode(cpunode *cpu)
{
	if(cpu == NULL)
		return;
	if(cpu->link != NULL)
		free(cpu->link);
	free(cpu);
}

/**
 * Counts the number of jobs contained in the job tree (downward).
 */
int CountJobs(const deque *jobs)
{
	return deque_size(jobs);
}

char* JobString(char* string, const job_t* job)
{
  if(string == NULL)
    return string;
  if(job == NULL)
    return strcpy(string,"NULL");
  else
    {
      switch(job->type)
	{
	case QM:
	  return strcpy(string, "Optimization");
	  break;
	case QMLIG:
	  return strcpy(string, "Optimization (Ligand only)");
	  break;
	case PRESCREENING:
	  return strcpy(string,"Prescreening");
	  break;
	case DOCK:
	  return strcpy(string,"Docking");
	  break;
	default:
	  return strcpy(string,"Unknown type");
	  break;
	}
    }
  return string;
}

/**
 * Determines whether or not the tree is empty.
 */
int IsEmpty(const RBTree *job_tree)
{
	return job_tree == NULL;
}

/**
 * Prints the address and name of the jobs contained
 * in the tree.
 */
void PrintJoblist(const RBTree *jobs)
{
	if(jobs == NULL)
		return;
	PrintJoblist(jobs->left);
	printf("FILL THIS IN!");//0x%x -> %s\n",jobs->data->name,jobs->data->name);
	PrintJoblist(jobs->right);
}

/**
 * Converts data to kcals
 */
void ConvertToKCals(double *data, int original_units)
{
	if(data == NULL)
		return;

	switch(original_units)
	{
	case INPUT_J:
		*data *= 0.000239005736;
		break;
	case INPUT_KJ:
		*data *= 0.239005736;
		break;
	case INPUT_CAL:
		*data *= 0.001;
		break;
	default:
		break;
	}
}



void SetupDOCK(struct DOCKresult *data)
{
	const char InitString[] = STR_INIT;
	const int InitStrLen = strlen(InitString) + 1;

    data->G_binding = (double) DOUBLE_INIT;
    data->Ki_DOCK = (double) DOUBLE_INIT;

    data->Ki_unit = calloc((size_t) InitStrLen,sizeof(char));
    strncpy(data->Ki_unit,InitString,InitStrLen);

    data->E_inter = (double) DOUBLE_INIT;
    data->E_nonbond = (double) DOUBLE_INIT;
    data->E_electrostat = (double) DOUBLE_INIT;
    data->E_internal = (double) DOUBLE_INIT;
    data->G_tors = (double) DOUBLE_INIT;
    data->E_unbound = (double) DOUBLE_INIT;
    data->rmsd_ref = (double) DOUBLE_INIT;
    data->time = (double) DOUBLE_INIT;
}


void FPrintDOCK(FILE *file, const struct DOCKresult *data)
{
    fprintf(file, "    AutoDock Results\n");
    fprintf(file, "%s\n",DASH_HDR);
    if (data->G_binding != DOUBLE_INIT)
        fprintf(file, "        Free energy of binding               = %6.2f kcal/mol\n",data->G_binding);
    if (data->Ki_DOCK != DOUBLE_INIT)
        fprintf(file, "        Estimated Ki                         = %6.2f %s\n",data->Ki_DOCK,data->Ki_unit);
    if (data->E_inter != DOUBLE_INIT)
        fprintf(file, "        Intermolecular interaction energy    = %6.2f kcal/mol\n",data->E_inter);
    if (data->E_nonbond != DOUBLE_INIT)
        fprintf(file, "        Nonbonded interaction energy         = %6.2f kcal/mol\n",data->E_nonbond);
    if (data->E_electrostat != DOUBLE_INIT)
        fprintf(file, "        Electrostatic energy                 = %6.2f kcal/mol\n",data->E_electrostat);
    if (data->E_internal != DOUBLE_INIT)
        fprintf(file, "        Internal energy                      = %6.2f kcal/mol\n",data->E_internal);
    if (data->G_tors != DOUBLE_INIT)
        fprintf(file, "        Torsional free energy                = %6.2f kcal/mol\n",data->G_tors);
    if (data->E_unbound != DOUBLE_INIT)
        fprintf(file, "        Energy of the unbound system         = %6.2f kcal/mol\n",data->E_unbound);
    if (data->rmsd_ref != DOUBLE_INIT)
        fprintf(file, "        RMSD from the input structure        = %6.2f Angstrom\n",data->rmsd_ref);
    if (data->time != DOUBLE_INIT)
        fprintf(file, "        Time to perform AutoDock calculation = %6.2f seconds\n",data->time);
    fprintf(file, "\n\n");

    fflush(file);

    return;
}


void SetupQM(struct QMresult *data)
{
    data->Hf = (double) DOUBLE_INIT;
    data->S = (double) DOUBLE_INIT;
    data->Cp = (double) DOUBLE_INIT;
    data->G = (double) DOUBLE_INIT;
    data->ZPE = (double) DOUBLE_INIT;
    data->E_dielec = (double) DOUBLE_INIT;
    data->mu_x = (double) DOUBLE_INIT;
    data->mu_y = (double) DOUBLE_INIT;
    data->mu_z = (double) DOUBLE_INIT;
    data->mu_total = (double) DOUBLE_INIT;
    data->time = (double) DOUBLE_INIT;

/*     data->method = NULL; */
/*     data->method = realloc(data->method,strlen(InitString)+1); */
/*     strncpy(data->method,InitString,strlen(InitString)+1); */

    data->num_SCFs = INT_INIT;
    data->COSMO_A = (double) DOUBLE_INIT;
    data->COSMO_V = (double) DOUBLE_INIT;
    data->vdW_A = (double) DOUBLE_INIT;
    data->G_prot = (double) DOUBLE_INIT;
    //data->G_liginput = (double) DOUBLE_INIT;
    data->G_ligand = (double) DOUBLE_INIT;
    data->G_binding = (double) DOUBLE_INIT;
    data->Ki_QM = (double) DOUBLE_INIT;
/*     data->Ki_type = NULL; */
/*     data->Ki_type = realloc(data->Ki_type,strlen(InitString)+1); */
/*     strncpy(data->Ki_type,InitString,strlen(InitString)+1); */

    data->method = strdup(STR_INIT);
    data->Ki_type = strdup(STR_INIT);

}


void FPrintQM(FILE *file, const struct QMresult *data)
{
    // NOTE: This should conditionally add all these things to a single character buffer and then
    //       print out if it's non-zero in length.
    fprintf(file, "    MOPAC Results\n");
    fprintf(file, "%s\n",DASH_HDR);
    if (strcmp(data->method,STR_INIT) != 0)
	    fprintf(file, "        Semi-empirical Method       = %s\n",data->method);
    if (data->Hf != DOUBLE_INIT)
	    fprintf(file, "        Heat of Formation           = % .5f kcal/mol\n",data->Hf);
    if (data->S != DOUBLE_INIT)
	    fprintf(file, "        Entropy                     = % .6G kcal/mol\n",data->S);
    if (data->Cp != DOUBLE_INIT)
	    fprintf(file, "        Heat Capacity (Cp)          = % .6G kcal/mol\n",data->Cp);
    if (data->G != DOUBLE_INIT)
	    fprintf(file, "        Free Energy                 = % .4f kcal/mol\n",data->G);
    if (data->ZPE != DOUBLE_INIT)
	    fprintf(file, "        Zero Point Energy           = % .3f kcal/mol\n",data->ZPE);
    if (data->E_dielec != DOUBLE_INIT)
	    fprintf(file, "        Dielectric Energy           = % .6f kcal/mol\n",data->E_dielec);
    if (data->COSMO_A != DOUBLE_INIT)
	    fprintf(file, "        COSMO Surface Area          = %.2f sq. Angstrom\n",data->COSMO_A);
    if (data->COSMO_V != DOUBLE_INIT)
	    fprintf(file, "        COSMO Volume                = %.2f cu. Angstrom\n",data->COSMO_V);
    if (data->vdW_A != DOUBLE_INIT)
	    fprintf(file, "        van der Waals Surface Area  = %.2f sq. Angstrom\n",data->vdW_A);
    if (data->mu_x != DOUBLE_INIT && data->mu_y != DOUBLE_INIT && data->mu_z != DOUBLE_INIT) {
	    fprintf(file, "        Dipole Moment\n");
	    fprintf(file, "            x-component             = % .3f Debye\n",data->mu_x);
	    fprintf(file, "            y-component             = % .3f Debye\n",data->mu_y);
	    fprintf(file, "            z-component             = % .3f Debye\n",data->mu_z);
	    fprintf(file, "            -----------------------------------------\n");
	    fprintf(file, "                                    = % .3f Debye\n",data->mu_total);
    }
    if (data->time != DOUBLE_INIT)
	    fprintf(file, "        Calculation Time            = %.2f (seconds)\n",data->time);
    if (data->num_SCFs != INT_INIT)
	    fprintf(file, "        Number of SCF iterations    = %d\n",data->num_SCFs);
    if (strcmp(data->Ki_type,STR_INIT) != 0)
	    fprintf(file, "        Ki Estimation Method        = %s\n",data->Ki_type);
    if (data->G_prot != DOUBLE_INIT)
	    fprintf(file, "        Free Energy of Receptor     = %.6f kcal/mol\n",data->G_prot);
    if (data->G_ligand != DOUBLE_INIT)
	    fprintf(file, "        Free Energy of Ligand       = %.6f kcal/mol\n",data->G_ligand);
    if (data->G_binding != DOUBLE_INIT) {
        if (data->G_binding > ZERO) {
	        fprintf(file, "        Free Energy of Binding      = %.6f kcal/mol (Does Not Bind)\n",data->G_binding);
	        fprintf(file, "        Inhibition Constant (Ki)    = %.6G M\n",INFINITY);
        }
        else {
	        fprintf(file, "        Free Energy of Binding      = %.6f kcal/mol\n",data->G_binding);
	        fprintf(file, "        Inhibition Constant (Ki)    = %.6G M\n",data->Ki_QM);
        }
    }
    if (data->Ki_QM != DOUBLE_INIT)
    fprintf(file, "\n\n");

    fflush(file);

    return;
}


struct ClusterRep *InitClusterRep(struct ClusterRep *data)
{
    data = (struct ClusterRep *) malloc(sizeof(struct ClusterRep));

/*     data->index = 1; */
/*     data->size  = 1; */
    data->index = INT_INIT;
    data->size  = INT_INIT;
//	const char InitString[] = STR_INIT;
//	const int InitStrLen = strlen(InitString) + 1;
//	data->DOCK_Ki_unit = calloc((size_t) InitStrLen,sizeof(char));
//    strncpy(data->DOCK_Ki_unit,InitString,InitStrLen);
//    data->QM_method = strdup(STR_INIT);
//    data->QM_Ki_type = strdup(STR_INIT);

	SetupQM(&data->optimized);
	SetupDOCK(&data->docked);
    data->next = NULL;
    return data;
}


void FreeDOCK(struct DOCKresult *data)
{
  free(data->Ki_unit);
}

void FreeQM(struct QMresult *data)
{
	free(data->method);
	free(data->Ki_type);
}

void FreeClusterRep(struct ClusterRep *data)
{
  if(data == NULL)
    return;
  if (data->next != NULL) {
    FreeClusterRep(data->next);
  }

  FreeDOCK(&data->docked);
  FreeQM(&data->optimized);
  free(data);
}


void FPrintClusterRep(FILE *file, const struct ClusterRep *data)
{
  if(data == NULL)
    return;
  if(file == NULL)
    file = stdout;

  fprintf(file, "%s\n",WIGL_HDR);
  fprintf(file, "Data for cluster number %d\n",data->index);
  fprintf(file, "=== Number of poses in cluster: %d\n\n",data->size);
  FPrintDOCK(file, &data->docked/*,data->DOCK_Ki_unit*/);
  FPrintQM(file, &data->optimized/*,data->QM_method,data->QM_Ki_type*/);
  fprintf(file, "%s\n",WIGL_HDR);

  fflush(file);
  return;
}

void SetupSTIC(struct STICelement *data)
{
	if(data == NULL)
		return;
	data->S = INT_INIT;
	data->T = INT_INIT;
	data->I = INT_INIT;
	data->C = INT_INIT;
	data->charge = (double) DOUBLE_INIT;
	data->total_charge = (double) DOUBLE_INIT;
	data->Hf = (double) DOUBLE_INIT;
	data->G = (double) DOUBLE_INIT;
	data->reps = NULL;
}

struct STICelement *InitSTIC(struct STICelement *data)
{
	data = (struct STICelement *) malloc(sizeof(struct STICelement));

	SetupSTIC(data);

	return data;
}


void FreeSTIC(struct STICelement *data)
{
  if(data == NULL)
    return;
  FreeClusterRep(data->reps);
  free(data);
  return;
}


void FPrintSTIC(FILE *file, struct STICelement *data)
{
    struct ClusterRep *clustptr;
    double coeff_dock[4096];
    double coeff_qm[4096];
    double Z_DOCK = ZERO;
    double Z_QM = ZERO;
    double Ki_docking = ZERO;
    double Ki_qm      = ZERO;
    double adjusted_value = ZERO;
    char units[16];

    BoltzmannAvgSTIC(data,coeff_dock,coeff_qm,&Z_DOCK,&Z_QM,&Ki_docking,&Ki_qm);

    if(file == NULL)
    	file = stdout;

    fprintf(file, "%s\n",EQ_HDR);
    fprintf(file, "Stereoisomer index      = %d\n",data->S);
    fprintf(file, "Tautomer index          = %d\n",data->T);
    fprintf(file, "Ionization state index  = %d\n",data->I);
    fprintf(file, "Conformer index         = %d\n",data->C);
    fprintf(file, "Charge on ligand        = %.2f\n",data->charge);
    if (data->total_charge != DOUBLE_INIT)
    {
        fprintf(file, "Total charge on complex = %.2f\n",data->total_charge);
    }
    else
    {
        fprintf(file, "Total charge on complex = unavailable\n");
    }
    fprintf(file, "Heat of Formation       = %.5f kcal/mol\n",data->Hf);
    if (data->G != DOUBLE_INIT)
    {
        fprintf(file, "Free Energy             = %.5f kcal/mol\n\n",data->G);
    }
    else
    {
        fprintf(file, "Free Energy             = unavailable\n\n");
    }
    if (Ki_docking != DOUBLE_INIT) 
    {
        adjusted_value = Ki_docking;
        AdjustUnits(&adjusted_value,units);
        fprintf(file, "Mean K_i (docked)       = %.6g %sM\n",adjusted_value,units);
    }
    else
    {
        fprintf(file, "Mean K_i (docked)       = unavailable\n");
    }
    if (Ki_qm != DOUBLE_INIT)
    {
        adjusted_value = Ki_qm;
        AdjustUnits(&adjusted_value,units);
        fprintf(file, "Mean K_i (optimized)    = %.6g %sM\n",adjusted_value,units);
    }
    else
    {
        fprintf(file, "Mean K_i (optimized)    = unavailable\n");
    }
    fprintf(file, "Number of clusters      = %d\n",NumClusterReps(data));

    clustptr = data->reps;
    while (clustptr != NULL) {
        FPrintClusterRep(file, clustptr);
        clustptr = clustptr->next;
        if (clustptr != NULL) fprintf(file,"\n");
    }

    fprintf(file, "%s\n",EQ_HDR);

    fflush(file);
    return;
}


int NumClusterReps(const struct STICelement *data)
{
    struct ClusterRep *clustptr;
    int numreps = 0;

    clustptr = data->reps;
    while (clustptr != NULL) {
        numreps++;
        clustptr = clustptr->next;
    }

    return numreps;
}


struct compoundData *InitCompound(struct compoundData *data)
{
    data = (struct compoundData *) malloc(sizeof(struct compoundData));

    data->ID = calloc((size_t) strlen(STR_INIT)+1,sizeof(char));
    strncpy(data->ID,STR_INIT,strlen(STR_INIT)+1);
    data->num_STIC = 0;
    return data;
}

void FreeCompound(struct compoundData *data)
{
  if(data != NULL)
    return;
  free(data->ID);
  free(data);
  return;
}


void FPrintCompounds(FILE *file, const struct compoundData *data)
{
	if(file == NULL)
		file = stdout;
	fprintf(file,"\n");
    fflush(file);
    if (data != NULL) {
        fprintf(file,"%s\n",PCNT_HDR);
        fprintf(file,"Compound ID = %s\n",data->ID);
        fprintf(file,"Number of STICs = %d\n",data->num_STIC);
        fprintf(file,"%s\n",PCNT_HDR);
     }

    fflush(file);
    return;
}

void FPrintSTICTree(FILE *file, const RBTree *data)
{
	if(file == NULL)
			file = stdout;
	if(data == NULL)
		return;
	if(data->type != STIC)
	{
		fprintf(stderr,"PrintSTICTree attempted on a non stic tree.\n");
		exit(RBTREE_TYPE_ERROR);
	}

	FPrintSTICTree(file,data->left);

	FPrintSTIC(file,data->data);

	FPrintSTICTree(file,data->right);
}

/**
 * Prints the Compound Tree to the specified file.
 * If file is NULL, the data is written to standard output.
 */
void FPrintCompoundTree(FILE *file, const RBTree *data)
{
	if(file == NULL)
		file = stdout;

	if(data == NULL)
		return;
	if(data->type != COMPOUND)
	{
		fprintf(stderr,"PrintCompoundTree attempted on a non compound tree.\n");
		exit(RBTREE_TYPE_ERROR);
	}
	FPrintCompoundTree(file,data->left);

	FPrintCompounds(file,((CompoundTree*)data->data)->data);
	FPrintSTICTree(file,((CompoundTree*)data->data)->stics);

	FPrintCompoundTree(file,data->right);
}

char *BufferAddInt(char *buffptr, size_t *numbytes, int value)
{
    size_t BUFFER_SIZE = (size_t) (MSG_BUFSIZ * MEM_BLOCK_SIZE);

    if ((*numbytes + sizeof(int)) > BUFFER_SIZE) {
        printf("ERROR - Unable to add value to transfer buffer!\n");
        printf("ERROR - Try increasing MSG_BUFSIZ and recompiling.\n");
        printf("ERROR - Current buffer size is %lu bytes\n",BUFFER_SIZE);
        exit(-1);
    }
    memcpy(buffptr,&value,sizeof(int));
    buffptr += sizeof(int);
    *numbytes += sizeof(int);

    return buffptr;
}


char *BufferGetInt(char *buffptr, int *value)
{
    memcpy(value,buffptr,sizeof(int));
    buffptr += sizeof(int);

    return buffptr;
}


char *BufferAddDbl(char *buffptr, size_t *numbytes, double value)
{
    size_t BUFFER_SIZE = (size_t) (MSG_BUFSIZ * MEM_BLOCK_SIZE);

    if ((*numbytes + sizeof(double)) > BUFFER_SIZE) {
        printf("ERROR - Unable to add value to transfer buffer!\n");
        printf("ERROR - Try increasing MSG_BUFSIZ and recompiling.\n");
        printf("ERROR - Current buffer size is %lu bytes\n",BUFFER_SIZE);
        exit(-1);
    }
    memcpy(buffptr,&value,sizeof(double));
    buffptr += sizeof(double);
    *numbytes += sizeof(double);

    return buffptr;
}


char *BufferGetDbl(char *buffptr, double *value)
{
    memcpy(value,buffptr,sizeof(double));
    buffptr += sizeof(double);

    return buffptr;
}


char *BufferAddStr(char *buffptr, size_t *numbytes, char *string)
{
    size_t BUFFER_SIZE = (size_t) (MSG_BUFSIZ * MEM_BLOCK_SIZE);
    size_t mystrlen;

    if ((*numbytes + sizeof(size_t)) > BUFFER_SIZE) {
        printf("ERROR - Unable to add value to transfer buffer!\n");
        printf("ERROR - Try increasing MSG_BUFSIZ and recompiling.\n");
        printf("ERROR - Current buffer size is %lu bytes\n",BUFFER_SIZE);
        exit(-1);
    }
    mystrlen = strlen(string) + 1;
    memcpy(buffptr,&mystrlen,sizeof(size_t));
    buffptr += sizeof(size_t);
    *numbytes += sizeof(size_t);
    if ((*numbytes + mystrlen*sizeof(char)) > BUFFER_SIZE) {
        printf("ERROR - Unable to add value to transfer buffer!\n");
        printf("ERROR - Try increasing MSG_BUFSIZ and recompiling.\n");
        printf("ERROR - Current buffer size is %lu bytes\n",BUFFER_SIZE);
        exit(-1);
    }
    memcpy(buffptr,string,mystrlen * sizeof(char));
    buffptr += mystrlen * sizeof(char);
    *numbytes += mystrlen * sizeof(char);
	
    return buffptr;
}


char *BufferGetStr(char *buffptr, char **string)
{
    size_t mystrlen;

    memcpy(&mystrlen,buffptr,sizeof(size_t));
    buffptr += sizeof(size_t);

    *string = realloc(*string,mystrlen * sizeof(char));
    memcpy(*string,buffptr,mystrlen * sizeof(char));
    string = &(*string);
    buffptr += mystrlen * sizeof(char);

    return buffptr;
}


int PackBufferSTIC(void **buffer, size_t *size, const char *jobname, int jobretval,
		   struct STICelement *data)
{
    int retval = 0;
    struct ClusterRep curr_rep, *curr_pointer = data->reps;
    tpl_node *tn = tpl_map(STIC_TPL_MAP,data,&curr_rep);

    if(tn != NULL)
    {
    	tpl_pack(tn,1);//Pack STIC data
    	//Iterate through cluster reps
    	while(curr_pointer!= NULL)
    	{
    		curr_rep = *curr_pointer;
    		tpl_pack(tn,2);
    		curr_pointer = curr_pointer->next;
    	}
    	tpl_dump(tn,TPL_MEM,buffer,size);
    	tpl_free(tn);tn = NULL;
    }
    else
    	retval = TPL_ERROR;

    return retval;
}


int UnpackBufferSTIC(void *buffer, size_t size,char **jobname, int *jobretval, struct STICelement *data)
{
	int retval = 0;
	struct ClusterRep curr_rep,*curr_pointer = NULL;
	tpl_node *tn = tpl_map(STIC_TPL_MAP,data,&curr_rep);

	if(tn != NULL && tpl_load(tn,TPL_MEM,buffer,size) == 0)
	{
		tpl_unpack(tn,1);//Pack STIC data
		//Unpack first Cluster rep (assuming there is one)
		if(tpl_unpack(tn,2) > 0)
		{
			curr_rep.next = NULL;
			data->reps = InitClusterRep(data->reps);
			curr_pointer = data->reps;
			*curr_pointer = curr_rep;
		}
		while(tpl_unpack(tn,2) > 0)
		{
			curr_rep.next = NULL;
			curr_pointer->next = InitClusterRep(curr_pointer->next);
			curr_pointer = curr_pointer->next;
			*curr_pointer = curr_rep;
		}
		tpl_free(tn);tn = NULL;
	}
	else
		retval = TPL_ERROR;

	return retval;
}

int UnpackFileSTIC(FILE *file,char **jobname, int *jobretval, struct STICelement *data)
{
	int retval = 0;
	struct ClusterRep curr_rep,*curr_pointer = NULL;
	tpl_node *tn = tpl_map(STIC_TPL_MAP,data,&curr_rep);

	if(tn != NULL && tpl_load(tn,TPL_FD,fileno(file)) == 0)
	{
		tpl_unpack(tn,1);//Pack STIC data
		//Unpack first Cluster rep (assuming there is one)
		if(tpl_unpack(tn,2) > 0)
		{
			curr_rep.next = NULL;
			data->reps = InitClusterRep(data->reps);
			curr_pointer = data->reps;
			*curr_pointer = curr_rep;
		}
		while(tpl_unpack(tn,2) > 0)
		{
			curr_rep.next = NULL;
			curr_pointer->next = InitClusterRep(curr_pointer->next);
			curr_pointer = curr_pointer->next;
			*curr_pointer = curr_rep;
		}
		tpl_free(tn);tn = NULL;
	}
	else
		retval = TPL_ERROR;

	return retval;
}


int MergeSTICData(char *compound_name, struct STICelement *data, RBTree **molecule_list)
{
	RBTree *CompoundNode, *STICNode, *STICList;
	struct STICelement *STICptr = NULL;
	struct ClusterRep *curr_rep,*data_rep;
	int stic_is_new, compound_is_new = 0;
	int retval = 0;

	STICNode = CompoundNode = NULL;

	//Are there compounds?
	//If there are compounds, get a pointer to the compounds
	//whose name matches compound_name
	//Else, insert a new compound
	if(*molecule_list != NULL)
	{
		struct compoundData *tmp = NULL;
                char *name_holder = NULL; // This is used to temporarily hold the calloc'ed ID from InitCompound.
		tmp = InitCompound(tmp);
                name_holder = tmp->ID;
		tmp->ID = compound_name;//ok, since tmp is temporary for the search
		CompoundNode = FindRBTNode(*molecule_list,tmp);
		compound_is_new = (CompoundNode == NULL);
		tmp->ID = name_holder;
                FreeCompound(tmp);
	}
	else
		compound_is_new = 1;

	if(compound_is_new)
	{
		size_t name_size = strlen(compound_name) + 1;
		struct compoundData *cmpd = NULL;
                CompoundTree *new_compound_leaf = malloc(sizeof(CompoundTree));

		CompoundNode = InitRBTree(CompoundNode,COMPOUND);
		CompoundNode->data = new_compound_leaf;

                // Initialize compound and copy its name into the structure.
		cmpd = InitCompound(cmpd);
                free(cmpd->ID); // not to leak this memory. It's calloc'ed in InitCompound.
		strncpy(cmpd->ID,compound_name,name_size);

                // Assign
                new_compound_leaf->stics = NULL;
                new_compound_leaf->data = cmpd;
		CompoundNode->data = new_compound_leaf;
		insertRBT(molecule_list,CompoundNode);
	}

	// Get STIC Pointer.
	// If CompoundNode already has a similar STIC
	// (based on RBTDirection sort criteria),
	// then update it.
	// Else create a new one with this data.
	STICList = ((CompoundTree*)CompoundNode->data)->stics;
	if(STICList == NULL)
		stic_is_new = 1;
	else
	{
		STICNode = FindRBTNode(STICList,data);
		stic_is_new = (STICNode == NULL);
	}

	if(stic_is_new)
	{
		((CompoundTree*)CompoundNode->data)->data->num_STIC++;
		STICptr = InitSTIC(STICptr);
		STICNode = InitRBTree(STICNode,STIC);
		STICNode->data = STICptr;
		*STICptr = *data;
                STICptr->reps = NULL;
		((CompoundTree*)CompoundNode->data)->stics = insertRBT(&((CompoundTree*)CompoundNode->data)->stics,STICNode);
	}
	else
	{
		STICptr = (struct STICelement*)STICNode->data;
	}

	// If the stic in the compound list has default values,
	// copy that STIC data from the new stic
	if (STICptr->S == INT_INIT) STICptr->S = data->S;
	if (STICptr->T == INT_INIT) STICptr->T = data->T;
	if (STICptr->I == INT_INIT) STICptr->I = data->I;
	if (STICptr->C == INT_INIT) STICptr->C = data->C;
	if (STICptr->charge == DOUBLE_INIT) STICptr->charge = data->charge;
	if (STICptr->total_charge == DOUBLE_INIT) STICptr->total_charge = data->total_charge;
	if (STICptr->Hf == DOUBLE_INIT) STICptr->Hf = data->Hf;
	if (STICptr->G == DOUBLE_INIT) STICptr->G = data->G;

	// Iterate through the new stic's cluster reps and insert/update
	// as needed
	data_rep = data->reps;
	while(data_rep != NULL)
	{
		//Find (or create) the correct cluster rep in the stic
		curr_rep = STICptr->reps;
		while(curr_rep != NULL && curr_rep->index != data_rep->index)
		{
			if(curr_rep->index == INT_INIT)
				break;
			if(curr_rep->next == NULL)
			{
				curr_rep->next = InitClusterRep(curr_rep->next);
				curr_rep = curr_rep->next;
				break;
			}
			curr_rep = curr_rep->next;
		}

		if(curr_rep == NULL)
		{
			curr_rep = InitClusterRep(curr_rep);
			STICptr->reps = curr_rep;
		}

		if(curr_rep->index == INT_INIT)
			curr_rep->index = data_rep->index;
		if(curr_rep->size == INT_INIT)
			curr_rep->size = data_rep->size;

		//Update DOCK data
		struct DOCKresult *dock_ptr, *dock_data;
		size_t newlen,oldlen;
		dock_ptr = &curr_rep->docked;
		dock_data = &data_rep->docked;
		if (dock_ptr->G_binding == DOUBLE_INIT) dock_ptr->G_binding = dock_data->G_binding;
		if (dock_ptr->Ki_DOCK == DOUBLE_INIT) dock_ptr->Ki_DOCK = dock_data->Ki_DOCK;
		if (strcmp(dock_ptr->Ki_unit,STR_INIT) == 0) {
			newlen = strlen(dock_data->Ki_unit) + 1;
			oldlen = strlen(dock_ptr->Ki_unit) + 1;
			memset(dock_ptr->Ki_unit,0,oldlen);
			dock_ptr->Ki_unit = realloc(dock_ptr->Ki_unit, newlen*sizeof(char));
			strncpy(dock_ptr->Ki_unit,dock_data->Ki_unit,newlen);
		}
		if (dock_ptr->E_inter == DOUBLE_INIT) dock_ptr->E_inter = dock_data->E_inter;
		if (dock_ptr->E_nonbond == DOUBLE_INIT) dock_ptr->E_nonbond = dock_data->E_nonbond;
		if (dock_ptr->E_electrostat == DOUBLE_INIT)
			dock_ptr->E_electrostat = dock_data->E_electrostat;
		if (dock_ptr->E_internal == DOUBLE_INIT) dock_ptr->E_internal = dock_data->E_internal;
		if (dock_ptr->G_tors == DOUBLE_INIT) dock_ptr->G_tors = dock_data->G_tors;
		if (dock_ptr->E_unbound == DOUBLE_INIT) dock_ptr->E_unbound = dock_data->E_unbound;
		if (dock_ptr->rmsd_ref == DOUBLE_INIT) dock_ptr->rmsd_ref = dock_data->rmsd_ref;
		if (dock_ptr->time == DOUBLE_INIT) dock_ptr->time = dock_data->time;

		// Update the quantum mechanical results
		struct QMresult *qm_ptr, *qm_data;
		qm_ptr = &curr_rep->optimized;
		qm_data = &data_rep->optimized;
		if (qm_ptr->Hf == DOUBLE_INIT) qm_ptr->Hf = qm_data->Hf;
		if (qm_ptr->S == DOUBLE_INIT) qm_ptr->S = qm_data->S;
		if (qm_ptr->Cp == DOUBLE_INIT) qm_ptr->Cp = qm_data->Cp;
		if (qm_ptr->G == DOUBLE_INIT) qm_ptr->G = qm_data->G;
		if (qm_ptr->ZPE == DOUBLE_INIT) qm_ptr->ZPE = qm_data->ZPE;
		if (qm_ptr->E_dielec == DOUBLE_INIT) qm_ptr->E_dielec = qm_data->E_dielec;
		if (qm_ptr->mu_x == DOUBLE_INIT) qm_ptr->mu_x = qm_data->mu_x;
		if (qm_ptr->mu_y == DOUBLE_INIT) qm_ptr->mu_y = qm_data->mu_y;
		if (qm_ptr->mu_z == DOUBLE_INIT) qm_ptr->mu_z = qm_data->mu_z;
		if (qm_ptr->mu_total == DOUBLE_INIT) qm_ptr->mu_total = qm_data->mu_total;
		if (qm_ptr->time == DOUBLE_INIT) qm_ptr->time = qm_data->time;
		if (strcmp(qm_ptr->method,STR_INIT) == 0) {
			newlen = strlen(qm_data->method) + 1;
			oldlen = strlen(qm_ptr->method) + 1;
			memset(qm_ptr->method, 0,oldlen);
			qm_ptr->method = realloc(qm_ptr->method, newlen*sizeof(char));
			strncpy(qm_ptr->method,qm_data->method,newlen);
		}
		if (qm_ptr->num_SCFs == INT_INIT) qm_ptr->num_SCFs = qm_data->num_SCFs;
		if (qm_ptr->COSMO_A == DOUBLE_INIT) qm_ptr->COSMO_A = qm_data->COSMO_A;
		if (qm_ptr->COSMO_V == DOUBLE_INIT) qm_ptr->COSMO_V = qm_data->COSMO_V;
		if (qm_ptr->vdW_A == DOUBLE_INIT) qm_ptr->vdW_A = qm_data->vdW_A;
		if (qm_ptr->G_prot == DOUBLE_INIT) qm_ptr->G_prot = qm_data->G_prot;
		//if (qm_ptr->G_liginput == DOUBLE_INIT) qm_ptr->G_liginput = qm_data->G_liginput;
		if (qm_ptr->G_ligand == DOUBLE_INIT) qm_ptr->G_ligand = qm_data->G_ligand;
		if (qm_ptr->G_binding == DOUBLE_INIT) qm_ptr->G_binding = qm_data->G_binding;
		if (qm_ptr->Ki_QM == DOUBLE_INIT) qm_ptr->Ki_QM = qm_data->Ki_QM;
		if (strcmp(qm_ptr->Ki_type,STR_INIT) == 0) {
			newlen = strlen(qm_data->Ki_type) + 1;
			oldlen = strlen(qm_ptr->Ki_type) + 1;
			memset(qm_ptr->Ki_type, 0,oldlen);
			qm_ptr->Ki_type = realloc(qm_ptr->Ki_type, newlen*sizeof(char));
			strncpy(qm_ptr->Ki_type,qm_data->Ki_type,newlen);
		}
		data_rep = data_rep->next;
	}

	fflush(stdout);

        CompoundNode->data = malloc(sizeof(CompoundTree));

	return retval;

}


job_t* InitJob(job_t *newjob)
{
	newjob = (job_t*)malloc(sizeof(job_t));
	newjob->name = (char*)calloc(FILENAME_MAX,sizeof(char));
	newjob->mopac_res = (char*)calloc(FILENAME_MAX,sizeof(char));
	SetupSTIC(&newjob->input_data);
	newjob->host_rank = -1;
	newjob->type = DOCK;
	return newjob;
}

void FreeJob(job_t *job)
{
	free(job->name);
	free(job->mopac_res);
	free(job);
}

/**
 * Creates a node for the new job and inserts it into the job
 * tree. If the tree does not exist yet, it is created with
 * newjob as its only member. Otherwise, newjob is inserted.
 *
 * Returns: pointer to the job node created for newjob
 */
deque_node* InsertJob(deque *jobs, char newjob[FILENAME_MAX], JobParameters *params, int energy_count)
{
	deque_node *jobnode = NULL;
	job_t *jobinfo = NULL;
	char *token, *saveptr;
	float temp;

	//Create new job node
	jobinfo = InitJob(jobinfo);
	jobinfo->type = (params && params->prescreen) ? PRESCREENING : DOCK;

	jobnode = InitDequeNode(jobnode);
	jobnode->data = jobinfo;

	// If enthalpy, and optionally free energy, are provided,
	// they should be the second and third tokens, respectively,
	// with a space delimiter.
	token = strtok_r(newjob," \n",&saveptr);
	if(token != NULL)
	{
		strcpy(jobinfo->name,token);
		token = strtok_r(NULL," \n",&saveptr);

		//Enthalpy data?
		if(token != NULL)
		{
			if(energy_count == 0)
				fprintf(stderr,"Warning - InsertJob did not expect energy data.");
			if(sscanf(token,"%f",&temp) == 1)
			{
				jobinfo->input_data.Hf = temp;
				if(params->input_units != INPUT_KCAL)
					ConvertToKCals(&jobinfo->input_data.Hf,params->input_units);
			}
			else
				fprintf(stderr, "Warning - InsertJob could not read expected floating point data for Enthalpy: %s",token);

			token = strtok_r(NULL," \n",&saveptr);
		}
		else if(energy_count > 0)
			fprintf(stderr,"Warning - InsertJob expected Enthalpy data.\n");

		//Free energy data?
		if(token != NULL)
		{
			if(energy_count != 2)
				fprintf(stderr,"Warning - InsertJob did not expect free energy data.");
			if(sscanf(token,"%f",&temp) == 1)
			{
				jobinfo->input_data.G = temp;
				if(params->input_units != INPUT_KCAL)
					ConvertToKCals(&jobinfo->input_data.G,params->input_units);
			}
			else
				fprintf(stderr, "Warning - InsertJob could not read expected floating point data for Free Energy: %s",token);

		}
		else if(energy_count > 1)
			fprintf(stderr,"Warning - InsertJob expected Free Energy data.\n");

	}

	//Insert new job node in tree
	if(jobs != NULL)
		push_back(jobs,jobnode);
	return jobnode;
}

/**
 * Creates (and returns) a tree containing jobs for
 * each entry in the specified input file. Lines
 * in the input file may be commented out by
 * beginning the line with a '#' character.
 */
deque* CreateJobList(char inpfile[FILENAME_MAX], JobParameters *params)
{
	//
	// THIS IS ONE SUPER LONG FUNCTION. PROBABLY SHOULD BE FACTORED...
	//
    char nextjob[FILENAME_MAX];
    char *token, *saveptr;// for string tokenizing
    FILE *jobfile;
    deque *jobs = NULL;
    int energy_count = 0;//energy_count = number of energy data expected on each line.
    int length = 0;

    jobs = InitDeque(jobs);

    memset(nextjob, 0,FILENAME_MAX);

    printf_master("Opening ligand file: %s ... ",inpfile);
    jobfile = fopen(inpfile,"r");

    if (jobfile == NULL) {
    	printf_master("FAILED!!\n");
        printf_master("ERROR - %s\n",strerror(errno));
        MPI_Abort(MPI_COMM_WORLD,DEFAULT_ERR_CODE);
    }
    else {
	printf_master ("SUCCESS!\n");
    }

    // Get first line. If it contains optional key words, parse them.
    // Otherwise move on.
    //
    // NOTE: This needs to be revamped to provide a more free-form input
    //       Should just look for keywords on each line and then start
    //       reading ligands when a new keywork (e.g lig_list) is found
    if(fgets(nextjob,FILENAME_MAX,jobfile)!=NULL)
    {
    	int get_next_line = TRUE;
    	char tmpstring[FILENAME_MAX];

    	length = strlen(nextjob);
    	if(length > 1)
    	{
    		strcpy(tmpstring,nextjob);
    		token = strtok_r(nextjob," \n",&saveptr);
    		while(token != NULL)
    		{
    			if(strcmp(token,"enthalpy") == 0)
    			{
    				energy_count++;
    			}
    			else if(token[0] == '#')
    				break;
    			else if(strcmp(token,"free_energy") == 0)
    			{
    				energy_count = 2;
    				params->UseFreeEnergy = TRUE;
    			}
    			else if(strcmp(token,"kJ") == 0)
    				params->input_units = INPUT_KJ;
    			else if(strcmp(token,"kcal") == 0)
    				params->input_units = INPUT_KCAL;
    			else if(strcmp(token,"cal") == 0)
    				params->input_units = INPUT_CAL;
    			else if(strcmp(token,"J") == 0)
    				params->input_units = INPUT_J;
    			else
    			{
    				get_next_line = FALSE;
    				strcpy(nextjob,tmpstring);
    				break;
    			}
    			token = strtok_r(NULL," \n",&saveptr);
    		}
    		if( get_next_line && fgets(nextjob,FILENAME_MAX,jobfile) == NULL)
    		{
    			fprintf(stderr, "ERROR - Broken ligand input file (%s). Expecting target energy line.\n",inpfile);
    			MPI_Abort(MPI_COMM_WORLD,DEFAULT_ERR_CODE);
    		}
    	}
    }

    // At this point nextjob has what should be the target_energy line.
    if(nextjob == NULL || (token = strtok_r(nextjob," \n",&saveptr)) == NULL || strcmp(token,"target_energy") != 0)
    {
    	printf("ERROR - Missing target_energy line. The first required line must be of the form:\n");
    	printf("ERROR - \n");
    	printf("ERROR - target_energy <enthalpy> [free_energy]\n");
    	printf("ERROR - \n");
    	printf("ERROR - Where free energy is an optional energy value.\n");
    	MPI_Abort(MPI_COMM_WORLD,DEFAULT_ERR_CODE);
    }
    else
    {
    	float temp_input;
    	token = strtok_r(NULL," \n",&saveptr);

    	if (token == NULL)
    	{
    		printf("ERROR - Missing target enthalpy from input file: %s\n",inpfile);
    		MPI_Abort(MPI_COMM_WORLD,DEFAULT_ERR_CODE);
    	}

    	if(sscanf(token,"%f",&temp_input) == 1)
    	{
    		params->target_Hf = temp_input;
    		if(params->input_units != INPUT_KCAL)
    			ConvertToKCals(&params->target_Hf,params->input_units);
    	}
    	else
    		fprintf(stderr, "Warning - Could not read Target Enthalpy: %s",token);

    	if(energy_count == 2)
    	{
    		token = strtok_r(NULL," \n",&saveptr);
    		if (token == NULL)
    		{
    			printf("ERROR - Missing target free energy from input file: %s\n",inpfile);
    			MPI_Abort(MPI_COMM_WORLD,DEFAULT_ERR_CODE);
    		}

    		if(sscanf(token,"%f",&temp_input) == 1)
    		{
    			params->target_G = temp_input;
    			if(params->input_units != INPUT_KCAL)
    				ConvertToKCals(&params->target_G,params->input_units);
    		}
    		else
    			fprintf(stderr, "Warning - Could not read Target Enthalpy: %s",token);

    	}
    }

    printf_master("Creating job list from file ... ");
    while (fgets(nextjob,FILENAME_MAX,jobfile)!=NULL) {
	length = strlen(nextjob);
	if (length > 1) {
            if (nextjob[0] != '#') {
	        nextjob[length-1] = '\0';
                InsertJob(jobs,nextjob, params, energy_count);
            }
	}
    }
    printf_master("done!\n");

    if (fclose(jobfile) != 0) {
	printf_master("ERROR: Failed to close file!!!\n");
        printf_master("ERROR: %s\n",strerror(errno));
	MPI_Abort(MPI_COMM_WORLD,DEFAULT_ERR_CODE);
    }

    return jobs;
}


