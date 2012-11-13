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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "defines.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "structs.h"
#include "RBTree.h"



void ResortRBTree(RBTree **dest,const RBTree *source,int sort_criterion)
{
	RBTree *new_node = NULL;
	if(source == NULL)
		return;

	new_node = InitRBTree(new_node,source->type);
	if(source->type == STIC_AVG)
	{
		struct STICAvg *new_stic = (struct STICAvg*)malloc(sizeof(struct STICAvg));
		*new_stic = *((struct STICAvg*)source->data);
		new_stic->Ki.sort_criterion = sort_criterion;
		new_node->data = new_stic;
	}
	else if(source->type == COMPOUND_AVG)
	{
		struct CompoundAvg *new_cmpd = (struct CompoundAvg*)malloc(sizeof(struct CompoundAvg));
		*new_cmpd = *((struct CompoundAvg*)source->data);
		new_cmpd->Ki.sort_criterion = sort_criterion;
		new_node->data = new_cmpd;
	}

	insertRBT(dest,new_node);
	ResortRBTree(dest,source->left,sort_criterion);
	ResortRBTree(dest,source->right,sort_criterion);
}


/*
 * End result is in Molars. Conversion from other units is done within
 * this function.
 */
void BoltzmannAvgSTIC(struct STICelement *stic, double *Z_list_DOCK, double *Z_list_QM, double *Z_DOCK, double *Z_QM, double *dock,double *qm)
{
	struct ClusterRep *curr_rep;
	double factor;
	size_t part_counter = 0;

	if(stic == NULL || dock == NULL || qm == NULL || Z_DOCK == NULL || Z_QM == NULL)
	{
		fprintf(stderr,"BoltzmannAvgSTIC: NULL POINTER\n");
		exit(NULL_POINTER);
	}

	// Initialize partition coefficient sums
	*Z_QM = 0;
	*Z_DOCK = 0;
	*dock = *qm = 0;

	//Is there anything to average?
	if(stic->reps == NULL)
		return;

	//Iterate through reps, sum Ki * boltzmann factor and sum boltzmann factors
	curr_rep = stic->reps;
	while(curr_rep != NULL)
	{
		// Dock
		// Autodock gives the actual G and Ki for products.
		// Therefore, Boltzmann average Ki and weight based on G
		if(curr_rep->docked.G_binding != DOUBLE_INIT && curr_rep->docked.Ki_DOCK != DOUBLE_INIT)
		{
		  double temp_ki = curr_rep->docked.Ki_DOCK;
			factor = curr_rep->docked.G_binding;
			factor = exp(factor/EFFERCIO_KT);
			//put Ki unit conversion (to Molars) in factor;
			switch(curr_rep->docked.Ki_unit[0])
			{
			case 'u':
				temp_ki *= 1e-6;
				break;
			case 'n':
				temp_ki *= 1e-9;
				break;
			case 'm':
				temp_ki *= 1e-3;
				break;
			case 'p':
				temp_ki *= 1e-12;
				break;
			default:
				fprintf(stderr,"BoltzmannAvgSTIC: Unknown units for Ki: %s\n",curr_rep->docked.Ki_unit);
			case 'M':
				break;
			}
			Z_list_DOCK[part_counter] = factor;
			*Z_DOCK += factor;
			*dock += factor*temp_ki;

		}

		// QM calculated Ki based on G by the formula Ki := exp(-(G_prod - G_react)/kT)
		// Therefore the Boltzman average of Ki, weighted by G,
		// is simplified using Ki in place of the exponential factor of the average
		if(curr_rep->optimized.Hf != DOUBLE_INIT)
		{
			double E_ligand, E_prod;
			if(curr_rep->optimized.G != DOUBLE_INIT)
			{
				E_ligand = stic->G;
				E_prod = curr_rep->optimized.G;
			}
			else
			{
				E_ligand = stic->Hf;
				E_prod = curr_rep->optimized.Hf;
			}

			if(E_prod != DOUBLE_INIT && E_ligand != DOUBLE_INIT && curr_rep->optimized.G_prot != DOUBLE_INIT)
			{
				factor = E_prod - E_ligand - curr_rep->optimized.G_prot;
				factor = exp(factor/EFFERCIO_KT);
				curr_rep->optimized.Ki_QM = factor;//Here, units are in M, while Autodock will probably use mM, uM or even nM

				Z_list_QM[part_counter] = curr_rep->optimized.Ki_QM;
				*Z_QM += curr_rep->optimized.Ki_QM;
				*qm += curr_rep->optimized.Ki_QM*curr_rep->optimized.Ki_QM;
			}
		}

		part_counter++;
		curr_rep = curr_rep->next;
	}

	// Divide data by sums
	if(*Z_DOCK == 0)
		*dock = DOUBLE_INIT;
	else
		*dock /= *Z_DOCK;

	if(*Z_QM == 0)
		*qm = DOUBLE_INIT;
	else
		*qm /= *Z_QM;

	part_counter = 0;
	curr_rep = stic->reps;
	while(curr_rep != 0)
	{
		if(*Z_DOCK != 0)
			Z_list_DOCK[part_counter] /= *Z_DOCK;
		else
			Z_list_DOCK[part_counter]  = 0;

		if(*Z_QM != 0)
			Z_list_QM[part_counter] /= *Z_QM;
		else
			Z_list_QM[part_counter] = 0;
		part_counter++;
		curr_rep = curr_rep->next;
	}

	return;

}

void BoltzmannAvgSTICTree(RBTree **averages, double *Z, const RBTree *stics)
{
	struct STICelement *curr_stic = NULL;
	struct STICAvg *new_avg = NULL;
	RBTree *new_node = NULL;
	double Z_CR_DOCK, Z_CR_QM, E;
	double *Z_list_DOCK, *Z_list_QM;
	size_t num_reps = 0;

	if(stics == NULL || averages == NULL)
		return;


	curr_stic = (struct STICelement *)stics->data;
	num_reps = NumClusterReps(curr_stic);
	Z_list_DOCK = (double*)malloc(num_reps*sizeof(double));
	Z_list_QM = (double*)malloc(num_reps*sizeof(double));

	new_node = InitRBTree(new_node,STIC_AVG);
	new_avg = (struct STICAvg*)malloc(sizeof(struct STICAvg));
	new_avg->stic = curr_stic;
	new_avg->Ki.sort_criterion = AVG_SORT_DOCK;
	new_node->data = new_avg;

	// Average cluster reps to get the STIC's Ki
	BoltzmannAvgSTIC(curr_stic,Z_list_DOCK,Z_list_QM,&Z_CR_DOCK,&Z_CR_QM,&new_avg->Ki.Ki_DOCK, &new_avg->Ki.Ki_QM);
	free(Z_list_DOCK);
	free(Z_list_QM);

	// Initialize the STIC's partition coefficient. This will be divided by Z_DOCK and Z_QM later
	if(curr_stic->G != DOUBLE_INIT)
		E = curr_stic->G;
	else
		E = curr_stic->Hf;
	if(E != DOUBLE_INIT)
	{
		new_avg->Z = exp(E/EFFERCIO_KT);
		*Z += new_avg->Z;
	}

	insertRBT(averages,new_node);

	BoltzmannAvgSTICTree(averages,Z,stics->left);
	BoltzmannAvgSTICTree(averages,Z,stics->right);
}

void GetCompoundSum(double *Ki_DOCK, double *Ki_QM,
		double *Z_DOCK, double *Z_QM, RBTree *avgs, int give_predicted)
{
	struct STICAvg *curr_node = NULL;
	if(avgs == NULL)// end of recursion
		return;

	curr_node = (struct STICAvg*)avgs->data;

	if(give_predicted || curr_node->Ki.Ki_DOCK != DOUBLE_INIT)
	{
		*Z_DOCK += curr_node->Z;
		*Ki_DOCK += curr_node->Z * curr_node->Ki.Ki_DOCK;
	}

	if(give_predicted || curr_node->Ki.Ki_QM != DOUBLE_INIT)
	{
		*Z_QM += curr_node->Z;
		*Ki_QM += curr_node->Z * curr_node->Ki.Ki_QM;
	}

	GetCompoundSum(Ki_DOCK, Ki_QM, Z_DOCK, Z_QM, avgs->left, give_predicted);
	GetCompoundSum(Ki_DOCK, Ki_QM, Z_DOCK, Z_QM, avgs->right, give_predicted);
}

void DivideSTICCoefficients(double Z_DOCK,RBTree *avg)
{
	struct STICAvg *data;
	if(avg == NULL)
		return;
	data = (struct STICAvg*)avg->data;
	if(data == NULL)
		return;

	data->Z /= Z_DOCK;

	DivideSTICCoefficients(Z_DOCK,avg->left);
	DivideSTICCoefficients(Z_DOCK,avg->right);
}

void GetCompoundAverage(double *Ki_DOCK, double *Ki_QM, RBTree *avg, int give_predicted)
{
	double Z_DOCK, Z_QM;
	Z_DOCK = Z_QM = 0;
	*Ki_DOCK = 0;*Ki_QM = 0;
	GetCompoundSum(Ki_DOCK,Ki_QM,&Z_DOCK,&Z_QM,avg, give_predicted);
	*Ki_DOCK /= Z_DOCK;
	*Ki_QM /= Z_QM;
	DivideSTICCoefficients(Z_DOCK,avg);
}

/**
 * Averages over stics to obtain Ki data for a particular compound.
 * The tree "average" will contain Ki values for QM and DOCK results,
 * but will be sorted by DOCK values.
 */
void BoltzmannAvgCompounds(RBTree **average, const RBTree *CompoundList)
{
	double Z_stics;// sum of partition functions for stics
	struct CompoundAvg* avg = NULL;
	RBTree *new_node = NULL,*stic_averages = NULL;

	if(CompoundList == NULL || CompoundList->data == NULL)// end of recursion
		return;

	// Generate coefficients
	Z_stics = 0;
	BoltzmannAvgSTICTree(&stic_averages,&Z_stics,((CompoundTree*)CompoundList->data)->stics);

	// Initialize Tree node for this compound
	new_node = InitRBTree(new_node,COMPOUND_AVG);
	avg = (struct CompoundAvg*)malloc(sizeof(struct CompoundAvg));
	avg->compound = (CompoundTree*)CompoundList->data;
	avg->Ki.Ki_DOCK = avg->Ki.Ki_QM = 0.0;
	avg->Ki.sort_criterion = AVG_SORT_DOCK;
	new_node->data = avg;

	// Average coefficients
	GetCompoundAverage(&avg->Ki.Ki_DOCK,&avg->Ki.Ki_QM,stic_averages,FALSE);
	FreeRBTree(stic_averages);

	// Insert data into tree
	insertRBT(average,new_node);

	// Recursively work down tree
	BoltzmannAvgCompounds(average,CompoundList->left);
	BoltzmannAvgCompounds(average,CompoundList->right);

}

void AdjustUnits(double *val, char *units)
{
  double pow_10;
  if(val == NULL)
    return;

  if(*val <= 0)
  {
	  units[0] = (char)0;
	  return;
  }

  pow_10 = log10(*val);
  if(pow_10 > 0)
    {
      strcpy(units,"");
      return;
    }
  else if(pow_10 > -4)
    {
      *val *= 1e3;
      strcpy(units,"m");
      return;
    }
  else if(pow_10 > -7)
    {
      *val *= 1e6;
      strcpy(units,"u");
      return;
    }
  else if(pow_10 > -10)
    {
      *val *= 1e9;
      strcpy(units,"n");
      return;
    }
  else
    {
      *val *= 1e12;
      strcpy(units,"p");
    }
  return;
}

void FWriteClusterRep(FILE *nmr_file,const RBTree *cr_avgs)
{
	struct ClusterRep *rep = NULL;
	struct ClusterRepAvg *avg = NULL;
	char units[16];
	double adjusted_value;
	if(cr_avgs == NULL)// end of recursion
		return;

	FWriteClusterRep(nmr_file,cr_avgs->left);

	if(nmr_file == NULL)
		nmr_file = stdout;

	avg = (struct ClusterRepAvg*)cr_avgs->data;
	rep = avg->rep;

	fprintf(nmr_file,"          Cluster Rep %d  ",rep->index);
	if(avg->Ki.sort_criterion == AVG_SORT_DOCK)
	  adjusted_value = avg->Ki.Ki_DOCK;
	else
	  adjusted_value = avg->Ki.Ki_QM;

	AdjustUnits(&adjusted_value,units);

	if(adjusted_value == DOUBLE_INIT)
	{
		fprintf(stdout,"WARNING - Uninitialized Cluster representative %d\n",rep->index);
		fprintf(nmr_file," FAILED\n");
	}
	else
		fprintf(nmr_file,"%.03f %sM (%.03f)\n",adjusted_value,units,avg->Z);


	FWriteClusterRep(nmr_file,cr_avgs->right);
}

void FWriteSTICAvg_rec(FILE *nmr_file,const RBTree *stic, const char *cmpd_id)
{
	struct STICelement *data =NULL;
	struct STICAvg *avg = NULL;
	double fval = 0, *Z_list_DOCK, *Z_list_QM, Z_QM, Z_DOCK;
	RBTree *cr_avgs = NULL;
	size_t num_reps = 0;
	char units[16];
	int compound_failed = FALSE;
	if(stic == NULL)// end of recursion
		return;

	FWriteSTICAvg_rec(nmr_file,stic->left,cmpd_id);
	if(nmr_file == NULL)
		nmr_file = stdout;

	avg = (struct STICAvg*)stic->data;
	data = avg->stic;
	fval = (avg->Ki.sort_criterion == AVG_SORT_DOCK) ? avg->Ki.Ki_DOCK : avg->Ki.Ki_QM;

	// Check to see if we should report data or FAILED
	if(avg->Ki.sort_criterion == AVG_SORT_DOCK && avg->Ki.Ki_DOCK == DOUBLE_INIT)
		compound_failed = TRUE;
	if(avg->Ki.sort_criterion == AVG_SORT_QM && avg->Ki.Ki_QM == DOUBLE_INIT)
		compound_failed = TRUE;

	AdjustUnits(&fval,units);
	if(!compound_failed)
	{
		fprintf(nmr_file,"     %s_S%03dT%03dI%03dC%03d  %.03f %sM (%.03f)\n",
				cmpd_id,data->S,data->T,data->I,data->C,fval,units,avg->Z);
	}
	else
	{
		fprintf(stdout,"WARNING - STIC averaging failed for  %s_S%03dT%03dI%03dC%03d\n",cmpd_id,data->S,data->T,data->I,data->C);
		fprintf(nmr_file,"     %s_S%03dT%03dI%03dC%03d  FAILED (expected %.03f)\n",cmpd_id,data->S,data->T,data->I,data->C,avg->Z);
	}
	fprintf(nmr_file,"     %s\n",WIGL_HDR);

	num_reps = NumClusterReps(data);
	if(!compound_failed && num_reps != 0)
	{
		double ki_dock,ki_qm;
		int i = 0;
		struct ClusterRep *curr_rep = data->reps;

		Z_list_DOCK = (double*)malloc(num_reps*sizeof(double));
		Z_list_QM = (double*)malloc(num_reps*sizeof(double));
		BoltzmannAvgSTIC(data,Z_list_DOCK,Z_list_QM,&Z_DOCK,&Z_QM,&ki_dock,&ki_qm);
		for(;i<num_reps;i++)
		{
			struct ClusterRepAvg *cra = (struct ClusterRepAvg*)malloc(sizeof(struct ClusterRepAvg));
			RBTree *new_node;
			new_node = InitRBTree(new_node,CR_AVG);

			cra->Ki.Ki_DOCK = curr_rep->docked.Ki_DOCK;
			switch(curr_rep->docked.Ki_unit[0])
			{
			case 'u':
				cra->Ki.Ki_DOCK *= 1e-6;
				break;
			case 'n':
				cra->Ki.Ki_DOCK *= 1e-9;
				break;
			case 'm':
				cra->Ki.Ki_DOCK *= 1e-3;
				break;
			case 'p':
				cra->Ki.Ki_DOCK *= 1e-12;
				break;
			default:
				fprintf(stderr,"Unknown units for Ki: %s\n",curr_rep->docked.Ki_unit);
			case 'M':
				break;
			}
			cra->Ki.Ki_QM = curr_rep->optimized.Ki_QM;
			cra->Ki.sort_criterion = avg->Ki.sort_criterion;
			if(cra->Ki.sort_criterion == AVG_SORT_DOCK)
				cra->Z = Z_list_DOCK[i];
			else
				cra->Z = Z_list_QM[i];

			cra->rep = curr_rep;

			new_node->data = cra;
			insertRBT(&cr_avgs,new_node);

			curr_rep = curr_rep->next;
		}
		free(Z_list_DOCK);
		free(Z_list_QM);

		if(cr_avgs != NULL)
		{
			FWriteClusterRep(nmr_file,cr_avgs);
			FreeRBTree(cr_avgs);
		}
	}
	fprintf(nmr_file,"\n");

	FWriteSTICAvg_rec(nmr_file,stic->right,cmpd_id);
}

void FWriteSTICAvg(FILE *nmr_file,CompoundTree *compound, int sort_criterion)
{
	double Z_stics,Ki_DOCK,Ki_QM;
	RBTree *stic_averages = NULL, *stic_qm_averages = NULL;
	if(nmr_file == NULL)
		nmr_file = stdout;

	// Generate coefficients
	Z_stics = 0;
	BoltzmannAvgSTICTree(&stic_averages,&Z_stics,compound->stics);
	if(Z_stics == 0)
	{
		Ki_DOCK = Ki_QM = 0;
	}
	else
		GetCompoundAverage(&Ki_DOCK,&Ki_QM,stic_averages,TRUE);


	if(sort_criterion == AVG_SORT_QM)
	{
		ResortRBTree(&stic_qm_averages,stic_averages,AVG_SORT_QM);
		FWriteSTICAvg_rec(nmr_file,stic_qm_averages,compound->data->ID);
	}
	else
		FWriteSTICAvg_rec(nmr_file,stic_averages,compound->data->ID);

	FreeRBTree(stic_averages);
	FreeRBTree(stic_qm_averages);
}


void FWriteCompoundAvg(FILE *summary_file,FILE *nmr_file,FILE **well_file,
		int *compound_counter,int max_compounds,
		const RBTree *averages,const char *prefix,int data_type)
{
	char cmpd_format[16];
	char units[16];
	const struct CompoundAvg *cmpd = NULL;
	double data = 0;
	int compound_failed = FALSE;

	if(averages == NULL)// end of recursion
		return;
	FWriteCompoundAvg(summary_file,nmr_file,well_file,compound_counter,max_compounds,averages->left,prefix,data_type);

	if(summary_file == NULL)
	{
		summary_file = stdout;
	}
	if(nmr_file == NULL)
	{
		// if both are null, only write nmr_file, as it is the most verbose.
		if(summary_file == stdout)
			summary_file = NULL;
		nmr_file = stdout;
	}

	*compound_counter += 1;
	if(*compound_counter % max_compounds == 0)
	{
		if(well_file != NULL && *well_file != NULL && *well_file != stdout)
		{
			fclose(*well_file);
			*well_file = NULL;
		}
	}

	cmpd = (struct CompoundAvg*)averages->data;
	if(data_type == AVG_SORT_DOCK)
            data = cmpd->Ki.Ki_DOCK;
	else
            data = cmpd->Ki.Ki_QM;

	compound_failed = (data == DOUBLE_INIT) || isnan(data);

	if(compound_failed)
	{
		fprintf(stdout,"WARNING - Averaging of compound %s FAILED.\n",cmpd->compound->data->ID);
		strcpy(cmpd_format,"%s FAILED\n");
	}
	else
	{
		strcpy(cmpd_format,"%s %.03f %sM\n");
	}

	AdjustUnits(&data,units);

	// If well file is null, try to open a file for it. If that fails, move on; the whole list will be eith
	if(well_file != NULL && *well_file == NULL)
	{
		char well_filename[FILENAME_MAX];
		sprintf(well_filename,"%s_%03d",prefix,((int) *compound_counter/max_compounds));
		*well_file = fopen(well_filename,"w");
	}
	if(well_file != NULL && *well_file != NULL)
	{
		if(!compound_failed)
			fprintf(*well_file,cmpd_format,cmpd->compound->data->ID,data,units);
		else
			fprintf(*well_file,cmpd_format,cmpd->compound->data->ID);
	}
	if(summary_file != NULL)
	{
		if(!compound_failed)
			fprintf(summary_file,cmpd_format,cmpd->compound->data->ID,data,units);
		else
			fprintf(summary_file,cmpd_format,cmpd->compound->data->ID);
	}
	// Print calls.
	if(!compound_failed)
		fprintf(nmr_file,cmpd_format,cmpd->compound->data->ID,data,units);
	else
		fprintf(nmr_file,cmpd_format,cmpd->compound->data->ID);

	fprintf(nmr_file,"%s\n",EQ_HDR);
	FWriteSTICAvg(nmr_file,cmpd->compound,cmpd->Ki.sort_criterion);
	fprintf(nmr_file,"\n");

	FWriteCompoundAvg(summary_file,nmr_file,well_file,compound_counter,max_compounds,averages->right,prefix,data_type);
}

void BoltzmannAvgCompoundTree(RBTree *CompoundList, const char *analysis_dir, int use_gibbs)
{

	RBTree *averages = NULL, *avg_QM = NULL;
	FILE *summary_file = NULL, *well_file = NULL, *nmr_file = NULL;
	char summary_filename[FILENAME_MAX], nmr_filename[FILENAME_MAX], *ID = NULL;
	char adk_dir[FILENAME_MAX],qm_dir[FILENAME_MAX];
	int compound_counter;
	static const char *master_name = "MASTER";
	static const char *nmr_name = "NMR_DATA";

	if(VerifyDir(analysis_dir,1,"MASTER") != 0)
	{
		fprintf(stderr,"ERROR - Could not create analysis dir. Using current directory instead.\n");
		analysis_dir = "./";
	}

	sprintf(adk_dir,"%s/adk",analysis_dir);
	sprintf(qm_dir,"%s/qm",analysis_dir);
	VerifyDir(adk_dir,1,"MASTER");
	VerifyDir(qm_dir,1,"MASTER");

	// Obtain a tree of CompoundAvg structs with Ki data, ordered by
	// Ki_DOCK.
	BoltzmannAvgCompounds(&averages,CompoundList);

	// Write the average to file.
	ID = ((CompoundTree*)CompoundList->data)->data->ID;
	sprintf(summary_filename,"%s/%s_DOCK",adk_dir,master_name);
	summary_file = fopen(summary_filename,"w");
	if(summary_file == NULL)
	{
		fprintf(stderr,"Could not open %s for writing\n Writing to stdout\n",summary_filename);
		fprintf(stderr,"Reason: %s\n",strerror(errno));
		summary_file = stdout;
	}

	sprintf(nmr_filename,"%s/%s_DOCK",adk_dir,nmr_name);
	nmr_file = fopen(nmr_filename,"w");
	if(nmr_file == NULL)
	{
		fprintf(stderr,"Could not open %s for writing\n Writing to stdout\n",summary_filename);
		fprintf(stderr,"Reason: %s\n",strerror(errno));
		nmr_file = stdout;
		if(summary_file == stdout)
			summary_file = NULL;
	}

	compound_counter = 0;
	FWriteCompoundAvg(summary_file,nmr_file,&well_file,&compound_counter,96,averages,summary_filename,AVG_SORT_DOCK);

	if(summary_file != NULL && summary_file != stdout)
		fclose(summary_file);
	if(nmr_file != NULL && nmr_file != stdout)
		fclose(nmr_file);
	if(well_file != NULL && well_file != stdout)
	{
		fclose(well_file);
		well_file = NULL;
	}

	// Resort for QM data
	ResortRBTree(&avg_QM,averages,AVG_SORT_QM);
	FreeRBTree(averages);

	// Write QM data, just like dock.
	sprintf(summary_filename,"%s/%s_QM",qm_dir,master_name);
	summary_file = fopen(summary_filename,"w");
	if(summary_file == NULL)
	{
		fprintf(stderr,"Could not open %s for writing\n Writing to stdout\n",summary_filename);
		fprintf(stderr,"Reason: %s\n",strerror(errno));
		summary_file = stdout;
	}

	sprintf(nmr_filename,"%s/%s_QM",qm_dir,nmr_name);
	nmr_file = fopen(nmr_filename,"w");
	if(nmr_file == NULL)
	{
		fprintf(stderr,"Could not open %s for writing\n Writing to stdout\n",summary_filename);
		fprintf(stderr,"Reason: %s\n",strerror(errno));
		nmr_file = stdout;
		if(summary_file == stdout)
			summary_file = NULL;
	}

	compound_counter = 0;
	FWriteCompoundAvg(summary_file,nmr_file,&well_file,&compound_counter,96,avg_QM,summary_filename,AVG_SORT_QM);

	if(summary_file != NULL && summary_file != stdout)
		fclose(summary_file);
	if(well_file != NULL && well_file != stdout)
		fclose(well_file);

	FreeRBTree(avg_QM);
}

