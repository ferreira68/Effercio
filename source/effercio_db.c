// $Id: effercio_db.c 6 2011-11-17 19:55:32Z dcoss $


#include "effercio_db.h"

static void FPrintSQLClusterRep(FILE *sql_file, int compound_id, struct ClusterRep *rep, struct STICelement *stic, JobParameters *params)
{
  struct QMresult *qm = NULL;
  struct DOCKresult *dock = NULL;

  if(rep == NULL || stic == NULL)
    return;
  if(params == NULL || params->qm_method == NULL)
    {
      fprintf(stderr,"ERROR - No QM method supplied to FPrintSQLClusterRep\n");
      return;
    }
  
  if(sql_file == NULL)
    sql_file = stdout;

  qm = &rep->optimized;
  dock = &rep->docked;

  // Upsert data
  fprintf(sql_file, "insert into charge_models(stic_id, method_id, cluster_rep_idx, g_binding  , ki_dock    , ki_units, e_inter, e_nonbond, e_electrostat, e_internal, g_tors, e_unbound, rmsd_def, run_time, hf, cp, g, zpe, e_dielec, mu_x, mu_y,mu_z, mu_total, num_scfs, cosmo_a, cosmo_v, vdw_a, g_prot, g_liginput, g_ligand, g_binding_qm, ki_qm, ki_type) select stic_id,method_id,%d,%lf,%lf,'%s',%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,'%s' from stics,methods where stics.compound_id in (select compound_id from compounds where sj_id = %d) and method like '%s';\n",rep->index,qm->G_binding, dock->Ki_DOCK,dock->Ki_unit,dock->E_internal,dock->E_nonbond,dock->E_electrostat,dock->E_internal,dock->G_tors,dock->E_unbound,dock->rmsd_ref,dock->time,qm->Hf,qm->Cp,qm->G,qm->ZPE,qm->E_dielec,qm->mu_x,qm->mu_y,qm->mu_z,qm->mu_total,qm->num_SCFs,qm->COSMO_A,qm->COSMO_V,qm->vdW_A,qm->G_prot,0.0,qm->G_ligand,qm->G_binding,qm->Ki_QM,qm->Ki_type,compound_id,params->qm_method);

  fprintf(sql_file,"insert into target_cm_results(target_id,cm_id) select target_id,cm_id from charge_models,targets where  cluster_rep_idx = %d and stic_id in (select stic_id from stics where compound_id in (select compound_id from compounds where sj_id = %d) and stereoisomer_idx = %d and tautomer_idx = %d and ionization_idx = %d and conformer_idx = %d) and target_name like '%s';\n",rep->index,compound_id,stic->S,stic->T,stic->I,stic->C,params->receptor_name);
}

static void FPrintSQLSTICTree(FILE *sql_file, int compound_id, RBTree *stics, JobParameters *params)
{
  struct STICelement *curr_stic = NULL;
  if(stics == NULL)
    return;
  if(sql_file == NULL)
    sql_file = stdout;

  FPrintSQLSTICTree(sql_file,compound_id,stics->left, params);
		    
  curr_stic = (struct STICelement*) stics->data;
  if(curr_stic != NULL)
    {
      struct ClusterRep *reps = curr_stic->reps;
      // Write SQL statements
      fprintf(sql_file,"insert into STICS(compound_id,stereoisomer_idx,tautomer_idx,ionization_idx,conformer_idx,charge) select compound_id,%d,%d,%d,%d,%lf from compounds where compounds.sj_id = %d  and not exists (select * from stics where compound_id in (select compound_id from compounds where sj_id = %d) and stereoisomer_idx = %d and tautomer_idx = %d and ionization_idx = %d and conformer_idx = %d);\n",curr_stic->S, curr_stic->T, curr_stic->I, curr_stic->C, curr_stic->charge, compound_id, compound_id, curr_stic->S, curr_stic->T, curr_stic->I, curr_stic->C);

      // Write Cluster rep SQL statements
      for(;reps != NULL;reps = reps->next)
	FPrintSQLClusterRep(sql_file, compound_id,reps, curr_stic, params);
    }

  FPrintSQLSTICTree(sql_file,compound_id,stics->right, params);
  
}

/**
 * Traverses a tree of compounds and writes SQL statements to store
 * STIC data in a SQL database.
 *
 * Designed to be safe enough, such that if an error occurs, Effercio can
 * continue rather than exit
 */
static void FPrintSQLCompoundTree(FILE *sql_file, RBTree *CompoundList, JobParameters *params)
{
  int compound_id = -1;
  CompoundTree *compound_tree = NULL;
  struct compoundData* compound = NULL;
  RBTree *stics = NULL;
  if(sql_file == NULL)
    sql_file = stdout;
  
  if(CompoundList == NULL || CompoundList->data == NULL)
    return;

  // Move into the left branch
  FPrintSQLCompoundTree(sql_file,CompoundList->left, params);

  compound_tree = (CompoundTree*) CompoundList->data;
  compound = (struct compoundData*) compound_tree->data;
  stics = (RBTree*) compound_tree->stics;
  if(compound->ID != NULL && stics != NULL)
    {
      int num_read = sscanf(compound->ID,"%d",&compound_id);
      if(num_read != 1)
	fprintf(stderr,"ERROR - Could not obtain SJ ID (integer) from compound: %s\n",compound->ID);
      else
	{
	  fprintf(sql_file,"insert into compounds(sj_id) values (%d);\n",compound_id);
	  FPrintSQLSTICTree(sql_file, compound_id, stics, params);// Traverse STIC tree
	}
    }
  // Move into the right branch
  FPrintSQLCompoundTree(sql_file,CompoundList->right, params);
  
}

void FPrintSQL(FILE *sql_file, RBTree *CompoundList, JobParameters *params)
{
  if(sql_file == NULL)
    sql_file = stdout;
  
  if(CompoundList == NULL)
    {
      fprintf(stderr,"WARNING - FPrintSQL was given a null pointer for the compound list.\n");
      return;
    }

  if(CompoundList->type != COMPOUND)
    {
      fprintf(stderr,"ERROR - PrintCompoundTree attempted on a non compound tree.\n");
      return;
    }

  fprintf(sql_file,"START TRANSACTION;\n");

  fprintf(sql_file,"insert into targets(target_name) select '%s' from (select 1 as value) as T where not exists (select * from targets where target_name like '%s');\n",params->receptor_name,params->receptor_name);
  fprintf(sql_file,"insert into methods(method) select '%s' from (select 1 as value) as T where not exists (select 1 from methods where method like '%s');\n",params->qm_method, params->qm_method);

  FPrintSQLCompoundTree(sql_file,CompoundList, params);
  fprintf(sql_file,"COMMIT;\n");
}
