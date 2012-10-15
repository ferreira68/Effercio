//
// $Id: structs.h 9 2012-03-01 20:01:25Z dcoss $
//

/**
 * cpunode and jobnodes were changed by David Coss, PhD on or after
 * 4 Mar 2011
 */

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
 *************************************************************************/
#ifndef EFFERCIO_STRUCTS
#define EFFERCIO_STRUCTS

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "tpl.h"

//Enumeration which indicates which process is to be run by the slave node.
enum JOB_TYPES{QUIT_SIGNAL = 0,PRESCREENING,DOCK,QM,QMLIG,NUM_JOB_TYPES};

//Error Codes used for uniform exit status
enum ERROR_CODES{UNKNOWN_ERROR = 1,RBTREE_SORT_ERROR,RBTREE_TYPE_ERROR,
		 NULL_POINTER,TPL_ERROR,FILESYS_ERROR,MOPAC_ERROR,DOCK_ERROR};

//Structure data type for stack used to track free CPUs.
typedef struct cpunode_t {
	int rank;
	struct cpunode_t *link;
} cpunode;


struct compoundData {
    char *ID;                       // Unique compound identifier
    int num_STIC;                   // Number of STIC enumerations for this compound
};


struct STICelement{
    int S;                    // index for stereoisomer
    int T;                    // index for tautomer
    int I;                    // index for ionization state
    int C;                    // index for conformer
    double charge;            // charge on this ligand
    double total_charge;      // total charge of the protein/ligand complex
    double Hf;                // Heat of formation for this STIC
    double G;                 // Free energy for this STIC
    struct ClusterRep *reps;
};

struct DOCKresult {
    double G_binding;      // Estimated Free Energy of Binding *Used for Boltzmann*
    double Ki_DOCK;        // Estimated Inhibition Constant, Ki
    double E_inter;        // Final Intermolecular Energy
    double E_nonbond;      // vdW + Hbond + desolv Energy
    double E_electrostat;  // Electrostatic Energy
    double E_internal;     // Final Total Internal Energy
    double G_tors;         // Torsional Free Energy
    double E_unbound;      // Unbound System's Energy
    double rmsd_ref;       // RMSD from reference structure
    double time;           // Time for docking run in seconds
    char *Ki_unit;
};


struct QMresult {
    double Hf;        // Heat of Formation *Used for Boltzmann*
    double S;         // Entropy
    double Cp;        // Heat capacity at constant pressure
    double G;         // Free Energy of Complex, from Mopac
    double ZPE;       // Zero-point energy
    double E_dielec;  // dielectric energy
    double mu_x,mu_y,mu_z;     // Components of the dipole moment {x,y,z}
    double mu_total;  // Total Dipole
    double time;      // Time for calclation in seconds
    int num_SCFs;     // Number of SCF iterations
    double COSMO_A;   // COSMO area for solvation model
    double COSMO_V;   // COSMO volume for solvation model
    double vdW_A;     // van der Waals surface area
    double G_prot;    // Free Energy of protein, from input deck
    //double G_liginput;// Free Energy of ligand input structure
    double G_ligand;  // Free Energy of ligand, from Mopac
    double G_binding; // Estimated Free Energy of Binding, calculated by G_binding := G_products - G_reactants *Used for Boltzmann*
    double Ki_QM;     // Estimated Inhibition Constant
    char *method;
    char *Ki_type;
};

struct ClusterRep {
    int index;                   // Index for cluster from Docking
    int size;                    // Number of structures in cluster
//    char *DOCK_Ki_unit;
//    char *QM_method;
//    char *QM_Ki_type;
    struct DOCKresult docked;   // Pointer to AutoDock results
    struct QMresult optimized;  // Pointer to MOPAC results
    struct ClusterRep *next;
};

//Data used for averaging and sorting STICs
enum AVG_SORT_CRIT{AVG_SORT_DOCK,AVG_SORT_QM};
struct AvgData{
	double Ki_DOCK,Ki_QM;
	int sort_criterion;
};
struct STICAvg {
	double Z;
	struct AvgData Ki;
	struct STICelement *stic;
};


enum input_units{INPUT_KCAL = 0,INPUT_KJ,INPUT_J,INPUT_CAL};

typedef struct
{
	char receptor_dir[FILENAME_MAX];
	char receptor_name[FILENAME_MAX];
	char ligand_dir[FILENAME_MAX];
	char results_dir[FILENAME_MAX];
	char optimized_dir[FILENAME_MAX];
	char clusters_dir[FILENAME_MAX];
	char analysis_dir[FILENAME_MAX];
	char **dock_params;
	int num_dock_params;
	char **mopac_header_params;
	int num_mopac_header_params;
	char **mopac_footer_params;
	int num_mopac_footer_params;
	char *qm_method;
	char username[LOGIN_NAME_MAX],node_tag[HOST_NAME_MAX];
	int transfer_node;
	int verbose;
	int use_mozyme;
	int doMOPAC;
	int restart_job;
	int prescreen;
	int UseFreeEnergy;
	int input_units;
	time_t start_time;
	time_t total_time;// Total time allowed to run, if specified.

	double target_Hf;// Target data provided from input deck
	double target_G;
}JobParameters;

//Job type. Contains information of jobs to be run by the slave nodes.
//Information includes name of the job (taken from the queue file).
//Type of job (see JOB_TYPE).
//And MPI rank of the CPU running the job.
typedef struct {
    char   *name; // Compound file name without extension
    char *mopac_res;// restart file for mopac
    int type;// Type of job to be run. This will determine the process taken by the slave node.
    int host_rank;// MPI rank of node running job, if it has been assigned. Otherwise, negative.
    struct STICelement input_data;// Data about the STIC provided by the ligand input file, i.e. enthalpy and free energy.
}job_t;



//Declaring functions avoids "assignment makes pointer from integer
//without a cast" compilation warns and potential problems later on.
struct compoundData *InitCompound(struct compoundData *data);
struct STICelement *InitSTIC(struct STICelement *data);
struct ClusterRep *InitClusterRep(struct ClusterRep *data);
cpunode* PopCPUNode(cpunode **stack);
cpunode* PushCPUNode(cpunode **stack,int rank);

#endif
