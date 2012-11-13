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

#define PROG_NAME "Effercio"
/*#define PROG_NAME "AutodockMPI"*/
#define MEM_BLOCK_SIZE 1024
#define OUTPUT_WIDTH 80
#define STRUCTS_MAX 25
#define TRUE   1
#define FALSE  0
#define MASTER 0
#define DEFAULT_ERR_CODE 1
#define THE_DIGITS "0123456789"
#define ZERO 0.0E0
#define ONE 1.0E0
#define FLOAT_TRANSFER_COUNT 2 //number of floats sent to slave nodes for stic data from input file.
#define CHARGE_TOL 2.5E-1
#define DAY_IN_SEC 8.64E4
#define FILE_BUFF_SIZE 2097152 // Added to help with compile. Need to remove once it's no longer needed
//#define DOUBLE_INIT 0xABCDEF
#define DOUBLE_INIT -9876543.21
#define INT_INIT -123456
#define STR_INIT "Undefined"
// Used for pretty output
#define DASH_HDR "-------------------------------------------------------------------------------"
#define EQ_HDR "==============================================================================="
#define WIGL_HDR "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
#define PLUS_HDR "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
#define AT_HDR "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
#define HASH_HDR "###############################################################################"
#define PCNT_HDR "\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%\%"
// AutoDock
#define DOCK_RESULTS_HDR "LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER"
#define CLUSTER_MARKER "USER    Cluster Rank ="
#define CLUSTER_SIZE_MARKER "Number of conformations in this cluster ="
#define MODEL_START "MODEL"
#define MODEL_END "ENDMDL"
#define G_binding_MARKER "Estimated Free Energy of Binding"
#define Ki_EST_MARKER "Estimated Inhibition Constant, Ki"
#define E_INTER_MARKER "Final Intermolecular Energy"
#define E_NONBOND_MARKER "vdW + Hbond + desolv Energy"
#define E_ELECTROSTAT_MARKER "Electrostatic Energy"
#define E_INTERNAL_MARKER "Final Total Internal Energy"
#define G_TORS_MARKER "Torsional Free Energy"
#define E_UNBOUND_MARKER "Unbound System's Energy"
#define DOCK_RMSD_MARKER "RMSD from reference structure       ="
#define DOCK_LIGAND_CHARGE "Total charge on ligand"
// MOPAC
#define MOPAC_DEBUG_TIME 10000
#define FREE_ENERGY_TIME_RATIO 0.86
#define MOPAC_KEYWDS_SOLV "EPS=78.40 RSOLV=1.30"
#define MOPAC_KEYWDS_CHRG "MULLIK"
#define MOPAC_CALC_CHARGE " COMPUTED CHARGE ON SYSTEM:"
#define MOPAC_EXCESS_CYCLES "EXCESS NUMBER OF OPTIMIZATION CYCLES"
#define MOPAC_LAMBDA_EXCESS "TOO MANY ITERATIONS IN LAMDA BISECT"
#define MOPAC_GEOMETRY_ERROR "TO CONTINUE CALCULATION SPECIFY \"GEO-OK\""
#define MOPAC_ENTHALPY "FINAL HEAT OF FORMATION"
#define MOPAC_VDW_AREA "VAN DER WAALS AREA"
#define MOPAC_E_DIELEC "DIELECTRIC ENERGY"
#define MOPAC_COSMO_AREA "COSMO AREA"
#define MOPAC_COSMO_VOL "COSMO VOLUME"
#define MOPAC_numSCF "SCF CALCULATIONS"
#define MOPAC_DIPOLE "DIPOLE           X         Y         Z       TOTAL"
#define MOPAC_TIME "TOTAL CPU TIME:"
#define MOPAC_SUCCESSFUL " == MOPAC DONE =="
#define MOPAC_ZPE_MARKER "ZERO POINT ENERGY"
#define MOPAC_THERMO_MARKER "CALCULATED THERMODYNAMIC PROPERTIES"
#define MOPAC_NO_TIME_LEFT "THERE IS NOT ENOUGH TIME FOR ANOTHER CYCLE"
#define MOPAC_GEOM_UNREC "STRUCTURE UNRECOGNIZABLE"
// TPL
#define STIC_TPL_MAP "A(S(iiiiffff))A(S(ii$(ffffffffffs)$(fffffffffffifffffffss)))"
#define COMPOUND_TPL_MAP "A(S(si))"
#define JOBLIST_TPL_MAP "A(S(ssii$(iiiiffff)))"
#define JOBS_RESTART_FILE "Effercio.jobs.restart"
#define BUSY_RESTART_FILE "Effercio.busy.restart"
#define COMPOUND_RESTART_FILE "Effercio.compound.restart"
#define STIC_RESTART_FILE "Effercio.stic.restart"

#define EFFERCIO_KT 0.59 //kcal/mol

#ifndef ARG_MAX
#define ARG_MAX 4096
#endif

