// $Id: help.c 9 2012-03-01 20:01:25Z dcoss $
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

#include "defines.h"
#include "memwatch.h"

void PrintHelp()
{
	const char *mytabs = "\t";

	printf_master("\nUsage: %s [options] <receptor_name> <file>\n\n",PROG_NAME);
	printf_master("%swhere <receptor_name> is the name of the receptor .pdbqt file\n",mytabs);
	printf_master("%sand <file> contains a list of .pdbqt files without extensions\n",mytabs);
	printf_master("%sOne file name per line.\n\n",mytabs);
	printf_master("\n%sExample: %s 3DAB drug_names\n",mytabs,PROG_NAME);
	printf_master("%swould search for receptor files named 3DAB.xxx in the default\n",mytabs);
	printf_master("%sdirectory (see below) and run against all files listed in\n",mytabs);
	printf_master("%sdrug_names in the default ligands directory.\n\n",mytabs);

	printf_master("Options:\n");

	printf_master("  -h\n%sPrint this help message\n",mytabs);

	printf_master("  -v, --verbose\n%sVerbosb output\n",mytabs);

	printf_master("  -l, --ligands\n%sDirectory containing ligand files in .pdbqt format\n",mytabs);

	printf_master("  -m, --maps\n%sDirectory containing the receptor .pdbqt and .map files\n",mytabs);

	printf_master("  -d, --docks\n%sDirectory to store the AutoDock results\n",mytabs);

	printf_master("  -c, --cluster-reps\n%sDirectory to store the cluster representatives\n",mytabs);

	printf_master("  -a, --analysis\n%sDirectory to store the analysis reports\n",mytabs);

	printf_master("  -q, --quantum\n%sTurn on MOPAC refinement of the cluster respresentatives\n",mytabs);

	printf_master("  -o, --optimized\n%sDirectory to store the MOPAC geometry optimizations\n",mytabs);

	printf_master("  -R, --restart\n%sRestart a previously unfinished %s run.\n",mytabs,PROG_NAME);
	printf_master("%sThis requires a file named %s.restart\n",mytabs,PROG_NAME);

	printf_master("  -t, --time\n%sTotal time for the job in seconds (default is one year)\n",mytabs);

	printf_master("  -p, --param\n%sUser-specified parameter to be used in the AutoDock runs\n",mytabs);
	printf_master("%sEXAMPLE: -p ga_num_evals=100000 or -param ga_nam_evals=100000\n",mytabs);
	printf_master("%s        This would insert the following line into the .dpf file\n",mytabs);
	printf_master("%s        > ga_num_evals 2500000\n",mytabs);
	printf_master("  -k, --mopac\n%sUser-specified keywords to be added to MOPAC header *and* footer.\n",mytabs);
	printf_master("%sEXAMPLE: -k \"FMAT FOCK\"\n",mytabs);
	printf_master("%s        This would add FMAT and FOCK to the MOPAC header and footer lines.\n",mytabs);
	printf_master("      --mopac-header\n%sUser-specified keywords to be added only to the MOPAC header.\n",mytabs);
	printf_master("      --mopac-footer\n%sUser-specified keywords to be added only to the MOPAC footer.\n",mytabs);
	printf_master("\n");
}
