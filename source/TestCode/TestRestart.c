/**
 * Test restoring a state using .restore files
 *
 * mpicc -DHAVE_CONFIG_H -I../../ -I.. -ggdb -O0 TestRestart.c ../state.c ../deque.c ../tpl.c ../structs.c ../io_utils.c ../RBTree.c ../analysis.c -lm
 */

#include "state.h"
#include <stdio.h>

int main () {
    cpunode *free_cpus = NULL;
    deque *jobs, *busy_list;
    RBTree *CompoundList = NULL;
    JobParameters params;
    int retval;

    jobs = InitDeque(NULL);
    busy_list = InitDeque(NULL);

    retval = RestoreState(&free_cpus, jobs, busy_list, &CompoundList, &params);
    if (retval) {
        printf("Restore successful\n");
    } else {
        printf("Error restoring state\n");
        return retval;
    }
    params.restart_job = 1;
    
    printf("Restored Jobs:\n");
    PrintJoblist(jobs->head);
    printf("\n\n");

    printf("Compound List:\n");
    FPrintCompoundTree(stdout, CompoundList);
    printf("\n\n");

    printf("For save test, renaming compound restart file\n");
    rename("Effercio.compound.restart", "Effercio.compound.restart.orig");

    SaveState(jobs, busy_list, CompoundList, &params);
    printf("Saved State\n");
    BoltzmannAvgCompoundTree(CompoundList, "/effercio/effercio_test/analysis", params.UseFreeEnergy);

    return retval;

}

