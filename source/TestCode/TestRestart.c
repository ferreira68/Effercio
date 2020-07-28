/**
 * Test restoring a state using .restore files
 *
 * mpicc -DHAVE_CONFIG_H -I../../ -I.. -ggdb -O0 -o TestRestart TestRestart.c ../state.c ../deque.c ../tpl.c ../structs.c ../io_utils.c ../RBTree.c ../analysis.c -lm
 */

#include "state.h"
#include <stdio.h>

int main () {
    cpunode *free_cpus = NULL;
    deque *jobs, *busy_list;
    RBTree *CompoundList = NULL;
    JobParameters params;
    int retval;
    struct STICelement results = {1,2,3,4,-42.0,-15,-1,-42,NULL};

    struct ClusterRep* testrep = InitClusterRep(NULL);
    results.reps = testrep;
    testrep->index = 12;
    testrep->size = 34;
    testrep->docked.G_binding = 99999;
    testrep->docked.Ki_unit = strdup("Foo");
    testrep->optimized.Hf = -42.0;
    testrep->optimized.method  = strdup("Foo");
    testrep->optimized.Ki_type = strdup("Foo");


    jobs = InitDeque(NULL);
    busy_list = InitDeque(NULL);

    retval = RestoreState(&free_cpus, jobs, busy_list, &CompoundList, &params);
    if (retval) {
        printf("Restore successful\n");
        retval = 0;
    } else {
        printf("Error restoring state\n");
        return 1;
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

    printf("Testing MergeSTICData\n");
    MergeSTICData(((CompoundTree*)CompoundList->data)->data->ID,&results,&CompoundList);

    printf("Testing BoltzmannAvgCompoundTree\n");
    BoltzmannAvgCompoundTree(CompoundList, "/effercio/effercio_test/analysis", params.UseFreeEnergy);

    FPrintCompoundTree(stdout, CompoundList);

    return retval;

}

