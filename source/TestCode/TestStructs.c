#include "mpi.h"
#include "../structs.h"
#include "../RBTree.h"
#include "../memwatch.h"
#include "../defines.h"
#include "../effercio_db.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


int main(void)
{
    // DOCK Routines tested and working on 29 Mar 2011 by DRC
    // No memory leaks
  /**/
  size_t i;
    struct DOCKresult dock;
    SetupDOCK(&dock);
    dock.Ki_DOCK = 42;
    strcpy(dock.Ki_unit,"nM");
    FPrintDOCK(stdout,&dock,"uM");


    // QM Routines tested and working on 15SEP10 by AMF
    // No memory leaks
    struct QMresult qm;

    SetupQM(&qm);
    qm.G_ligand = 123;
    FPrintQM(stdout,&qm,"SEQM","uM");


    // ClusterRep Routines tested and working on 15SEP10 by AMF
    // No memory leaks
    struct ClusterRep *clustptr;

    clustptr = InitClusterRep(clustptr);
    clustptr->docked = dock;
    clustptr->optimized = qm;
    FPrintClusterRep(stdout,clustptr);
    // FreeClusterRep(clustptr);

   
    // STIC Routines tested and working on 15SEP10 by AMF
    // No memory leaks
    struct STICelement STICdata, *STICptr, *unpacked;

    STICptr = InitSTIC(STICptr);
    STICptr->reps = clustptr;
    FPrintSTIC(stdout,STICptr);

    void *buffer;
    size_t buffer_size;
    PackBufferSTIC(&buffer,&buffer_size,"the job",0,STICptr);
    FILE* test_file = fopen("packed_stic.tpl","w");
    fwrite(buffer,1,buffer_size,test_file);
    fclose(test_file);
    
    printf("UNPACKING\n");
    unpacked = InitSTIC(unpacked);
    UnpackBufferSTIC(buffer,buffer_size,"the job",0,unpacked);
    FPrintSTIC(stdout,unpacked);

    // Test sql printing
    RBTree *CompoundList = NULL;
    JobParameters params;
    char method[] = "PM6";
    FILE *sql_file = fopen("sql_file.sql","w");
    if(sql_file == NULL)
      {
	fprintf(stderr,"Could not open \"sql_file.sql\".\n");
	exit(-1);
      }
    strcat(params.receptor_name,"RCPTR");
    params.qm_method = method;
    MergeSTICData("123456789",unpacked,&CompoundList);
    FPrintSQL(sql_file,CompoundList,&params);
    fclose(sql_file);
    /**/


    // MergeSTIC routine tested and working on 16SEP10 by AMF
    // No memory leaks
    /*/
    struct compoundData *test_cmpd_ptr;
    test_cmpd_ptr = InitCompound(test_cmpd_ptr);
    struct STICelement *test_STIC_ptr, *test_STIC_ptr2;
    test_STIC_ptr = InitSTIC(test_STIC_ptr);
    char *Name1 = "Compound 1";
    char *Name2 = "Compound 2";
    RBTree *CompoundList = NULL;

    
    test_STIC_ptr->S = 100;
    test_STIC_ptr->T = 101;
    test_STIC_ptr->I = 102;
    test_STIC_ptr->C = 103;
    test_STIC_ptr->charge = 104;
    test_STIC_ptr->total_charge = 105;
    test_STIC_ptr->G = -106;
    test_STIC_ptr->Hf = 107;
    test_STIC_ptr->reps = InitClusterRep(test_STIC_ptr->reps);
    test_STIC_ptr->reps->index = 106;
    test_STIC_ptr->reps->size = 107;
    test_STIC_ptr->reps->docked.G_binding = -108;
    test_STIC_ptr->reps->docked.Ki_DOCK = 109;
    strcpy(test_STIC_ptr->reps->docked.Ki_unit,"mM");
    test_STIC_ptr->reps->docked.E_inter = 110;
    test_STIC_ptr->reps->docked.E_nonbond = 111;
    test_STIC_ptr->reps->docked.E_electrostat = 112;
    test_STIC_ptr->reps->docked.E_internal = 113;
    test_STIC_ptr->reps->docked.G_tors = 114;
    test_STIC_ptr->reps->docked.E_unbound = 115;
    test_STIC_ptr->reps->docked.rmsd_ref = 116;
    test_STIC_ptr->reps->docked.time = 117;
    test_STIC_ptr->reps->optimized.Hf = -118;
    test_STIC_ptr->reps->optimized.S = 119;
    test_STIC_ptr->reps->optimized.Cp = 120;
    test_STIC_ptr->reps->optimized.G = 121;
    test_STIC_ptr->reps->optimized.ZPE = 122;
    test_STIC_ptr->reps->optimized.E_dielec = 123;
    test_STIC_ptr->reps->optimized.mu_x = 124;
    test_STIC_ptr->reps->optimized.mu_y = 125;
    test_STIC_ptr->reps->optimized.mu_z = 126;
    test_STIC_ptr->reps->optimized.mu_total = 127;
    test_STIC_ptr->reps->optimized.time = 128;
    strcpy(test_STIC_ptr->reps->optimized.method,"method   ");
    test_STIC_ptr->reps->optimized.num_SCFs = 129;
    test_STIC_ptr->reps->optimized.COSMO_A = 130;
    test_STIC_ptr->reps->optimized.COSMO_V = 131;
    test_STIC_ptr->reps->optimized.vdW_A = 132;
    test_STIC_ptr->reps->optimized.G_prot = 133;
    test_STIC_ptr->reps->optimized.G_ligand = 135;
    test_STIC_ptr->reps->optimized.G_binding = -136;
    test_STIC_ptr->reps->optimized.Ki_QM = 137;
    strcpy(test_STIC_ptr->reps->optimized.Ki_type,"Ki_type  ");
    test_STIC_ptr->reps->next = InitClusterRep(test_STIC_ptr->reps->next);
    test_STIC_ptr->reps->next->index = 206;
    test_STIC_ptr->reps->next->size = 207;
    test_STIC_ptr->reps->next->docked.G_binding = -208;
    test_STIC_ptr->reps->next->docked.Ki_DOCK = 209;
    strcpy(test_STIC_ptr->reps->next->docked.Ki_unit,"mM");
    test_STIC_ptr->reps->next->docked.E_inter = 210;
    test_STIC_ptr->reps->next->docked.E_nonbond = 211;
    test_STIC_ptr->reps->next->docked.E_electrostat = 212;
    test_STIC_ptr->reps->next->docked.E_internal = 213;
    test_STIC_ptr->reps->next->docked.G_tors = 214;
    test_STIC_ptr->reps->next->docked.E_unbound = 215;
    test_STIC_ptr->reps->next->docked.rmsd_ref = 216;
    test_STIC_ptr->reps->next->docked.time = 217;
    test_STIC_ptr->reps->next->optimized.Hf = -218;
    test_STIC_ptr->reps->next->optimized.S = 219;
    test_STIC_ptr->reps->next->optimized.Cp = 220;
    test_STIC_ptr->reps->next->optimized.G = 221;
    test_STIC_ptr->reps->next->optimized.ZPE = 222;
    test_STIC_ptr->reps->next->optimized.E_dielec = 223;
    test_STIC_ptr->reps->next->optimized.mu_x = 224;
    test_STIC_ptr->reps->next->optimized.mu_y = 225;
    test_STIC_ptr->reps->next->optimized.mu_z = 226;
    test_STIC_ptr->reps->next->optimized.mu_total = 227;
    test_STIC_ptr->reps->next->optimized.time = 228;
    strcpy(test_STIC_ptr->reps->next->optimized.method,"method   ");
    test_STIC_ptr->reps->next->optimized.num_SCFs = 229;
    test_STIC_ptr->reps->next->optimized.COSMO_A = 230;
    test_STIC_ptr->reps->next->optimized.COSMO_V = 231;
    test_STIC_ptr->reps->next->optimized.vdW_A = 232;
    test_STIC_ptr->reps->next->optimized.G_prot = 233;
    test_STIC_ptr->reps->next->optimized.G_ligand = 235;
    test_STIC_ptr->reps->next->optimized.G_binding = -236;
    test_STIC_ptr->reps->next->optimized.Ki_QM = 237;
    strcpy(test_STIC_ptr->reps->next->optimized.Ki_type,"Ki_type  ");
    MergeSTICData(Name1,test_STIC_ptr,&CompoundList);
    MergeSTICData(Name2,test_STIC_ptr,&CompoundList);
    FreeSTIC(test_STIC_ptr);

    test_STIC_ptr2 = InitSTIC(test_STIC_ptr2);

    test_STIC_ptr2->S = 300;
    test_STIC_ptr2->T = 301;
    test_STIC_ptr2->I = 302;
    test_STIC_ptr2->C = 303;
    test_STIC_ptr2->charge = 304;
    test_STIC_ptr2->total_charge = 305;
    test_STIC_ptr2->G = 306;
    test_STIC_ptr2->Hf = 307;
    test_STIC_ptr2->reps = InitClusterRep(test_STIC_ptr2->reps);
    test_STIC_ptr2->reps->index = 306;
    test_STIC_ptr2->reps->size = 307;
    test_STIC_ptr2->reps->docked.G_binding = 308;
    test_STIC_ptr2->reps->docked.Ki_DOCK = 309;
    strcpy(test_STIC_ptr2->reps->docked.Ki_unit,"mM");
    test_STIC_ptr2->reps->docked.E_inter = 310;
    test_STIC_ptr2->reps->docked.E_nonbond = 311;
    test_STIC_ptr2->reps->docked.E_electrostat = 312;
    test_STIC_ptr2->reps->docked.E_internal = 313;
    test_STIC_ptr2->reps->docked.G_tors = 314;
    test_STIC_ptr2->reps->docked.E_unbound = 315;
    test_STIC_ptr2->reps->docked.rmsd_ref = 316;
    test_STIC_ptr2->reps->docked.time = 317;
    test_STIC_ptr2->reps->optimized.Hf = 318;
    test_STIC_ptr2->reps->optimized.S = 319;
    test_STIC_ptr2->reps->optimized.Cp = 320;
    test_STIC_ptr2->reps->optimized.G = 321;
    test_STIC_ptr2->reps->optimized.ZPE = 322;
    test_STIC_ptr2->reps->optimized.E_dielec = 323;
    test_STIC_ptr2->reps->optimized.mu_x = 324;
    test_STIC_ptr2->reps->optimized.mu_y = 325;
    test_STIC_ptr2->reps->optimized.mu_z = 326;
    test_STIC_ptr2->reps->optimized.mu_total = 327;
    test_STIC_ptr2->reps->optimized.time = 328;
    strcpy(test_STIC_ptr2->reps->optimized.method,"method   ");
    test_STIC_ptr2->reps->optimized.num_SCFs = 329;
    test_STIC_ptr2->reps->optimized.COSMO_A = 330;
    test_STIC_ptr2->reps->optimized.COSMO_V = 331;
    test_STIC_ptr2->reps->optimized.vdW_A = 332;
    test_STIC_ptr2->reps->optimized.G_prot = 333;
    test_STIC_ptr2->reps->optimized.G_ligand = 335;
    test_STIC_ptr2->reps->optimized.G_binding = 336;
    test_STIC_ptr2->reps->optimized.Ki_QM = 337;
    strcpy(test_STIC_ptr2->reps->optimized.Ki_type,"Ki_type  ");
    test_STIC_ptr2->reps->next = InitClusterRep(test_STIC_ptr2->reps->next);
    test_STIC_ptr2->reps->next->index = 406;
    test_STIC_ptr2->reps->next->size = 407;
    test_STIC_ptr2->reps->next->docked.G_binding = 408;
    test_STIC_ptr2->reps->next->docked.Ki_DOCK = 409;
    strcpy(test_STIC_ptr2->reps->next->docked.Ki_unit,"mM");
    test_STIC_ptr2->reps->next->docked.E_inter = 410;
    test_STIC_ptr2->reps->next->docked.E_nonbond = 411;
    test_STIC_ptr2->reps->next->docked.E_electrostat = 412;
    test_STIC_ptr2->reps->next->docked.E_internal = 413;
    test_STIC_ptr2->reps->next->docked.G_tors = 414;
    test_STIC_ptr2->reps->next->docked.E_unbound = 415;
    test_STIC_ptr2->reps->next->docked.rmsd_ref = 416;
    test_STIC_ptr2->reps->next->docked.time = 417;
    test_STIC_ptr2->reps->next->optimized.Hf = 418;
    test_STIC_ptr2->reps->next->optimized.S = 419;
    test_STIC_ptr2->reps->next->optimized.Cp = 420;
    test_STIC_ptr2->reps->next->optimized.G = 421;
    test_STIC_ptr2->reps->next->optimized.ZPE = 422;
    test_STIC_ptr2->reps->next->optimized.E_dielec = 423;
    test_STIC_ptr2->reps->next->optimized.mu_x = 424;
    test_STIC_ptr2->reps->next->optimized.mu_y = 425;
    test_STIC_ptr2->reps->next->optimized.mu_z = 426;
    test_STIC_ptr2->reps->next->optimized.mu_total = 427;
    test_STIC_ptr2->reps->next->optimized.time = 428;
    strcpy(test_STIC_ptr2->reps->next->optimized.method,"method   ");
    test_STIC_ptr2->reps->next->optimized.num_SCFs = 429;
    test_STIC_ptr2->reps->next->optimized.COSMO_A = 430;
    test_STIC_ptr2->reps->next->optimized.COSMO_V = 431;
    test_STIC_ptr2->reps->next->optimized.vdW_A = 432;
    test_STIC_ptr2->reps->next->optimized.G_prot = 433;
    test_STIC_ptr2->reps->next->optimized.G_ligand = 435;
    test_STIC_ptr2->reps->next->optimized.G_binding = 436;
    test_STIC_ptr2->reps->next->optimized.Ki_QM = 437;

    //MergeSTICData(Name1,test_STIC_ptr2,&CompoundList);
    //MergeSTICData(Name2,test_STIC_ptr2,&CompoundList);
    //PrintCompoundTree(CompoundList);

    size_t i = 0;
    for(;i<9;i++)
      {
	struct STICelement * newstic;
	newstic = InitSTIC(newstic);
	newstic->S = 100+  i;
	newstic->T = 100;
	newstic->I = 100;
	newstic->C = 100;
	newstic->reps = InitClusterRep(newstic->reps);
	newstic->reps->optimized.G = i*2;
	MergeSTICData("Compound 1",newstic,&CompoundList);
      }
    MergeSTICData(Name1,test_STIC_ptr2,&CompoundList);
    MergeSTICData(Name2,test_STIC_ptr2,&CompoundList);
    PrintCompoundTree(CompoundList);

    FreeRBTree(CompoundList);
    FreeCompound(test_cmpd_ptr);
    FreeSTIC(test_STIC_ptr2);
    /**/

    //Test that would run before I started working on structs --drc 29 Mar 2011
/*
    struct STICelement *TestSTIC;
    struct compoundData *MolList;
    char *MolName1,*MolName2;

    TestSTIC = InitSTIC(TestSTIC);
    MolList  = InitCompound(MolList);
    MolName1 = strdup("SJ172550-1");
    MolName2 = strdup("SJ001113-1");

    TestSTIC->S = 1;
    TestSTIC->T = 2;
    TestSTIC->I = 3;
    TestSTIC->C = 4;
    TestSTIC->Hf = -89.39203;

    MergeSTICData(MolName1,TestSTIC,&MolList);

    TestSTIC->C = 5;
    MergeSTICData(MolName1,TestSTIC,&MolList);


    TestSTIC->S = 4;
    TestSTIC->T = 3;
    TestSTIC->I = 2;
    TestSTIC->C = 1;
    TestSTIC->Hf = -89.39203;

    MergeSTICData(MolName2,TestSTIC,&MolList);

    TestSTIC->S = 1;
    TestSTIC->T = 2;
    TestSTIC->I = 3;
    TestSTIC->C = 4;
    TestSTIC->Hf = -89.39203;
    TestSTIC->reps->optimized.Hf = -99999;
    MergeSTICData(MolName1,TestSTIC,&MolList);
    
    PrintCompounds(MolList);

    FreeSTIC(TestSTIC);
    FreeCompound(MolList);
    free(MolName1);
    free(MolName2);
    /**/

/*
    STICptr = InitSTIC(STICptr);
    size_t BUFFER_SIZE = (size_t) (BUFSIZ * MEM_BLOCK_SIZE);
    char *buffer = malloc(BUFFER_SIZE);
    size_t numbytes = 0;
    char *jobname = "Testing";
    int retval = 0;
    PackBufferSTIC(buffer,jobname,retval,STICptr,&numbytes);
    struct STICelement *tempSTIC;
    tempSTIC = InitSTIC(tempSTIC);
    char *mystr = malloc(sizeof(char *));
    int myint;
    UnpackBufferSTIC(buffer,&mystr,&myint,tempSTIC);
    PrintSTIC(tempSTIC);
    FreeSTIC(tempSTIC);
    free(buffer);
    FreeSTIC(STICptr);
    free(mystr);
*/

    return 0;
}
