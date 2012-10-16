//
// $Id: RBTree.h 9 2012-03-01 20:01:25Z dcoss $
//
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
 * This code sets up and operates on a Red-Black Tree, allow efficient   *
 * log(N) inserting, deleting and searching of job nodes in Effercio    *
 *************************************************************************/
#ifndef RBTREE_H
#define RBTREE_H 

#include <stdio.h>
#include "structs.h"

/**
 * Tree node for Compounds for a particular Effercio run.
 * Insert sort by Compound ID
*/
typedef struct {
	struct compoundData *data;
	struct RBTreeNode *stics;
}CompoundTree;

struct CompoundAvg {
	struct  AvgData Ki;
	CompoundTree* compound;// Pointer to the actual compound data (Compound data and STICs).
};

struct ClusterRepAvg{
	struct AvgData Ki;
	double Z;
	struct ClusterRep* rep;
};

typedef enum{JOB=0,COMPOUND,STIC,STIC_AVG,COMPOUND_AVG,CR_AVG,NUM_TREE_TYPES}RBTree_type;

typedef enum {RBT_BLACK = 0,RBT_RED}RBTree_Color;
typedef struct RBTreeNode {
	RBTree_Color color;//red or black.
	void *data;
	RBTree_type type;
	struct RBTreeNode *left,*right,*parent;
}RBTree;


RBTree*  InitRBTree(RBTree *node,int node_type);
RBTree* FirstRBTNode(RBTree* tree);
RBTree* insertRBT(RBTree **ptrRoot,RBTree *x);
RBTree* FindRBTNode(RBTree *tree, void *new_data);

#endif //RBTREE_H
