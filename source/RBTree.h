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
 * log(N) inserting, deleting and searching of job nodes in Effercio     *
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
void FreeRBTree(RBTree *node);

#endif //RBTREE_H
