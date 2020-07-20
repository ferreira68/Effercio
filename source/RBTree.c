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
#include <stdlib.h>
#include <string.h>
#include "RBTree.h"

RBTree*  InitRBTree(RBTree *node,int node_type)
{
	// node = (RBTree*)malloc(sizeof(RBTree));
	node = (RBTree*)calloc(1,sizeof(RBTree));
	node->type = node_type;
	node->data = NULL;
	node->left = node->parent = node->right = NULL;
	node->color = RBT_BLACK;


	return node;
}

void FreeRBTree(RBTree *node)
{
	if(node == NULL)
		return;

	switch(node->type)
	{
	case COMPOUND:
		FreeCompound(((CompoundTree*)node->data)->data);
		FreeRBTree(((CompoundTree*)node->data)->stics);
	case JOB:case STIC:
		case STIC_AVG:case COMPOUND_AVG: case CR_AVG://<--- For these, since they point to data contained in other structures, the programmer is required to free their data before calling this function. CAVEAT PROGRAMMER!
		free(node->data);
		break;
	default:
		fprintf(stderr,"WARNING - FreeRBTree: Unknown tree node type: %d.",node->type);
		break;
	}
	node->data = NULL;
	FreeRBTree(node->left);
	FreeRBTree(node->right);

}

static RBTree* RightRotate(RBTree *tree)
{
	RBTree *left, *right;
	if(tree == NULL)
		return NULL;;
	if(tree->left == NULL)
		return NULL;;

	//Rotate
	left = tree->left;
	right = left->right;
	left->right = tree;
	tree->left = right;

	//Update parenthood
	left->parent = tree->parent;
	if(tree->parent != NULL)//double-check that tree is not a root
	{
		if(tree->parent->left == tree)
			tree->parent->left = left;
		else
			tree->parent->right = left;
	}

	tree->parent = left;
	if(right != NULL)
		right->parent = tree;
	return left;
}


static RBTree* LeftRotate(RBTree *tree)
{
	RBTree *left, *right;
	if(tree == NULL)
		return NULL;
	if(tree->right == NULL)
		return NULL;

	//Rotate
	right = tree->right;
	left = right->left;
	tree->right = left;
	right->left = tree;

	//Update parenthood
	right->parent = tree->parent;
	if(tree->parent != NULL)//double-check that tree is not a root
		if(tree->parent != NULL)
		{
			if(tree->parent->left == tree)
				tree->parent->left = right;
			else
				tree->parent->right = right;
		}
	tree->parent = right;
	if(left != NULL)
		left->parent = tree;
	return right;
}



static short RBTDirection_STIC(const struct STICelement *lstic, const struct STICelement *rstic)
{
	if(lstic->S > rstic->S)
		return -1;
	if(lstic->S < rstic->S)
		return 1;
	if(lstic->C > rstic->C)
	  return -1;
	if(lstic->C < rstic->C)
	  return 1;
	if(lstic->T > rstic->T)
	  return -1;
	if(lstic->T < rstic->T)
	  return 1;
	if(lstic->I > rstic->I)
	  return -1;
	if(lstic->I < rstic->I)
	  return 1;
	return 0;

}

static short AvgSort(struct AvgData *lhs, struct AvgData *rhs)
{
	if(lhs->sort_criterion == AVG_SORT_DOCK)
	{
		if(lhs->Ki_DOCK > rhs->Ki_DOCK)
			return -1;
		//if(lhs->Ki_DOCK < rhs->Ki_DOCK)
		return 1;
		//return 0;
	}
	else
	{
		if(lhs->Ki_QM > rhs->Ki_QM)
			return -1;
		//if(lhs->Ki_QM < rhs->Ki_QM)
		return 1;
		//return 0;
	}
}

static short DoubleSort(double lhs, double rhs)
{
  if(lhs > rhs)
    return -1;
  return 1;
}

static short RBTDirection_compound(const struct compoundData *lhs, const struct compoundData *rhs)
{
	const char *lname = lhs->ID;
	const char *rname = rhs->ID;

	//This order is necessary for alphabetical insertion.
	return (short)strcmp(rname,lname);
}

static short RBTDirection_job(const job_t *lhs, const job_t *rhs)
{
	const char *lname = lhs->name;
	const char *rname = rhs->name;

	return (short) -1*strcmp(rname,lname);
}

static short RBTDirection(const RBTree *lhs, const RBTree *rhs)
{
	if(lhs == NULL || rhs == NULL)
	{
		fprintf(stderr,"RBTDirection: NULL POINTER!");
		exit(NULL_POINTER);
	}

	if(lhs->type != rhs->type)
	{
		fprintf(stderr,"RBTDirection: Cannot sort different tree node types.");
		exit(RBTREE_SORT_ERROR);
	}

	switch(lhs->type)
	{
	case JOB:
		return RBTDirection_job((job_t*)lhs->data,(job_t*)rhs->data);
	case COMPOUND:
		return RBTDirection_compound(((CompoundTree*)lhs->data)->data,((CompoundTree*)rhs->data)->data);
	case STIC:
		return RBTDirection_STIC((struct STICelement*)lhs->data,(struct STICelement*)rhs->data);
	case STIC_AVG:	case CR_AVG:
		{
		  double lhs_Z, rhs_Z;
			if(lhs->type == STIC_AVG)
			{
				lhs_Z = ((struct STICAvg*)lhs->data)->Z;
				rhs_Z = ((struct STICAvg*)rhs->data)->Z;
			}
			else if(lhs->type == CR_AVG)
			{
				lhs_Z = ((struct ClusterRepAvg*)lhs->data)->Z;
				rhs_Z = ((struct ClusterRepAvg*)rhs->data)->Z;
			}
		  return -1.0*DoubleSort(lhs_Z,rhs_Z);
		}
	case COMPOUND_AVG:
	  {
	    struct AvgData *lhs_avg, *rhs_avg;
	    lhs_avg = &((struct CompoundAvg*)lhs->data)->Ki;
	    rhs_avg = &((struct CompoundAvg*)rhs->data)->Ki;
	    return AvgSort(lhs_avg,rhs_avg);
	  }
	default:
		break;
	}
	fprintf(stderr,"RBTDirection: Cannot sort different tree node types.");
	exit(RBTREE_SORT_ERROR);
}

static void insertBST(RBTree *rootnode,RBTree *x)
{
	int direction;
	if(rootnode == NULL)
	{
		rootnode = x;
		x->parent = NULL;
		return;
	}

	direction = RBTDirection(rootnode,x);
	if(direction < 0)//goes to the left
	{
		if(rootnode->left == NULL)
		{
			rootnode->left = x;
			x->parent = rootnode;
			return;
		}
		return insertBST(rootnode->left,x);
	}
	else if(direction > 0)//goes to the right
	{
		if(rootnode->right == NULL)
		{
			rootnode->right = x;
			x->parent = rootnode;
			return;
		}
		return insertBST(rootnode->right,x);
	}
	else
	{
		fprintf(stderr,"Warning: insertBST cannot insert when two nodes are the same.");fflush(stderr);
	}
	return;
}

/**
 * Inserts x into the Red-Black tree, whose root pointer
 * is located in rootaddr.
 *
 * Returns: new root node
 */
RBTree* insertRBT(RBTree **ptrRoot,RBTree *x)
{
	RBTree *root = *ptrRoot;
	insertBST(root,x);
	if(root == NULL || root == x)
	{
		x->color = RBT_BLACK;
		*ptrRoot = x;
		return x;
	}
	x->color = RBT_RED;
	while(x != root && x->parent->color == RBT_RED)
	{
		if(x->parent == x->parent->parent->left)
		{
			RBTree *uncle = x->parent->parent->right;
			if(uncle != NULL && uncle->color == RBT_RED)
			{
				x->parent->color = uncle->color = RBT_BLACK;
				x->parent->parent->color = RBT_RED;
				x = x->parent->parent;
			}
			else
			{
				if(x == x->parent->right)
				{
					x = x->parent;
					LeftRotate(x);
				}
				x->parent->color = RBT_BLACK;
				x->parent->parent->color = RBT_RED;
				RightRotate(x->parent->parent);
			}
		}
		else
		{
			RBTree *uncle = x->parent->parent->left;
			if(uncle != NULL && uncle->color == RBT_RED)
			{
				x->parent->color = uncle->color = RBT_BLACK;
				x->parent->parent->color = RBT_RED;
				x = x->parent->parent;
			}
			else
			{
				if(x == x->parent->left)
				{
					x = x->parent;
					RightRotate(x);
				}
				x->parent->color = RBT_BLACK;
				x->parent->parent->color = RBT_RED;
				LeftRotate(x->parent->parent);
			}
		}
	}

	//Make sure the real root of the tree is black
	while(root->parent != NULL)
		root = root->parent;
	root->color = RBT_BLACK;
	*ptrRoot = root;
	return root;
}


RBTree* FindRBTNode(RBTree *tree, void *key)
{
	short direction;
	if(tree == NULL)
		return NULL;

	switch(tree->type)
	{
	case JOB:
		direction = RBTDirection_job((job_t*)tree->data,(job_t*)key);
		break;
	case STIC:
		direction = RBTDirection_STIC((struct STICelement*)tree->data,(struct STICelement*)key);
		break;
	case COMPOUND:
		direction = RBTDirection_compound(((CompoundTree*)tree->data)->data,(struct compoundData*)key);
		break;
	case STIC_AVG:
	  return -1*DoubleSort(((struct STICAvg*)tree->data)->Z,((struct STICAvg*)key)->Z);
	case COMPOUND_AVG:
		return AvgSort(&((struct CompoundAvg*)tree->data)->Ki,&((struct CompoundAvg*)key)->Ki);
	case CR_AVG:
	  return -1*DoubleSort(((struct ClusterRepAvg*)tree->data)->Z,((struct ClusterRepAvg*)key)->Z);
	default:
		fprintf(stderr,"FindRBTNode: Cannot search for the following tree type: %d\n",tree->type);
		exit(RBTREE_TYPE_ERROR);
	}

	if(direction == 0)
		return tree;
	if(direction > 0)
		return FindRBTNode(tree->right,key);

	return FindRBTNode(tree->left,key);
}

static RBTree* RBTSuccessor(RBTree *node)
{
	RBTree *temp;
	if(node == NULL)
		return NULL;
	temp = node->right;

	//If the node is already a right most child node,
	//move up the tree and to the right.
	if(temp == NULL)
	{
		if(node->parent == NULL)//If node is root, then it has no successor
			return NULL;
		temp = node->parent;
		while(temp->right == node)
		{
			node = temp;
			temp = temp->parent;
		}
		if(temp->parent == NULL)
			return NULL;
		return temp;
	}

	//If the node has right-children,
	//find the left most right child.
	while(temp->left != NULL)
	{
		temp = temp->left;
	}
	return temp;
}

static int RBTIsBlack(RBTree *node)
{
	return node == NULL || node->color == RBT_BLACK;
}

static void CleanRBT(RBTree *root, RBTree *node)
{
	if(node == NULL || root == NULL)
		return;

	while(node->color == RBT_BLACK && node != root)
	{
		if(node == node->parent->left)
		{
			RBTree *sibling = node->parent->right;
			if(sibling->color == RBT_RED)
			{
				sibling->color = RBT_BLACK;
				node->parent->color = RBT_RED;
				LeftRotate(node->parent);
				sibling = node->parent->right;
			}
			if(RBTIsBlack(sibling->right) && RBTIsBlack(sibling->left))
			{
				sibling->color = RBT_RED;
				node = node->parent;
			}
			else
			{
				if(RBTIsBlack(sibling->right))
				{
					sibling->left->color = RBT_BLACK;
					sibling->color = RBT_RED;
					RightRotate(sibling);
					sibling = node->parent->right;
				}
				sibling->color = node->parent->color;
				node->parent->color = sibling->right->color = RBT_BLACK;
				LeftRotate(node->parent);
				node = root;
			}
		}
		else
		{
			RBTree *sibling = node->parent->left;
			if(sibling->color == RBT_RED)
			{
				sibling->color = RBT_BLACK;
				node->parent->color = RBT_RED;
				RightRotate(node->parent);
				sibling = node->parent->left;
			}
			if(RBTIsBlack(sibling->right) && RBTIsBlack(sibling->left))
			{
				sibling->color = RBT_RED;
				node = node->parent;
			}
			else
			{
				if(RBTIsBlack(sibling->left))
				{
					sibling->right->color = RBT_BLACK;
					sibling->color = RBT_RED;
					LeftRotate(sibling);
					sibling = node->parent->left;
				}
				sibling->color = node->parent->color;
				node->parent->color = sibling->left->color = RBT_BLACK;
				RightRotate(node->parent);
				node = root;
			}
		}
	}
	node->color = RBT_BLACK;
	if(root->parent != 0)
	{
		root = node;
		while(root->parent != 0)
			root = root->parent;
		root->color = RBT_BLACK;
	}

}

/**
 * Will remove node from the tree. The node is
 * then returned, so that it can be freed from
 * memory or used in some other manner.
 */
RBTree* removeRBT(RBTree **tree, RBTree *node)
{
	RBTree *x,*y, *sentinal = 0, *root = *tree;

	if(*tree == node && node->left == NULL && node->right == NULL)//Is this node the only node in the tree.
	{
		*tree = NULL;
		return node;
	}

	if(node->left == NULL || node->right == NULL)
		y = node;
	else
		y = RBTSuccessor(node);

	if(y->left == NULL)
		x = y->right;
	else
		x = y->left;


	if(x == NULL)
	{
		sentinal = (RBTree*)malloc(sizeof(RBTree));
		sentinal->left = sentinal->right = NULL;
		sentinal->color = RBT_BLACK;
		x = sentinal;
	}
	x->parent = y->parent;

	/*if(root == y->parent)
	{
		root->left = x;
	}
	else */if(y->parent != NULL)
	{
		if(y == y->parent->left)
			y->parent->left = x;
		else
			y->parent->right = x;
	}
	if(y != node)
	{
		if(y->color == RBT_BLACK)
			CleanRBT(root,x);
		y->left = node->left;
		y->right = node->right;
		y->parent = node->parent;
		y->color = node->color;
		if(node->left != NULL)
			node->left->parent = y;
		if(node->right != NULL)
			node->right->parent = y;
		if(node->parent != NULL)
		{
			if(node == node->parent->left)
				node->parent->left = y;
			else
				node->parent->right = y;
		}
	}
	else
		if(y->color == RBT_BLACK)
			CleanRBT(root,x);

	//Ensure tree still holds the root node address.
	root = x;
	while(root->parent != NULL)
		root = root->parent;
	*tree = root;

	//For the purposes of cleaning the RB Tree, a dummy sentinal node
	//may be used for x. If this is the case, x should be a null pointer
	if(x == sentinal)
	{
		if(x == x->parent->left)
			x->parent->left = 0;
		if(x == x->parent->right)
			x->parent->right = 0;
	}
	if(sentinal != NULL)
		free(sentinal);

	return node;
}


int CountNodes(const RBTree *tree)
{
	if(tree == NULL)
		return 0;

	return 1 + CountNodes(tree->left) + CountNodes(tree->right);
}

RBTree* FirstRBTNode(RBTree* tree)
{
	if(tree == NULL || tree->left == NULL)
		return tree;
	return FirstRBTNode(tree->left);
}

RBTree* LastRBTNode(RBTree *tree)
{
	if(tree == NULL || tree->right == NULL)
		return tree;
	return FirstRBTNode(tree->right);
}
