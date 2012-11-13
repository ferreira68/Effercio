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

#include <stdio.h>
#include "deque.h"

deque* InitDeque(deque *newdeque)
{
	newdeque = (deque*)malloc(sizeof(deque));
	newdeque->head = newdeque->tail = NULL;

	return newdeque;
}

deque_node* InitDequeNode(deque_node *node)
{
	node = (deque_node*)malloc(sizeof(deque_node));
	node->type = DEQUE_JOB;
	node->data = NULL;
	node->left = node->right = NULL;
	return node;
}

deque_node* push_front(deque *the_deque, deque_node *node)
{
	deque_node *new_head = node;
	if(the_deque == NULL)
	{
		fprintf(stderr,"deque::push_front: null deque pointer.");
		exit(NULL_POINTER);
	}
	if(node == NULL)
		return NULL;

	//Is the deque empty?
	if(the_deque->head == NULL)
	{
		the_deque->head = the_deque->tail = node;
		return node;
	}

	while(node->right != NULL)
		node = node->right;
	node->right = the_deque->head;
	the_deque->head->left = node;
	the_deque->head = new_head;
	return new_head;
}

deque_node* push_back(deque *the_deque, deque_node *node)
{
	deque_node *new_tail = node;
	if(the_deque == NULL)
	{
		fprintf(stderr,"deque::push_back: null deque pointer.");
		exit(NULL_POINTER);
	}
	if(node == NULL)
		return NULL;

	//Is the deque empty?
	if(the_deque->head == NULL)
	{
		the_deque->head = the_deque->tail = node;
		return node;
	}

	while(node->left != NULL)
		node = node->left;
	node->left = the_deque->tail;
	the_deque->tail->right = node;
	the_deque->tail = new_tail;

	return new_tail;
}

deque_node* pop_front(deque *the_deque)
{
	deque_node *return_node;
	if(the_deque == NULL)
	{
		fprintf(stderr,"deque::pop_front: null deque pointer.");
		exit(NULL_POINTER);
	}

	return_node = the_deque->head;
	if(return_node != NULL)
	{
		the_deque->head = return_node->right;
		return_node->left = return_node->right = NULL;
	}
	if(the_deque->head == NULL)
		the_deque->tail = NULL;

	return return_node;
}

deque_node* pop_back(deque *the_deque)
{
	deque_node *return_node;
	if(the_deque == NULL)
	{
		fprintf(stderr,"deque::pop_back: null deque pointer.");
		exit(NULL_POINTER);
	}

	return_node = the_deque->tail;
	if(return_node != NULL)
	{
		the_deque->tail = return_node->left;
		return_node->left = return_node->right = NULL;
	}
	if(the_deque->tail == NULL)
			the_deque->head = NULL;

	return return_node;
}

deque_node* at(deque *the_deque, size_t i)
{
	deque_node *curr = the_deque->head;
	size_t index = 0;

	if(the_deque == NULL)
	{
		fprintf(stderr,"deque::at: null deque pointer.");
		exit(NULL_POINTER);
	}

	while(curr != NULL)
	{
		if(index == i)
			return curr;
		curr = curr->right;
		index++;
	}

	//If we're here, i is greater than the size of the deque
	return NULL;
}

deque_node* remove_node(deque *the_deque, deque_node *node)
{
	deque_node *left, *right;
	if(node == NULL || the_deque == NULL)
		return NULL;

	left = node->left;
	right = node->right;
	if(left != NULL)
		left->right = right;
	if(right != NULL)
		right->left = left;
	if(the_deque->head == node)
		the_deque->head = right;
	if(the_deque->tail == node)
		the_deque->tail = left;

	node->left = node->right = NULL;
	return node;
}

size_t deque_size(deque *the_deque)
{
	deque_node *curr;
	size_t size = 0;
	if(the_deque == NULL)
		return size;
	curr = the_deque->head;
	while(curr != NULL)
	{
		size++;
		curr = curr->right;
	}
	return size;
}

void clear_deque(deque *the_deque)
{
	if(the_deque == NULL)
		return;

	while(the_deque->head != NULL)
		free(pop_front(the_deque));
}


