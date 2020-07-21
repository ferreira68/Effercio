#ifndef DEQUE_H
#define DEQUE_H 1

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


#include "structs.h"

typedef enum {DEQUE_JOB = 0,DEQUE_DOUBLE,DEQUE_NUM_TYPES}deque_type;

typedef struct deque_node_t{
	deque_type type;
	void *data;
	struct deque_node_t *left, *right;
}deque_node;

typedef struct {
	deque_node *head, *tail;
}deque;

/**
 * Initializes a new Deque, with null pointers. An empty deque.
 */
deque* InitDeque(deque *newdeque);

deque_node* InitDequeNode(deque_node *node);

/**
 * Adds data to the front of the Deque. The head of the deque will
 * equal data.
 *
 * Returns a pointer to the new head of the deque, which equals data.
 */
deque_node* push_front(deque *the_deque, deque_node *data);


/**
 * Adds data to the back of the Deque. The tail of the deque will
 * equal data.
 *
 * Returns a pointer to the new tail of the deque, which equals data.
 */
deque_node* push_back(deque *the_deque, deque_node *data);

/**
 * Returns a pointer to the head of the Deque and removes it from the deque.
 *
 * Returns pointer to the head of the deque.
 */
deque_node* pop_front(deque *the_deque);

/**
 * Returns a pointer to the tail of the deque and removes it from the deque.
 *
 */
deque_node* pop_back(deque *the_deque);

deque_node* remove_node(deque *the_deque, deque_node *node);
void clear_deque(deque *the_deque);

/**
 * Returns the i-th element in the deque from the head. Head is the 0-th
 * element
 */
deque_node* at(deque *the_deque, size_t i);

/**
 * Returns the size of the deque
 */
size_t deque_size(const deque *the_deque);

#endif

