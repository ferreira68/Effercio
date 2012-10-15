//
// $Id: deque.h 6 2011-11-17 19:55:32Z dcoss $
//
#ifndef DEQUE_H
#define DEQUE_H 1

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

/**
 * Returns the i-th element in the deque from the head. Head is the 0-th
 * element
 */
deque_node* at(deque *the_deque, size_t i);

/**
 * Returns the size of the deque
 */
size_t deque_size(deque *the_deque);

#endif

