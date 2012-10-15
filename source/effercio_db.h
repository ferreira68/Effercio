// $Id: effercio_db.h 6 2011-11-17 19:55:32Z dcoss $


#ifndef EFFERCIO_DB_H
#define EFFERCIO_DB_H 1

#include "structs.h"
#include "RBTree.h"

#include <stdio.h>


/**
 * Writes SQL statements to load the entire compound tree 
 *
 * If sql_file is NULL, standard output is used.
 * If CompoundList is NULL, the function returns without printing anything,
 *    except a warning in standard error.
 */

void FPrintSQL(FILE *sql_file, RBTree *CompoundList, JobParameters *params);


#endif //EFFERCIO_DB_H
