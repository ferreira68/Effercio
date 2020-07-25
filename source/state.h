#ifndef STATE_H
#define STATE_H

#include "deque.h"
#include "RBTree.h"
#include "structs.h"

void SaveState(deque *jobs, deque *busy_list, RBTree *CompoundList, JobParameters *params);

int RestoreState(cpunode **free_cpus, deque *jobs, deque *busy_list, RBTree **CompoundList, JobParameters *params);

#endif

