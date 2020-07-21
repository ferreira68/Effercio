#ifndef DYNAMIC_ARRAYS
#define DYNAMIC_ARRAYS 

void add_strtoarray(char ***str_array, int *num_elements, const char *my_string);
void add_inttoarray(int **int_array, int *num_elements, const int my_int);
void FreeStringArray(char ***str_array, int *num_elements);

#endif
