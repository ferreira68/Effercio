// $Id: DynamicArrays.c 6 2011-11-17 19:55:32Z dcoss $
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


int main(void)
{
    char **test;
    int count = 0;
    int i;
    char *str1 = "This is a test";
    char *str2 = "This is another test";
    char *str3 = "Si hoc legere scis nimium eruditiones habes";

    test = malloc(sizeof(char **) + sizeof(char *) + 1);
    (*test) = NULL;

    add_strtoarray(&test,&count,str1);
    add_strtoarray(&test,&count,str2);
    add_strtoarray(&test,&count,str3);

    for (i=0;i<count;i++) {
	printf("DEBUG - test[%d] = %s (%d chars)\n",i,test[i],strlen(test[i]));
    }

    for (i=0;i<count;i++) {
	free(test[i]);
        test[i] = NULL;
    }

    free(test);
    test = NULL;
}
