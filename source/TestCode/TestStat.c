// $Id: TestStat.c 6 2011-11-17 19:55:32Z dcoss $
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>

int main(void)
{

   char *file1name = "asldfkaslfdkjasl";
   char *file2name = "testval";
   struct dirent *file1, *file2;
   struct stat status_buffer;

   

   int i = stat(file1name,&status_buffer);
   int j = stat(file2name,&status_buffer);

   printf("i = %d    j = %d\n",i,j);
}
