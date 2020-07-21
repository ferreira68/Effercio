#ifndef IO_UTILS
#define IO_UTILS


void PrintHelp();
void PrintBanner();
void printf_master(const char *format, ...);
void printf_verbose(int verbose, const char *format, ...);

int VerifyDir(const char dirname[], int verbose, const char nodestr[]);
int CopyFile(const char source_name[], const char destination_name[]);
int MoveFile(const char from_name[], const char to_name[]);

void GetValueByColumn(const char *string, const char *key, int data_column, double *retval);
void GetStringByColumn(const char *string, const char *key, int data_column, char **retstr);

#endif

