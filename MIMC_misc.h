#define _MIMC2_MISC_

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifndef _STDINT_H_
#include <stdint.h>
#endif


//#ifndef _STRING_H_
//#include <string.h>
//#endif

char isLeapYear(int yint);
double get_datenum(char *string_d0);
float get_dt(char *string_d0,char *string_d1);
void getTimeStampStr(char *filename, char *str_timestamp);
