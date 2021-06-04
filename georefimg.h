#ifndef _STDINT_H_
#include <stdint.h>
#endif

#ifndef _GMA_
#include "GMA.h"
#endif


#define _GEOREFIMG_

typedef struct georefimg_uint8{
    double *x;
    double *y;
    GMA_uint8 *z;
} georefimg_uint8;

typedef struct georefimg_uint16{
    double *x;
    double *y;
    GMA_uint16 *z;
} georefimg_uint16;

typedef struct georefimg_uint32{
    double *x;
    double *y;
    GMA_uint32 *z;
} georefimg_uint32;

typedef struct georefimg_int32{
    double *x;
    double *y;
    GMA_int32 *z;
} georefimg_int32;

typedef struct georefimg_float{
    double *x;
    double *y;
    GMA_float *z;
} georefimg_float;

typedef struct georefimg_double{
    double *x;
    double *y;
    GMA_double *z;
} georefimg_double;

// Initialization routines
georefimg_uint8* georefimg_uint8_create(unsigned int size_row, unsigned int size_col);
georefimg_uint16* georefimg_uint16_create(unsigned int size_row, unsigned int size_col);
georefimg_uint32* georefimg_uint32_create(unsigned int size_row, unsigned int size_col);
georefimg_int32* georefimg_int32_create(unsigned int size_row, unsigned int size_col);
georefimg_float* georefimg_float_create(unsigned int size_row, unsigned int size_col);
georefimg_double* georefimg_double_create(unsigned int size_row, unsigned int size_col);

// De-allocation routines
void georefimg_uint8_destroy(georefimg_uint8* in);
void georefimg_uint16_destroy(georefimg_uint16* in);
void georefimg_uint32_destroy(georefimg_uint32* in);
void georefimg_int32_destroy(georefimg_int32* in);
void georefimg_float_destroy(georefimg_float* in);
void georefimg_double_destroy(georefimg_double* in);

// Routines that load values from file
georefimg_uint8* georefimg_uint8_load(char *filename);
georefimg_uint16* georefimg_uint16_load(char *filename);
georefimg_uint32* georefimg_uint32_load(char *filename);
georefimg_int32* georefimg_int32_load(char *filename);
georefimg_float* georefimg_float_load(char *filename);
georefimg_double* georefimg_double_load(char *filename);

// Routines that save values to file
void georefimg_uint8_save(char *filename, georefimg_uint8 *out);
void georefimg_uint16_save(char *filename, georefimg_uint16 *out);
void georefimg_uint32_save(char *filename, georefimg_uint32 *out);
void georefimg_int32_save(char *filename, georefimg_int32 *out);
void georefimg_float_save(char *filename, georefimg_float *out);
void georefimg_double_save(char *filename, georefimg_double *out);

// Routines that display values to file
void georefimg_uint8_print(georefimg_uint8 *in);
void georefimg_uint16_print(georefimg_uint16 *in);
void georefimg_uint32_print(georefimg_uint32 *in);
void georefimg_int32_print(georefimg_int32 *in);
void georefimg_float_print(georefimg_float *in);
void georefimg_double_print(georefimg_double *in);





