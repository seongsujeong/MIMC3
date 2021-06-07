/*
GMA.h
- Defines the matrix arrays in different data types
- Defines the basic I/O routines for each defined GMA structs

Change log:
- 03/09/2017 : Struct definitions completed
- 03/10/2017 : Basic I/O definition completed

Code by Seongsu Jeong (jeong.134@osu.edu)
Byrd Polar and Climate Research Center, The Ohio State University

*/

#ifndef _STDINT_H_
#include <stdint.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifndef _TIFFIO_
#include "tiffio.h"
#endif

#define _GMA_

typedef struct
{
    float real;
    float imag;
} cpxf; //complex float

typedef struct
{
    double real;
    double imag;
} cpxd; //complex double



typedef struct
{
    int32_t ncols;
    int32_t nrows;
    uint8_t **val;
    uint8_t *data;
} GMA_uint8;

typedef struct
{
    int32_t ncols;
    int32_t nrows;
    uint16_t **val;
    uint16_t *data;
} GMA_uint16;

typedef struct
{
    int32_t ncols;
    int32_t nrows;
    uint32_t **val;
    uint32_t *data;
} GMA_uint32;


typedef struct
{
    int32_t ncols;
    int32_t nrows;
    int32_t **val;
    int32_t *data;
} GMA_int32;

typedef struct
{
    int32_t ncols;
    int32_t nrows;
    float **val;
    float *data;
    
} GMA_float;

typedef struct
{
    int32_t ncols;
    int32_t nrows;
    double **val;
    double *data;
} GMA_double;

typedef struct
{
    int32_t ncols;
    int32_t nrows;
    cpxf **val;
    cpxf *data;
    
} GMA_cpxf;

typedef struct
{
    int32_t ncols;
    int32_t nrows;
    cpxd **val;
    cpxd *data;
    
} GMA_cpxd;



// routines that create and destroy the structs
GMA_uint8* GMA_uint8_create(int32_t size_row, int32_t size_col);
GMA_uint16* GMA_uint16_create(int32_t size_row, int32_t size_col);
GMA_uint32* GMA_uint32_create(int32_t size_row, int32_t size_col);
GMA_int32* GMA_int32_create(int32_t size_row, int32_t size_col);
GMA_float* GMA_float_create(int32_t size_row, int32_t size_col);
GMA_double* GMA_double_create(int32_t size_row, int32_t size_col);
GMA_double* GMA_double_create(int32_t size_row, int32_t size_col);

void GMA_uint8_destroy(GMA_uint8* in);
void GMA_uint16_destroy(GMA_uint16* in);
void GMA_uint32_destroy(GMA_uint32* in);
void GMA_int32_destroy(GMA_int32* in);
void GMA_float_destroy(GMA_float* in);
void GMA_double_destroy(GMA_double* in);

//routines that loads and saves the GMA
GMA_uint8* GMA_uint8_load(char *filename);
GMA_uint16* GMA_uint16_load(char *filename);
GMA_uint32* GMA_uint32_load(char *filename);
GMA_int32* GMA_int32_load(char *filename);
GMA_float* GMA_float_load(char *filename);
GMA_double* GMA_double_load(char *filename);

//routines that load pixel values from tiff image file
GMA_float* GMA_float_load_tiff(char *filename);

void GMA_uint8_save(char *filename, GMA_uint8 *out);
void GMA_uint16_save(char *filename, GMA_uint16 *out);
void GMA_uint32_save(char *filename, GMA_uint32 *out);
void GMA_int32_save(char *filename, GMA_int32 *out);
void GMA_float_save(char *filename, GMA_float *out);
void GMA_double_save(char *filename, GMA_double *out);

//routines to display the values
void GMA_uint8_print(GMA_uint8 *in);
void GMA_uint16_print(GMA_uint16 *in);
void GMA_uint32_print(GMA_uint32 *in);
void GMA_int32_print(GMA_int32 *in);
void GMA_float_print(GMA_float *in);
void GMA_double_print(GMA_double *in);




