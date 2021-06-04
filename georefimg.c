#ifndef _GEOREFIMG_
#include "georefimg.h"
#endif

#ifndef _GMA_
#include "GMA.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif



//Intialization routines
georefimg_uint8* georefimg_uint8_create(unsigned int size_row, unsigned int size_col)
{
    georefimg_uint8 *out;
    out=malloc(sizeof(georefimg_uint8));
	out->x=malloc(sizeof(double)*size_col);
	out->y=malloc(sizeof(double)*size_row);
	out->z=GMA_uint8_create(size_row,size_col);
    return out;
}

georefimg_uint16* georefimg_uint16_create(unsigned int size_row, unsigned int size_col)
{
    georefimg_uint16 *out;
    out=malloc(sizeof(georefimg_uint16));
	out->x=malloc(sizeof(double)*size_col);
	out->y=malloc(sizeof(double)*size_row);
	out->z=GMA_uint16_create(size_row,size_col);
    return out;
}

georefimg_uint32* georefimg_uint32_create(unsigned int size_row, unsigned int size_col)
{
    georefimg_uint32 *out;
    out=malloc(sizeof(georefimg_uint32));
	out->x=malloc(sizeof(double)*size_col);
	out->y=malloc(sizeof(double)*size_row);
	out->z=GMA_uint32_create(size_row,size_col);
    return out;
}

georefimg_int32* georefimg_int32_create(unsigned int size_row, unsigned int size_col)
{
    georefimg_int32 *out;
    out=malloc(sizeof(georefimg_int32));
	out->x=malloc(sizeof(double)*size_col);
	out->y=malloc(sizeof(double)*size_row);
	out->z=GMA_int32_create(size_row,size_col);
    return out;
}

georefimg_float* georefimg_float_create(unsigned int size_row, unsigned int size_col)
{
    georefimg_float *out;
    out=malloc(sizeof(georefimg_float));
	out->x=malloc(sizeof(double)*size_col);
	out->y=malloc(sizeof(double)*size_row);
	out->z=GMA_float_create(size_row,size_col);
    return out;
}

georefimg_double* georefimg_double_create(unsigned int size_row, unsigned int size_col)
{
    georefimg_double *out;
    out=malloc(sizeof(georefimg_double));
	out->x=malloc(sizeof(double)*size_col);
	out->y=malloc(sizeof(double)*size_row);
	out->z=GMA_double_create(size_row,size_col);
    return out;
}


// De-allocation routines
void georefimg_uint8_destroy(georefimg_uint8* in){
    free(in->x);
    free(in->y);
    GMA_uint8_destroy(in->z);
    free(in);
}

void georefimg_uint16_destroy(georefimg_uint16* in){
    free(in->x);
    free(in->y);
    GMA_uint16_destroy(in->z);
    free(in);
}

void georefimg_uint32_destroy(georefimg_uint32* in){
    free(in->x);
    free(in->y);
    GMA_uint32_destroy(in->z);
    free(in);
}

void georefimg_int32_destroy(georefimg_int32* in){
    free(in->x);
    free(in->y);
    GMA_int32_destroy(in->z);
    free(in);
}

void georefimg_float_destroy(georefimg_float* in){
    free(in->x);
    free(in->y);
    GMA_float_destroy(in->z);
    free(in);
}

void georefimg_double_destroy(georefimg_double* in){
    free(in->x);
    free(in->y);
    GMA_double_destroy(in->z);
    free(in);
}


// Routines that load values from file
georefimg_uint8* georefimg_uint8_load(char *filename){
    unsigned int cntr;
	unsigned int size_row,size_col;
	georefimg_uint8 *out;
	FILE *fin;
	fin=fopen(filename,"r+");
	fread(&size_row,sizeof(unsigned int),1,fin);//read the size of the array
	fread(&size_col,sizeof(unsigned int),1,fin);
	out=georefimg_uint8_create(size_row,size_col);//allocate the space for the output GEOIMG_FLOAT
	for(cntr=0;cntr<size_row;cntr++) fread(out->z->val[cntr],sizeof(uint8_t),size_col,fin);
	fread(out->y,sizeof(double),size_row,fin);//read the grid data
	fread(out->x,sizeof(double),size_col,fin);
	fclose(fin);
	return out;
}

georefimg_uint16* georefimg_uint16_load(char *filename){
    unsigned int cntr;
	unsigned int size_row,size_col;
	georefimg_uint16 *out;
	FILE *fin;
	fin=fopen(filename,"r+");
	fread(&size_row,sizeof(unsigned int),1,fin);//read the size of the array
	fread(&size_col,sizeof(unsigned int),1,fin);
	out=georefimg_uint16_create(size_row,size_col);//allocate the space for the output GEOIMG_FLOAT
	for(cntr=0;cntr<size_row;cntr++) fread(out->z->val[cntr],sizeof(uint16_t),size_col,fin);
	fread(out->y,sizeof(double),size_row,fin);//read the grid data
	fread(out->x,sizeof(double),size_col,fin);
	fclose(fin);
	return out;
}

georefimg_uint32* georefimg_uint32_load(char *filename){
    unsigned int cntr;
	unsigned int size_row,size_col;
	georefimg_uint32 *out;
	FILE *fin;
	fin=fopen(filename,"r+");
	fread(&size_row,sizeof(unsigned int),1,fin);//read the size of the array
	fread(&size_col,sizeof(unsigned int),1,fin);
	out=georefimg_uint32_create(size_row,size_col);//allocate the space for the output GEOIMG_FLOAT
	for(cntr=0;cntr<size_row;cntr++) fread(out->z->val[cntr],sizeof(uint32_t),size_col,fin);
	fread(out->y,sizeof(double),size_row,fin);//read the grid data
	fread(out->x,sizeof(double),size_col,fin);
	fclose(fin);
	return out;
}

georefimg_int32* georefimg_int32_load(char *filename){
    unsigned int cntr;
	unsigned int size_row,size_col;
	georefimg_int32 *out;
	FILE *fin;
	fin=fopen(filename,"r+");
	fread(&size_row,sizeof(unsigned int),1,fin);//read the size of the array
	fread(&size_col,sizeof(unsigned int),1,fin);
	out=georefimg_int32_create(size_row,size_col);//allocate the space for the output GEOIMG_FLOAT
	for(cntr=0;cntr<size_row;cntr++) fread(out->z->val[cntr],sizeof(int32_t),size_col,fin);
	fread(out->y,sizeof(double),size_row,fin);//read the grid data
	fread(out->x,sizeof(double),size_col,fin);
	fclose(fin);
	return out;
}

georefimg_float* georefimg_float_load(char *filename){
    unsigned int cntr;
	unsigned int size_row,size_col;
	georefimg_float *out;
	FILE *fin;
	fin=fopen(filename,"r+");
	fread(&size_row,sizeof(unsigned int),1,fin);//read the size of the array
	fread(&size_col,sizeof(unsigned int),1,fin);
	out=georefimg_float_create(size_row,size_col);//allocate the space for the output GEOIMG_FLOAT
	for(cntr=0;cntr<size_row;cntr++) fread(out->z->val[cntr],sizeof(float),size_col,fin);
	fread(out->y,sizeof(double),size_row,fin);//read the grid data
	fread(out->x,sizeof(double),size_col,fin);
	fclose(fin);
	return out;
}

georefimg_double* georefimg_double_load(char *filename){
    unsigned int cntr;
	unsigned int size_row,size_col;
	georefimg_double *out;
	FILE *fin;
	fin=fopen(filename,"r+");
	fread(&size_row,sizeof(unsigned int),1,fin);//read the size of the array
	fread(&size_col,sizeof(unsigned int),1,fin);
	out=georefimg_double_create(size_row,size_col);//allocate the space for the output GEOIMG_FLOAT
	for(cntr=0;cntr<size_row;cntr++) fread(out->z->val[cntr],sizeof(double),size_col,fin);
	fread(out->y,sizeof(double),size_row,fin);//read the grid data
	fread(out->x,sizeof(double),size_col,fin);
	fclose(fin);
	return out;
}


// Routines that save values to file
void georefimg_uint8_save(char *filename, georefimg_uint8 *out)
{
	unsigned int cntr; 
	FILE *fout;
	fout=fopen(filename,"w+");
	fwrite(&(out->z->nrows),sizeof(unsigned int),1,fout); //write up the array dimension
	fwrite(&(out->z->ncols),sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->z->nrows;cntr++) fwrite(out->z->val[cntr],sizeof(uint8_t),out->z->ncols,fout); //write up the main data
	fwrite(out->y,sizeof(double),out->z->nrows,fout);//write up the grid data
	fwrite(out->x,sizeof(double),out->z->ncols,fout);
	fclose(fout);
}

void georefimg_uint16_save(char *filename, georefimg_uint16 *out)
{
	unsigned int cntr; 
	FILE *fout;
	fout=fopen(filename,"w+");
	fwrite(&(out->z->nrows),sizeof(unsigned int),1,fout); //write up the array dimension
	fwrite(&(out->z->ncols),sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->z->nrows;cntr++) fwrite(out->z->val[cntr],sizeof(uint16_t),out->z->ncols,fout); //write up the main data
	fwrite(out->y,sizeof(double),out->z->nrows,fout);//write up the grid data
	fwrite(out->x,sizeof(double),out->z->ncols,fout);
	fclose(fout);
}

void georefimg_uint32_save(char *filename, georefimg_uint32 *out)
{
	unsigned int cntr; 
	FILE *fout;
	fout=fopen(filename,"w+");
	fwrite(&(out->z->nrows),sizeof(unsigned int),1,fout); //write up the array dimension
	fwrite(&(out->z->ncols),sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->z->nrows;cntr++) fwrite(out->z->val[cntr],sizeof(uint32_t),out->z->ncols,fout); //write up the main data
	fwrite(out->y,sizeof(double),out->z->nrows,fout);//write up the grid data
	fwrite(out->x,sizeof(double),out->z->ncols,fout);
	fclose(fout);
}

void georefimg_int32_save(char *filename, georefimg_int32 *out)
{
	unsigned int cntr; 
	FILE *fout;
	fout=fopen(filename,"w+");
	fwrite(&(out->z->nrows),sizeof(unsigned int),1,fout); //write up the array dimension
	fwrite(&(out->z->ncols),sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->z->nrows;cntr++) fwrite(out->z->val[cntr],sizeof(int32_t),out->z->ncols,fout); //write up the main data
	fwrite(out->y,sizeof(double),out->z->nrows,fout);//write up the grid data
	fwrite(out->x,sizeof(double),out->z->ncols,fout);
	fclose(fout);
}

void georefimg_float_save(char *filename, georefimg_float *out)
{
	int32_t cntr; 
	FILE *fout;
	fout=fopen(filename,"w+");
	fwrite(&(out->z->nrows),sizeof(unsigned int),1,fout); //write up the array dimension
	fwrite(&(out->z->ncols),sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->z->nrows;cntr++) fwrite(out->z->val[cntr],sizeof(float),out->z->ncols,fout); //write up the main data
	fwrite(out->y,sizeof(double),out->z->nrows,fout);//write up the grid data
	fwrite(out->x,sizeof(double),out->z->ncols,fout);
	fclose(fout);
}

void georefimg_double_save(char *filename, georefimg_double *out)
{
	unsigned int cntr; 
	FILE *fout;
	fout=fopen(filename,"w+");
	fwrite(&(out->z->nrows),sizeof(unsigned int),1,fout); //write up the array dimension
	fwrite(&(out->z->ncols),sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->z->nrows;cntr++) fwrite(out->z->val[cntr],sizeof(double),out->z->ncols,fout); //write up the main data
	fwrite(out->y,sizeof(double),out->z->nrows,fout);//write up the grid data
	fwrite(out->x,sizeof(double),out->z->ncols,fout);
	fclose(fout);
}