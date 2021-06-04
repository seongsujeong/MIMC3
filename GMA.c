#ifndef _GMA_
#include "GMA.h"
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

//Initialization routine
GMA_uint8* GMA_uint8_create(int32_t size_row, int32_t size_col)
{
    int32_t cnt;
    GMA_uint8 *out;
    out=(GMA_uint8*)malloc(sizeof(GMA_uint8));
	out->nrows=size_row;
	out->ncols=size_col;
	out->val=malloc(sizeof(uint8_t*)*size_row);
	for(cnt=0;cnt<size_row;cnt++) out->val[cnt]=malloc(sizeof(uint8_t)*size_col);
	return out;
}

GMA_uint16* GMA_uint16_create(int32_t size_row, int32_t size_col)
{
    int32_t cnt;
    GMA_uint16 *out;
    out=(GMA_uint16*)malloc(sizeof(GMA_uint16));
	out->nrows=size_row;
	out->ncols=size_col;
	out->val=malloc(sizeof(uint16_t*)*size_row);
	out->data=(uint16_t*)malloc(sizeof(uint16_t)*size_row*size_col);
	for(cnt=0;cnt<size_row;cnt++)
	{
		out->val[cnt]=&(out->data[cnt*size_col]);
	}
	return out;
}

GMA_uint32* GMA_uint32_create(int32_t size_row, int32_t size_col)
{
    int32_t cnt;
    GMA_uint32 *out;
    out=(GMA_uint32*)malloc(sizeof(GMA_uint32));
	out->nrows=size_row;
	out->ncols=size_col;
	out->val=malloc(sizeof(uint32_t*)*size_row);
	out->data=(uint32_t*)malloc(sizeof(uint32_t)*size_row*size_col);
	for(cnt=0;cnt<size_row;cnt++)
	{
		out->val[cnt]=&(out->data[cnt*size_col]);
	}
	return out;
}

GMA_int32* GMA_int32_create(int32_t size_row, int32_t size_col)
{
    int32_t cnt;
    GMA_int32 *out;
	out=(GMA_int32*)malloc(sizeof(GMA_int32));
	out->nrows=size_row;
	out->ncols=size_col;
	out->val=(int32_t**)malloc(sizeof(int32_t*)*size_row);
	out->data=(int32_t*)malloc(sizeof(int32_t)*size_row*size_col);
	for(cnt=0;cnt<size_row;cnt++)
	{
		out->val[cnt]=&(out->data[cnt*size_col]);
	}
	return out;

}

GMA_float* GMA_float_create(int32_t size_row, int32_t size_col)
{
	int32_t cnt;
    GMA_float *out;
	out=(GMA_float*)malloc(sizeof(GMA_float));
	out->nrows=size_row;
	out->ncols=size_col;
	out->val=(float**)malloc(sizeof(float*)*size_row);
	out->data=(float*)malloc(sizeof(float)*size_row*size_col);
	for(cnt=0;cnt<size_row;cnt++)
	{
		out->val[cnt]=&(out->data[cnt*size_col]);
	}
	return out;

}

GMA_double* GMA_double_create(int32_t size_row, int32_t size_col)
{
	int32_t cnt;
    GMA_double *out;
	out=(GMA_double*)malloc(sizeof(GMA_double));
	out->nrows=size_row;
	out->ncols=size_col;
	out->val=(double**)malloc(sizeof(double*)*size_row);
	//out->val[0]=(float*)malloc(sizeof(float)*size_row*size_col);
	out->data=(double*)malloc(sizeof(double)*size_row*size_col);
	for(cnt=0;cnt<size_row;cnt++)
	{
		//out->val[cnt]=out->val[0]+sizeof(float*)*size_col*cnt;
		out->val[cnt]=&(out->data[cnt*size_col]);
	}
	return out;
}


//De-allocation routine
void GMA_uint8_destroy(GMA_uint8* in)
{
    int32_t cnt;
    for(cnt=0;cnt<in->nrows;cnt++) free(in->val[cnt]);
    free(in);
}

void GMA_uint16_destroy(GMA_uint16* in)
{
    int32_t cnt;
    for(cnt=0;cnt<in->nrows;cnt++) free(in->val[cnt]);
    free(in);
}


void GMA_uint32_destroy(GMA_uint32* in)
{
    int32_t cnt;
    for(cnt=0;cnt<in->nrows;cnt++) free(in->val[cnt]);
    free(in);
}

void GMA_int32_destroy(GMA_int32* in)
{
	/*
    int32_t cnt;
    for(cnt=0;cnt<in->nrows;cnt++) free(in->val[cnt]);
    free(in);
	*/
	free(in->val);
	free(in->data);
	free(in);
}

void GMA_float_destroy(GMA_float* in)
{
	/*
	int32_t cnt;
    for(cnt=0;cnt<in->nrows;cnt++) free(in->val[cnt]);
    free(in->val);
	free(in);
	*/
	free(in->val);
	free(in->data);
	free(in);
	
}

void GMA_double_destroy(GMA_double* in)
{
	/*
    int32_t cnt;
    for(cnt=0;cnt<in->nrows;cnt++) free(in->val[cnt]);
    free(in);
	*/
	free(in->val);
	free(in->data);
	free(in);
}


//routines to load values from a file
GMA_uint8* GMA_uint8_load(char *filename)
{
    FILE *fin;
	int32_t cntr,size_row,size_col;
	GMA_uint8 *out;
	fin=fopen(filename,"r");
	fread(&size_row,sizeof(unsigned int),1,fin);
	fread(&size_col,sizeof(unsigned int),1,fin);
    out=GMA_uint8_create(size_row,size_col);
	for(cntr=0;cntr<size_row;cntr++) fread((void*)(out->val[cntr]),sizeof(uint8_t),size_col,fin);
	return out;
}

GMA_uint16* GMA_uint16_load(char *filename)
{
    FILE *fin;
	int32_t cntr,size_row,size_col;
	GMA_uint16 *out;
	fin=fopen(filename,"r");
	fread(&size_row,sizeof(unsigned int),1,fin);
	fread(&size_col,sizeof(unsigned int),1,fin);
    out=GMA_uint16_create(size_row,size_col);
	for(cntr=0;cntr<size_row;cntr++) fread((void*)(out->val[cntr]),sizeof(uint16_t),size_col,fin);
	return out;
}

GMA_uint32* GMA_uint32_load(char *filename)
{
    FILE *fin;
	int32_t cntr,size_row,size_col;
	GMA_uint32 *out;
	fin=fopen(filename,"r");
	fread(&size_row,sizeof(unsigned int),1,fin);
	fread(&size_col,sizeof(unsigned int),1,fin);
    out=GMA_uint32_create(size_row,size_col);
	for(cntr=0;cntr<size_row;cntr++) fread((void*)(out->val[cntr]),sizeof(uint32_t),size_col,fin);
	return out;
}

GMA_int32* GMA_int32_load(char *filename)
{
    FILE *fin;
	int32_t cntr,size_row,size_col;
	GMA_int32 *out;
	fin=fopen(filename,"r");
	fread(&size_row,sizeof(unsigned int),1,fin);
	fread(&size_col,sizeof(unsigned int),1,fin);
    out=GMA_int32_create(size_row,size_col);
	for(cntr=0;cntr<size_row;cntr++) fread((void*)(out->val[cntr]),sizeof(int32_t),size_col,fin);
	return out;
}

GMA_float* GMA_float_load(char *filename)
{
    FILE *fin;
	int32_t cntr,size_row,size_col;
	GMA_float *out;
	fin=fopen(filename,"r");
	fread(&size_row,sizeof(unsigned int),1,fin);
	fread(&size_col,sizeof(unsigned int),1,fin);
    out=GMA_float_create(size_row,size_col);
	for(cntr=0;cntr<size_row;cntr++) fread((void*)(out->val[cntr]),sizeof(float),size_col,fin);
	return out;
}

GMA_double* GMA_double_load(char *filename)
{
    FILE *fin;
	int32_t cntr,size_row,size_col;
	GMA_double *out;
	fin=fopen(filename,"r");
	fread(&size_row,sizeof(unsigned int),1,fin);
	fread(&size_col,sizeof(unsigned int),1,fin);
    out=GMA_double_create(size_row,size_col);
	for(cntr=0;cntr<size_row;cntr++) fread((void*)(out->val[cntr]),sizeof(double),size_col,fin);
	return out;
}

GMA_float* GMA_float_load_tiff(char *filename)
{
	GMA_float *out;
	//tiff loading routine from stack overflow
	TIFF* tif = TIFFOpen(filename, "r");
    if (tif)
	{
        uint32 imagelength,nsamples,imgdepth;
        tsize_t scanline;
        tdata_t buf;
        uint32 row;
        uint32 col;
		uint32 imageWidth;
		uint8 bytesperpixel;
		uint16_t dn;
		uint16_t *dnp_uint16;
		uint8_t *dnp_uint8;

		//uint16_t config;

        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
		TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
		TIFFGetField(tif, TIFFTAG_DATATYPE, &imgdepth);
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
        scanline = TIFFScanlineSize(tif);
        buf = _TIFFmalloc(scanline);
		bytesperpixel=(uint32)scanline/imageWidth;
		//bytesperpixel=(uint32)scanline/
		
		if(bytesperpixel==1)
		{
			dnp_uint8=(uint8_t*)buf;
		}
		else
		{
			dnp_uint16=(uint16_t*)buf;
		}


		printf("Loading TIFF - row=%d, col=%d, nsamples=%d, bytes per pixel=%d\n",imagelength,scanline,nsamples,bytesperpixel);
		out=GMA_float_create((int32_t)imagelength,(int32_t)imageWidth);
        
		if(bytesperpixel==1)
		{
			for (row = 0; row < imagelength; row++)
        	{
				TIFFReadScanline(tif, buf, row, 0);
				for(col=0;col<imageWidth;col++)
				{		
					out->data[row*imageWidth+col]=(float)dnp_uint8[col];
					
				}
        	}
		}
		else
		{
			for (row = 0; row < imagelength; row++)
        	{
				TIFFReadScanline(tif, buf, row, 0);
				for(col=0;col<imageWidth;col++)
				{		
					out->data[row*imageWidth+col]=(float)dnp_uint16[col];
				}
        	}
		}
        _TIFFfree(buf);
        TIFFClose(tif);
    }
	printf("Reference pixel DN=%f\n",out->val[23][34]);
	return out;
}

//Routines to save values to a file
void GMA_uint8_save(char *filename, GMA_uint8 *out)
{
    int32_t cntr,cntc;
	FILE *fout;
	fout=fopen(filename,"wb");
	fwrite(&out->nrows,sizeof(unsigned int),1,fout);  //write up the dimensions
	fwrite(&out->ncols,sizeof(unsigned int),1,fout);
	
	for(cntr=0;cntr<out->nrows;cntr++)
	{
		for(cntc=0;cntc<out->ncols;cntc++)
		{
			fwrite(&out->val[cntr][cntc],1,sizeof(uint8_t),fout);
		}
	}
	
	//fwrite(out->data,out->nrows*out->ncols,sizeof(uint8_t),fout);
	fclose(fout);
}



void GMA_uint16_save(char *filename, GMA_uint16 *out)
{
    int32_t cntr,cntc;
	FILE *fout;
	fout=fopen(filename,"wb");
	fwrite(&out->nrows,sizeof(unsigned int),1,fout);  //write up the dimensions
	fwrite(&out->ncols,sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->nrows;cntr++)
	{
		for(cntc=0;cntc<out->ncols;cntc++)
		{
			fwrite(&out->val[cntr][cntc],1,sizeof(uint16_t),fout);
		}
	}
	fclose(fout);
}

void GMA_uint32_save(char *filename, GMA_uint32 *out)
{
    int32_t cntr,cntc;
	FILE *fout;
	fout=fopen(filename,"wb");
	fwrite(&out->nrows,sizeof(unsigned int),1,fout);  //write up the dimensions
	fwrite(&out->ncols,sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->nrows;cntr++)
	{
		for(cntc=0;cntc<out->ncols;cntc++)
		{
			fwrite(&out->val[cntr][cntc],1,sizeof(uint32_t),fout);
		}
	}
	fclose(fout);
}

void GMA_int32_save(char *filename, GMA_int32 *out)
{
    int32_t cntr,cntc;
	FILE *fout;
	fout=fopen(filename,"wb");
	fwrite(&out->nrows,sizeof(unsigned int),1,fout);  //write up the dimensions
	fwrite(&out->ncols,sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->nrows;cntr++)
	{
		for(cntc=0;cntc<out->ncols;cntc++)
		{
			fwrite(&out->val[cntr][cntc],1,sizeof(int32_t),fout);
		}
	}
	fclose(fout);
}

void GMA_float_save(char *filename, GMA_float *out)
{
    int32_t cntr,cntc;
	FILE *fout;
	fout=fopen(filename,"wb");
	fwrite(&out->nrows,sizeof(unsigned int),1,fout);  //write up the dimensions
	fwrite(&out->ncols,sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->nrows;cntr++)
	{
		for(cntc=0;cntc<out->ncols;cntc++)
		{
			fwrite(&out->val[cntr][cntc],1,sizeof(float),fout);
		}
	}
	fclose(fout);
}

void GMA_double_save(char *filename, GMA_double *out)
{
    int32_t cntr,cntc;
	FILE *fout;
	fout=fopen(filename,"wb");
	fwrite(&out->nrows,sizeof(unsigned int),1,fout);  //write up the dimensions
	fwrite(&out->ncols,sizeof(unsigned int),1,fout);
	for(cntr=0;cntr<out->nrows;cntr++)
	{
		for(cntc=0;cntc<out->ncols;cntc++)
		{
			fwrite(&out->val[cntr][cntc],1,sizeof(double),fout);
		}
	}
	fclose(fout);
}

//routines to display the values
void GMA_uint8_print(GMA_uint8 *in){
    int32_t cntr,cntc;
    printf("GMA_uint8 - %d by %d\n",in->nrows,in->ncols);
    for(cntr=0;cntr<in->nrows;cntr++)
    {
        for(cntc=0;cntc<in->ncols;cntc++)
        {
            printf("%d\t",in->val[cntr][cntc]);
        }
        printf("\n");
    }
    printf("\n");
}

void GMA_uint16_print(GMA_uint16 *in){
    int32_t cntr,cntc;
    printf("GMA_uint16 - %d by %d\n",in->nrows,in->ncols);
    for(cntr=0;cntr<in->nrows;cntr++)
    {
        for(cntc=0;cntc<in->ncols;cntc++)
        {
            printf("%d\t",in->val[cntr][cntc]);
        }
        printf("\n");
    }
    printf("\n");
}

void GMA_uint32_print(GMA_uint32 *in){
    int32_t cntr,cntc;
    printf("GMA_uint32 - %d by %d\n",in->nrows,in->ncols);
    for(cntr=0;cntr<in->nrows;cntr++)
    {
        for(cntc=0;cntc<in->ncols;cntc++)
        {
            printf("%d\t",in->val[cntr][cntc]);
        }
        printf("\n");
    }
    printf("\n");
}

void GMA_int32_print(GMA_int32 *in){
    int32_t cntr,cntc;
    printf("GMA_int32 - %d by %d\n",in->nrows,in->ncols);
    for(cntr=0;cntr<in->nrows;cntr++)
    {
        for(cntc=0;cntc<in->ncols;cntc++)
        {
            printf("%d\t",in->val[cntr][cntc]);
        }
        printf("\n");
    }
    printf("\n");
}

void GMA_float_print(GMA_float *in){
    int32_t cntr,cntc;
    printf("GMA_float - %d by %d\n",in->nrows,in->ncols);
    for(cntr=0;cntr<in->nrows;cntr++)
    {
        for(cntc=0;cntc<in->ncols;cntc++)
        {
            //printf("%.2f\t",in->val[cntr][cntc]);
			printf("%f\t",in->val[cntr][cntc]);
        }
        printf("\n");
    }
    printf("\n");
}

void GMA_double_print(GMA_double *in){
    int32_t cntr,cntc;
    printf("GMA_double - %d by %d\n",in->nrows,in->ncols);
    for(cntr=0;cntr<in->nrows;cntr++)
    {
        for(cntc=0;cntc<in->ncols;cntc++)
        {
            printf("%f\t",in->val[cntr][cntc]);
        }
        printf("\n");
    }
    printf("\n");
}
