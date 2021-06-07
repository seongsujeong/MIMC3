/*
A main program to run single image matching.
*/

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _STDLIB_H_
#include <stdlib.h>
#endif

#ifndef _MATH_H_
#include <math.h>
#endif

#ifndef _MIMC2_MODULE_
#include "MIMC_module.h"
#endif

#ifndef _MIMC2_MISC_
#include "MIMC_misc.h"
#endif

#ifndef _GMA_
#include "GMA.h"
#endif

#ifndef _TIME_H_
#include <time.h>
#endif

#ifndef _OMP_H
#include <omp.h>
#endif

#ifndef _TIFFIO_
#include "tiffio.h"
#endif

//list of global variables
float dt;   //temporal baseline
int32_t num_dp; //number of multiple matching attempts
int32_t num_grid,dimx_vmap,dimy_vmap; //# of grids, dimension in mapx and mapy direction
param param_mimc2;  //parameters that controlls the software
GMA_float **kernel; //kernel for the convolution filter
int main(int argc, char* argv[])
{
    char string_ver[]="3.0.5";

    printf("\n\nMIMC version %s\n",string_ver);
    printf("(C)Seongsu Jeong, Byrd Polar and Climate Research Center, The Ohio State University\n");
    printf("Official MIMC website: iceflow.wordpress.com\n");
    printf("Please address all the related concerns through the website above.\n\n");

    /*
    input:
    - argv[1]: i0
    - argv[2]: i1
    - argv[3]: xyuvav
    - argv[4]: outout directory
    */

    
    unsigned int cnt,cntc,cnt1,cnt2;
    
    //GMA_uint16 *i0, *i1;
    GMA_float *i0, *i1, *i0c, *i1c;
    GMA_double *xyuvav;
    GMA_double *t0t1;
    GMA_int32 **uv_pivot;
    GMA_float *dp[32];
    GMA_float **vxyexyqual;
    //char gma_out[512];
    //param param_mimc2;
    //time_t t0=0,t1=0;
    double t0=0, t1=0;
    num_dp=32;
    
    int32_t ocw;
    int32_t cnt_grid;

    char timestampstr_i0[15];
    char timestampstr_i1[15];
    
    //outout filenames
    char filename_x[1024];
    char filename_y[1024];
    char filename_vx[1024];
    char filename_vy[1024];
    char filename_ex[1024];
    char filename_ey[1024];
    char filename_qual[1024];
    char filename_flagcp[1024];
    char filename_meta[1024];
    char filename_vmap[1024];
    char string_command[4096];
    


    getTimeStampStr(argv[1],timestampstr_i0);
    getTimeStampStr(argv[2],timestampstr_i1);    

    sprintf(filename_x,"%s/vmap_%s_%s_x.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    sprintf(filename_y,"%s/vmap_%s_%s_y.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    sprintf(filename_vx,"%s/vmap_%s_%s_vx.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    sprintf(filename_vy,"%s/vmap_%s_%s_vy.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    sprintf(filename_ex,"%s/vmap_%s_%s_ex.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    sprintf(filename_ey,"%s/vmap_%s_%s_ey.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    sprintf(filename_qual,"%s/vmap_%s_%s_qual.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    sprintf(filename_flagcp,"%s/vmap_%s_%s_flagcp.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    sprintf(filename_meta,"%s/vmap_%s_%s_meta.txt",argv[4],timestampstr_i0,timestampstr_i1);
    sprintf(filename_vmap,"%s/vmap_%s_%s.tar",argv[4],timestampstr_i0,timestampstr_i1);

    //check the existence of the vmap
    //TODO: Skip the processing if all result files exists

    dt=get_dt(timestampstr_i0,timestampstr_i1);
    
    printf("Input files:\n",argv[1]);
    printf("dt=%f days\n\n",dt);
    printf("Earlier image: %s\n",argv[1]);
    printf("Latter image: %s\n",argv[2]);
    printf("xyuvav matrix: %s\n",argv[3]);
    
    printf("output files:\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
           filename_x,filename_y,
           filename_vx,filename_vy,
           filename_ex,filename_ey,
           filename_qual,filename_flagcp,
           filename_meta,filename_vmap);

    sprintf(string_command,"ls -l %s",filename_vmap);
    int status_exixtence_vmap=system(string_command);
    if(status_exixtence_vmap==0)
    {
        printf("vmap already exists. Skipping\n");
        return -1;
    }
    
    //define the param
    printf("Definining params...");
    param_mimc2.vec_ocw[0]=7;
    param_mimc2.vec_ocw[1]=15;
    param_mimc2.vec_ocw[2]=30;
    param_mimc2.vec_ocw[3]=40;
    
    //Faster feature tracking for testing the portprocessing algorithm
    //param_mimc2.vec_ocw[0]=7;
    //param_mimc2.vec_ocw[1]=10;
    //param_mimc2.vec_ocw[2]=15;
    //param_mimc2.vec_ocw[3]=11;
    
    //Larger window size for highres image
    //param_mimc2.vec_ocw[0]=14;
    //param_mimc2.vec_ocw[1]=30;
    //param_mimc2.vec_ocw[2]=60;
    //param_mimc2.vec_ocw[3]=80;
    

    param_mimc2.AW_CRE=10.0;
    param_mimc2.AW_SF=1.8;
    param_mimc2.mpp=15; //meter per pixel
    param_mimc2.spacing_grid=300; //meter per pixel

    //param_mimc2.radius_neighbor=1500;
    //param_mimc2.radius_neighbor_dpf1_px=1000;
    //param_mimc2.radius_neighbor_ps=1500;
    param_mimc2.radius_neighbor=5.0;
    param_mimc2.radius_neighbor_dpf1=1000/300;
    param_mimc2.radius_neighbor_ps=5.0;
    //note change the unit of those radius values as spacing
    
    param_mimc2.num_cp_max=500;
    param_mimc2.num_cp_min=50;
    param_mimc2.ratio_cp=0.03;
    param_mimc2.thres_spd_cp=10;

    printf("successful\n");

    
    //designate the kernel
    printf("Definining Kernels...");
    kernel=(GMA_float**)malloc(sizeof(GMA_float*)*3);
    kernel[0]=GMA_float_create(1,3);    //horizontal gradient
    kernel[1]=GMA_float_create(3,1);    //vertical gradient
    kernel[2]=GMA_float_create(3,3);    //Laplacian
    
    kernel[0]->val[0][0]=-1;
    kernel[0]->val[0][1]=0;
    kernel[0]->val[0][2]=1;

    kernel[1]->val[0][0]=-1;
    kernel[1]->val[1][0]=0;
    kernel[1]->val[2][0]=1;

    kernel[2]->val[0][0]=-1.0/8;
    kernel[2]->val[0][1]=-1.0/8;
    kernel[2]->val[0][2]=-1.0/8;
    kernel[2]->val[1][0]=-1.0/8;
    kernel[2]->val[1][1]=1.0;
    kernel[2]->val[1][2]=-1.0/8;
    kernel[2]->val[2][0]=-1.0/8;
    kernel[2]->val[2][1]=-1.0/8;
    kernel[2]->val[2][2]=-1.0/8;

    printf("successful\n");
    
    //load xyuvav
    printf("Loading xyuvav...");
    //xyuvav=GMA_double_load("../MIMC2_C_GMA/xyuvav_rect_interp.GMA");
    xyuvav=GMA_double_load(argv[3]);
    num_grid=xyuvav->nrows;
    printf("successful\n");

    

    num_grid=xyuvav->nrows;
    //calculate dimx_vmap and dimy_vmap
    for(cnt_grid=1;cnt_grid<num_grid;cnt_grid++)
    {
        if((int)xyuvav->val[cnt_grid][2]==(int)xyuvav->val[0][2])
        {
            break;
        }
    }
    dimx_vmap=cnt_grid;
    dimy_vmap=num_grid/dimx_vmap;

    param_mimc2.mpp=(float)((xyuvav->val[1][0]-xyuvav->val[0][0])/(xyuvav->val[1][2]-xyuvav->val[0][2]));
    param_mimc2.spacing_grid=(float)(xyuvav->val[1][2]-xyuvav->val[0][2]);
    param_mimc2.meter_per_spacing=(float)(xyuvav->val[1][0]-xyuvav->val[0][0]);

    printf("MPP=%f, grid spacing=%f, meter per spacing=%fm\nDimension of the vmap: %d by %d (mapy / mapx)\n",param_mimc2.mpp,param_mimc2.spacing_grid,param_mimc2.meter_per_spacing,dimy_vmap,dimx_vmap);
    
    //load the images
    printf("Loading i0 and i1...");
    i0=GMA_float_load_tiff(argv[1]);
    i1=GMA_float_load_tiff(argv[2]);
    
    printf("successful\n");
    
    
    /*
    //Measure the CP offset of the image pair
    printf("Measuring CP offset - ");
    
    
    for(cnt_grid=0;cnt_grid<num_grid;cnt_grid++)
    {
        flag_cp->val[cnt_grid][0]=0;
    }
    offset_cp[0]=0;
    offset_cp[1]=1;
    
    if(get_offset_image(i0,i1,kernel,xyuvav,offset_cp,flag_cp)<0)
    {
        sprintf(string_command,"touch %s",filename_vmap);
        system(string_command);
        return -1;
    }
    
    sprintf(output_dp,"%s/flag_cp.gma",argv[4]);
    GMA_uint8_save(output_dp,flag_cp);
    */
    char output_dp[2048];
    int32_t offset_cp[2],offset_cp_reverse[2];
    offset_cp[0]=0;
    offset_cp[1]=1;
    offset_cp_reverse[0]=-offset_cp[0];
    offset_cp_reverse[1]=-offset_cp[1];
    GMA_uint8 *flag_cp=GMA_uint8_create(num_grid,1);
    sprintf(output_dp,"%s/FT_result/flag_cp.gma",argv[4]);
    flag_cp=GMA_uint8_load(output_dp);

    //Perform main feature tracking
    
    for(cnt=0;cnt<32;cnt++)
    {
        sprintf(output_dp,"%s/FT_result/dp_%02d.gma",argv[4],cnt);
        dp[cnt]=GMA_float_load(output_dp);
    }

    /*
    printf("Initiating matching\n");
    char output_dp[2048];
    for(cnt=0;cnt<4;cnt++)
    {
        ocw=param_mimc2.vec_ocw[cnt];
        uv_pivot=get_uv_pivot(xyuvav, dt, param_mimc2, ocw, i1);
        printf("Feature tracking - ocw=%d, Original forward - ",ocw);
        t0=omp_get_wtime();
        dp[cnt*2]=matching_ncc_dlc_2(i0,i1,xyuvav,offset_cp,uv_pivot,ocw,param_mimc2.AW_CRE,param_mimc2.AW_SF);
        t1=omp_get_wtime();
        printf("Elapsed time: %f\n",(t1-t0));
        
        sprintf(output_dp,"%s/dp_%d.gma",argv[4],cnt*2);
        GMA_float_save(output_dp,dp[cnt*2]);

        
        
        //reverse uv_pivot
        //printf("Reversing the pivot...");
        for(cnt1=0;cnt1<xyuvav->nrows;cnt1++)
        {
            for(cnt2=0;cnt2<uv_pivot[cnt1]->nrows;cnt2++)
            {
                uv_pivot[cnt1]->val[cnt2][0]=-(uv_pivot[cnt1]->val[cnt2][0]);
                uv_pivot[cnt1]->val[cnt2][1]=-(uv_pivot[cnt1]->val[cnt2][1]);
            }
        }
        //printf("Successful\n");
        
        printf("Feature tracking - ocw=%d, Swapped forward - ",ocw);
        t0=omp_get_wtime();
        dp[cnt*2+1]=matching_ncc_dlc_2(i1,i0,xyuvav,offset_cp_reverse,uv_pivot,ocw,param_mimc2.AW_CRE,param_mimc2.AW_SF);
        t1=omp_get_wtime();
        printf("Elapsed time: %f\n",(t1-t0));
        
        
        //reverse dp
        for(cnt1=0;cnt1<xyuvav->nrows;cnt1++)
        {
            dp[cnt*2+1]->val[cnt1][0]=-(dp[cnt*2+1]->val[cnt1][0]);
            dp[cnt*2+1]->val[cnt1][1]=-(dp[cnt*2+1]->val[cnt1][1]);
        }

        sprintf(output_dp,"%s/dp_%d.gma",argv[4],cnt*2+1);
        GMA_float_save(output_dp,dp[cnt*2+1]);

        
        for(cnt1=0;cnt1<xyuvav->nrows;cnt1++)
        {
            GMA_int32_destroy(uv_pivot[cnt1]);
        }
        
    }


    //Perform multiple matching for the filtered images
    i0c=GMA_float_create(i0->nrows,i0->ncols);
    i1c=GMA_float_create(i1->nrows,i1->ncols);
    for(cnt=0;cnt<3;cnt++)  //loop for the kernel for the spatial filter
    {
        printf("Generating filtered image (%d of %d)...",cnt,3);
        GMA_float_conv2(i0,kernel[cnt],i0c);
        GMA_float_conv2(i1,kernel[cnt],i1c);
        printf("successful\n");
        for(cntc=0;cntc<4;cntc++)   //loop for the reference chip size
        {
            ocw=param_mimc2.vec_ocw[cntc];
            uv_pivot=get_uv_pivot(xyuvav, dt, param_mimc2, ocw, i1);
            printf("Feature tracking - ocw=%d, Original forward - ",ocw);
            t0=omp_get_wtime();
            dp[cnt*8+cntc*2+8]=matching_ncc_dlc_2(i0c,i1c,xyuvav,offset_cp,uv_pivot,ocw,param_mimc2.AW_CRE,param_mimc2.AW_SF);
            t1=omp_get_wtime();
            printf("Elapsed time: %f\n",(t1-t0));
            
            sprintf(output_dp,"%s/dp_%d.gma",argv[4],cnt*8+cntc*2+8);
            GMA_float_save(output_dp,dp[cnt*8+cntc*2+8]);

            //reverse uv_pivot
            for(cnt1=0;cnt1<xyuvav->nrows;cnt1++)
            {
                for(cnt2=0;cnt2<uv_pivot[cnt1]->nrows;cnt2++)
                {
                    uv_pivot[cnt1]->val[cnt2][0]=-(uv_pivot[cnt1]->val[cnt2][0]);
                    uv_pivot[cnt1]->val[cnt2][1]=-(uv_pivot[cnt1]->val[cnt2][1]);
                }
            }
            printf("Feature tracking - ocw=%d, Swapped forward - ",ocw);
            t0=omp_get_wtime();
            dp[cnt*8+cntc*2+1+8]=matching_ncc_dlc_2(i1c,i0c,xyuvav,offset_cp_reverse,uv_pivot,ocw,param_mimc2.AW_CRE,param_mimc2.AW_SF);
            t1=omp_get_wtime();
            printf("Elapsed time: %f\n",(t1-t0));
            
            
            //reverse dp
            for(cnt1=0;cnt1<xyuvav->nrows;cnt1++)
            {
                dp[cnt*8+cntc*2+1+8]->val[cnt1][0]=-(dp[cnt*8+cntc*2+1+8]->val[cnt1][0]);
                dp[cnt*8+cntc*2+1+8]->val[cnt1][1]=-(dp[cnt*8+cntc*2+1+8]->val[cnt1][1]);
            }
            
            sprintf(output_dp,"%s/dp_%d.gma",argv[4],cnt*8+cntc*2+1+8);
            GMA_float_save(output_dp,dp[cnt*8+cntc*2+1+8]);


            for(cnt1=0;cnt1<xyuvav->nrows;cnt1++)
            {
                GMA_int32_destroy(uv_pivot[cnt1]);
            }
            
        }
    }
    //Feature tracking complete. Proceeding to the postprocessing
    */


    printf("Initiating postprocessing\n");
    vxyexyqual=mimc2_postprocess(dp,xyuvav,dt);
    
    //correct translative offset in the FT result
    float sdu=0.0;
    float sdv=0.0;
    float du_cp=0.0;
    float dv_cp=0.0;
    float ducp,dvcp;
    int32_t num_cp=0;
    for(cnt1=0;cnt1<dimy_vmap;cnt1++)
    {
        for(cnt2=0;cnt2<dimx_vmap;cnt2++)
        {
            cnt_grid=cnt1*dimx_vmap+cnt2;
            ducp=vxyexyqual[0]->val[cnt1][cnt2];
            dvcp=vxyexyqual[1]->val[cnt1][cnt2];
            if(!isnan(ducp) && !isnan(dvcp))
            {
                sdu+=ducp;
                sdv+=dvcp;
                num_cp++;
            }
        }
    }
    du_cp=sdu/(float)num_cp;
    dv_cp=sdv/(float)num_cp;

    for(cnt1=0;cnt1<dimy_vmap;cnt1++)
    {
        for(cnt2=0;cnt2<dimx_vmap;cnt2++)
        {
            vxyexyqual[0]->val[cnt1][cnt2]=vxyexyqual[0]->val[cnt1][cnt2]-du_cp;
            vxyexyqual[1]->val[cnt1][cnt2]=vxyexyqual[1]->val[cnt1][cnt2]-dv_cp;
        }
    }



    //convert the image displacement to the movement in map coordinates
    float factor_image_to_map=param_mimc2.mpp/dt*365;
    for(cnt1=0;cnt1<dimy_vmap;cnt1++)
    {
        for(cnt2=0;cnt2<dimx_vmap;cnt2++)
        {
            vxyexyqual[0]->val[cnt1][cnt2]=vxyexyqual[0]->val[cnt1][cnt2]*factor_image_to_map;
            vxyexyqual[1]->val[cnt1][cnt2]=-vxyexyqual[1]->val[cnt1][cnt2]*factor_image_to_map;
            vxyexyqual[2]->val[cnt1][cnt2]=sqrt(vxyexyqual[2]->val[cnt1][cnt2])*factor_image_to_map;
            vxyexyqual[3]->val[cnt1][cnt2]=sqrt(vxyexyqual[3]->val[cnt1][cnt2])*factor_image_to_map;
        }
    }


    //extract x and y grid and save it to GMA
    GMA_double *grid_x=GMA_double_create(1,dimx_vmap);
    GMA_double *grid_y=GMA_double_create(1,dimy_vmap);
    for(cnt=0;cnt<dimx_vmap;cnt++)
    {
        grid_x->val[0][cnt]=xyuvav->val[cnt][0];
    }
    for(cnt=0;cnt<dimy_vmap;cnt++)
    {
        grid_y->val[0][cnt]=xyuvav->val[cnt*dimx_vmap][1];
    }

    printf("Saving the output\n");
    GMA_uint8 *flag_cp_2d=GMA_uint8_create(dimy_vmap,dimx_vmap);
    for(cnt=0;cnt<dimy_vmap;cnt++)
    {
        for(cnt1=0;cnt1<dimx_vmap;cnt1++)
        {
            flag_cp_2d->val[cnt][cnt1]=flag_cp->val[cnt*dimx_vmap+cnt1][0];
        }
    }

    
    GMA_double_save(filename_x,grid_x);
    GMA_double_save(filename_y,grid_y);
    GMA_float_save(filename_vx,vxyexyqual[0]);
    GMA_float_save(filename_vy,vxyexyqual[1]);
    GMA_float_save(filename_ex,vxyexyqual[2]);
    GMA_float_save(filename_ey,vxyexyqual[3]);
    GMA_float_save(filename_qual,vxyexyqual[4]);
    GMA_uint8_save(filename_flagcp,flag_cp_2d);
    
    //save the metadata
    FILE *fout;
	fout=fopen(filename_meta,"w");
    fprintf(fout,"MIMC_version=%s\n",string_ver);
    fprintf(fout,"name_i0=%s\n",argv[1]);
    fprintf(fout,"name_i1=%s\n",argv[2]);
    fprintf(fout,"cp_offset_int_u=%d\n",offset_cp[0]);
    fprintf(fout,"cp_offset_int_v=%d\n",offset_cp[1]);
    fprintf(fout,"cp_offset_subint_u=%f\n",du_cp);
    fprintf(fout,"cp_offset_subint_v=%f\n",dv_cp);
    fclose(fout);

    //sprintf(filename_x,"vmap_%s_%s_x.GMA",timestampstr_i0,timestampstr_i1);
    //sprintf(filename_y,"vmap_%s_%s_y.GMA",timestampstr_i0,timestampstr_i1);
    //sprintf(filename_vx,"vmap_%s_%s_vx.GMA",timestampstr_i0,timestampstr_i1);
    //sprintf(filename_vy,"vmap_%s_%s_vy.GMA",timestampstr_i0,timestampstr_i1);
    //sprintf(filename_ex,"vmap_%s_%s_ex.GMA",timestampstr_i0,timestampstr_i1);
    //sprintf(filename_ey,"vmap_%s_%s_ey.GMA",timestampstr_i0,timestampstr_i1);
    //sprintf(filename_qual,"vmap_%s_%s_qual.GMA",timestampstr_i0,timestampstr_i1);
    //sprintf(filename_flagcp,"vmap_%s_%s_flagcp.GMA",timestampstr_i0,timestampstr_i1);
    //sprintf(filename_meta,"vmap_%s_%s_meta.txt",timestampstr_i0,timestampstr_i1);
    //sprintf(filename_vmap,"vmap_%s_%s.tar",timestampstr_i0,timestampstr_i1);


    //                                 1  2  3  4  5  6  7  8  9 10    11
    sprintf(string_command,"tar -cvzf %s %s %s %s %s %s %s %s %s %s -C %s",
    filename_vmap,filename_x,filename_y,filename_vx,filename_vy,filename_ex,filename_ey,filename_qual,filename_flagcp,filename_meta,argv[4]);
  //1             2          3          4           5           6           7           8             9               10            11
    system(string_command);
    //printf("%s\n",string_command);

    //sprintf(filename_x,"%s/vmap_%s_%s_x.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    //sprintf(filename_y,"%s/vmap_%s_%s_y.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    //sprintf(filename_vx,"%s/vmap_%s_%s_vx.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    //sprintf(filename_vy,"%s/vmap_%s_%s_vy.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    //sprintf(filename_ex,"%s/vmap_%s_%s_ex.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    //sprintf(filename_ey,"%s/vmap_%s_%s_ey.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    //sprintf(filename_qual,"%s/vmap_%s_%s_qual.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    //sprintf(filename_flagcp,"%s/vmap_%s_%s_flagcp.GMA",argv[4],timestampstr_i0,timestampstr_i1);
    //sprintf(filename_meta,"%s/vmap_%s_%s_meta.txt",argv[4],timestampstr_i0,timestampstr_i1);
    //sprintf(filename_vmap,"%s/vmap_%s_%s.tar",argv[4],timestampstr_i0,timestampstr_i1);
    //                          1  2  3  4  5  6  7  8  9
    sprintf(string_command,"rm %s %s %s %s %s %s %s %s %s",
    filename_x,filename_y,filename_vx,filename_vy,filename_ex,filename_ey,filename_qual,filename_flagcp,filename_meta);
  //1           2         3           4           5           6           7             8               9
    system(string_command);
    //printf("%s\n",string_command);



    //finalize the program - destroy the allocated arrays

    //destroy the uv_pivot array
    printf("Deallocating memory...");
    printf("vxyexyqual and flag_cp...");
    GMA_uint8_destroy(flag_cp_2d);
    GMA_float_destroy(vxyexyqual[0]);
    GMA_float_destroy(vxyexyqual[1]);
    GMA_float_destroy(vxyexyqual[2]);
    GMA_float_destroy(vxyexyqual[3]);
    GMA_float_destroy(vxyexyqual[4]);
    free(vxyexyqual);
    
    GMA_uint8_destroy(flag_cp);
    
    printf("i0...");
    GMA_float_destroy(i0);
    GMA_float_destroy(i0c);
    printf("i1...");
    GMA_float_destroy(i1);
    GMA_float_destroy(i1c);
    printf("xyuvav...");
    GMA_double_destroy(xyuvav);
    printf("1st kernel...");
    GMA_float_destroy(kernel[0]);
    printf("2nd kernel...");
    GMA_float_destroy(kernel[1]);
    printf("3rd kernel...");
    GMA_float_destroy(kernel[2]);
    printf("kernel array...");
    free(kernel);
    
    printf("dp array...");
    for(cnt=0;cnt<32;cnt++)
    {
        GMA_float_destroy(dp[cnt]);
    }
    printf("Successful\n");

    
    
    printf("Processing completed\n");

    return 0;
}



