#ifndef _MIMC2_MODULE_
#include "MIMC_module.h"
#endif

#ifndef _GMA_
#include "GMA.h"
#endif

#ifndef _MATH_H_
#include <math.h>
#endif

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _TIME_H_
#include <time.h>
#endif

#define MIN_DN 0.0000000001

#ifndef _OMP_H
#include <omp.h>
#endif


extern float dt;   //temporal baseline
extern int32_t num_dp; //number of multiple matching attempts
extern int32_t num_grid,dimx_vmap,dimy_vmap; //# of grids, dimension in mapx and mapy direction
extern param param_mimc2;  //parameters that controlls the software
extern GMA_float **kernel;
int get_offset_image(GMA_float *i0,GMA_float *i1,GMA_float **kernel, GMA_double *xyuvav, int32_t *offset, GMA_uint8 *flag_cp)
{
    //uint8_t isMeasurementSuccessful=1;
    int32_t uvoi[2];
    int32_t cnt,cntrow,cntrow1,cnt1,cnt2,cnt3,cntu,cntv;
    int32_t num_cp,num_cp_candidate;
    uint8_t *flag_cp_candidate=(uint8_t *)malloc(sizeof(uint8_t)*xyuvav->nrows);
    GMA_int32 *uv_pivot_cp;
    GMA_double *xyuvav_cp,*xyuvav_cp_sub;
    GMA_float *dp_cp[16], *dp_cp_sub[16];
    GMA_float *duv_cp;
    GMA_float *stack_dp_dp;
    GMA_float **imgchip_i0;
    GMA_float **imgchip_i1;
    int32_t *array_index_segment;
    int32_t ocw_chip=param_mimc2.vec_ocw[2]+param_mimc2.AW_CRE+2;
    int32_t num_dp_original=num_dp;
    //int32_t num_cp_least=50;
    char flag_sucessful_ft_cp=0;
    num_dp=16;//temporary value; don't forget to revert this value to the original
    GMA_float **mvn_dp_cp_sub;
    
    //determine the number of control points
    //printf("%f, %d\n",xyuvav->nrows*param_mimc2.ratio_cp, param_mimc2.num_cp_max);
    if(xyuvav->nrows*param_mimc2.ratio_cp > param_mimc2.num_cp_max)
    {
        num_cp=param_mimc2.num_cp_max;
    }
    else
    {
        num_cp=(int32_t)(xyuvav->nrows*param_mimc2.ratio_cp);
    }

    printf("Threshold for the # of CP=%d\n",num_cp);

    //determine the candidate for control points
    printf("Finding candidates for the control points...");
    num_cp_candidate=0;
    for(cntrow=0;cntrow<xyuvav->nrows;cntrow++)
    {
        float spd_apv_sq=xyuvav->val[cntrow][4]*xyuvav->val[cntrow][4]+xyuvav->val[cntrow][5]*xyuvav->val[cntrow][5];
        if(spd_apv_sq < param_mimc2.thres_spd_cp*param_mimc2.thres_spd_cp)
        {
            flag_cp_candidate[cntrow]=1;
            num_cp_candidate++;
        }
        else
        {
            flag_cp_candidate[cntrow]=0;
        }
    }

    //exclude the CP candidate if the corresponding image chip is not valid
    GMA_float *sarea_i0=GMA_float_create(param_mimc2.vec_ocw[2]*2+1,param_mimc2.vec_ocw[2]*2+1);
    GMA_float *sarea_i1=GMA_float_create(param_mimc2.vec_ocw[2]*2+1,param_mimc2.vec_ocw[2]*2+1);
    int32_t sum_invalid_i0;
    int32_t sum_invalid_i1;
    int32_t thres_numpx=(param_mimc2.vec_ocw[2]*2+1)*(param_mimc2.vec_ocw[2]*2+1)/2;
    for(cntrow=0;cntrow<xyuvav->nrows;cntrow++)
    {
        if(flag_cp_candidate[cntrow])
        {
            sum_invalid_i0=0;
            sum_invalid_i1=0;
            uvoi[0]=(int32_t)xyuvav->val[cntrow][2];
            uvoi[1]=(int32_t)xyuvav->val[cntrow][3];
            for(cnt1=-param_mimc2.vec_ocw[2];cnt1<=param_mimc2.vec_ocw[2];cnt1++)
            {
                for(cnt2=-param_mimc2.vec_ocw[2];cnt2<=param_mimc2.vec_ocw[2];cnt2++)
                {
                    //if(!isnan(i0->val[cnt1][cnt2]))
                    if(i0->val[cnt1+uvoi[1]][cnt2+uvoi[0]]<0.00001)
                    {
                        sum_invalid_i0++;
                    }

                    if(i1->val[cnt1+uvoi[1]][cnt2+uvoi[0]]<0.00001)
                    {
                        sum_invalid_i1++;
                    }
                }
            }
            if(sum_invalid_i0>thres_numpx || sum_invalid_i0>thres_numpx)
            {
                flag_cp_candidate[cntrow]=0;
                num_cp_candidate--;
            }
        }
    }
    GMA_float_destroy(sarea_i0);
    GMA_float_destroy(sarea_i1);
    
    //dump the flag_cp_candidate for the debugging purpose
    //GMA_uint8 *fcc=GMA_uint8_create(xyuvav->nrows,1);
    //fcc->data=flag_cp_candidate;
    //GMA_uint8_save("../MIMC2_C_GMA/flag_cp_candidate.GMA",fcc);

    printf("Completed: %d grids\n",num_cp_candidate);
    
    if(num_cp_candidate<param_mimc2.num_cp_min)
    {
        printf("Failed to extract sufficient amount of CP candidates (%d<%d)\n",num_cp_candidate,param_mimc2.num_cp_min);
        return -1;
    }
    
    if(num_cp>num_cp_candidate)
    {
        num_cp=(int32_t)((float)num_cp_candidate*0.75);
        printf("Not too much CP candidates(%d) -  #CP threshold was adjusted to %d\n",num_cp_candidate,num_cp);
    }
    
    //extract the CP grids from xyuvav
    cntrow1=0;
    xyuvav_cp=GMA_double_create(num_cp_candidate,7);
    duv_cp=GMA_float_create(num_cp_candidate,2);
    for(cntrow=0;cntrow<xyuvav->nrows;cntrow++)
    {
        if(flag_cp_candidate[cntrow])
        {
            xyuvav_cp->val[cntrow1][0]=xyuvav->val[cntrow][0];
            xyuvav_cp->val[cntrow1][1]=xyuvav->val[cntrow][1];
            xyuvav_cp->val[cntrow1][2]=xyuvav->val[cntrow][2];
            xyuvav_cp->val[cntrow1][3]=xyuvav->val[cntrow][3];
            xyuvav_cp->val[cntrow1][4]=xyuvav->val[cntrow][4];
            xyuvav_cp->val[cntrow1][5]=xyuvav->val[cntrow][5];
            xyuvav_cp->val[cntrow1][6]=(double)cntrow;  //row index of the grid chosen for the CP
            cntrow1++;
        }
    }

    //build up the uv_pivot for the extracted xyuvav_cp
    //NOTE: All CP measurement make use of the common uv_pivot
    int32_t num_pivot=(param_mimc2.AW_CRE*2+1)*(param_mimc2.AW_CRE*2+1);
    uv_pivot_cp=GMA_int32_create(num_pivot,2);
    int32_t cnt_pivot=0;
    for(cnt1=-param_mimc2.AW_CRE;cnt1<=param_mimc2.AW_CRE;cnt1++)
    {
        for(cnt2=-param_mimc2.AW_CRE;cnt2<=param_mimc2.AW_CRE;cnt2++)
        {
            uv_pivot_cp->val[cnt_pivot][0]=cnt1;
            uv_pivot_cp->val[cnt_pivot][1]=cnt2;
            cnt_pivot++;
        }
    }

    free(flag_cp_candidate);

    //randomly mix the order of the rows in xyuvav_cp
    GMA_double_randperm_row(xyuvav_cp);
    
    //determine the range of the segment
    //TODO: soft-code the numver of the segments, so that each segments has the grids whose mumbers are almost same as the CP threshold
    
    int32_t num_segment=(num_cp_candidate<param_mimc2.num_cp_min) ? 1 : num_cp_candidate/num_cp;
    
    printf("# segments=%d\n",num_segment);
    array_index_segment=(int32_t*)malloc(sizeof(int32_t)*(num_segment+1));
    array_index_segment[0]=0;
    
    for(cnt=1;cnt<=num_segment;cnt++)
    {
        array_index_segment[cnt]=(int32_t)(num_cp_candidate*((float)cnt/(float)num_segment));
    }

    //Go through fature tracking/dpf1 determination in each segment
    float sduv[2];
    sduv[0]=0.0;
    sduv[1]=0.0;
    int32_t num_cp_current=0;

    //    int32_t flag_is362273=-1;    //for debugging
    //    char filename_xyuvav_sub[1024];

    for(cnt=0;cnt<num_segment;cnt++)
    {   
        printf("CP measurement attempt: %d/%d - ",cnt+1,num_segment);
        int32_t num_grid_sub=array_index_segment[cnt+1]-array_index_segment[cnt];

        //allocate dp_cp_sub to store the feature tracking results of the segment
        for(cnt1=0;cnt1<16;cnt1++)  //NOTE: the "16" means the # of matching attempts per grid in CP measurement
        {
            dp_cp_sub[cnt1]=GMA_float_create(num_grid_sub,3);
        }
        
        //build up xyuvav_cp_sub
        xyuvav_cp_sub=GMA_double_create(num_grid_sub,7);
        
        //printf("Extracting xyuvav of the segment\n");
        for(cnt1=0;cnt1<num_grid_sub;cnt1++)
        {
            xyuvav_cp_sub->val[cnt1][0]=xyuvav_cp->val[array_index_segment[cnt]+cnt1][0];
            xyuvav_cp_sub->val[cnt1][1]=xyuvav_cp->val[array_index_segment[cnt]+cnt1][1];
            xyuvav_cp_sub->val[cnt1][2]=xyuvav_cp->val[array_index_segment[cnt]+cnt1][2];
            xyuvav_cp_sub->val[cnt1][3]=xyuvav_cp->val[array_index_segment[cnt]+cnt1][3];
            xyuvav_cp_sub->val[cnt1][4]=xyuvav_cp->val[array_index_segment[cnt]+cnt1][4];
            xyuvav_cp_sub->val[cnt1][5]=xyuvav_cp->val[array_index_segment[cnt]+cnt1][5];
            xyuvav_cp_sub->val[cnt1][6]=xyuvav_cp->val[array_index_segment[cnt]+cnt1][6];
            /*
            if((int32_t)(xyuvav_cp_sub->val[cnt1][6]+0.5)==362272)
            {
                printf("Heads up. Current segment contains the 362273rd segment\n");
                flag_is362273=cnt1;

            }
            sprintf(filename_xyuvav_sub,"../MIMC2_C_GMA/xyuvav_sub_%02d.GMA",cnt);
            GMA_double_save(filename_xyuvav_sub,xyuvav_cp_sub);
            */

        }
        //GMA_double_save("../MIMC2_C_GMA/xyuvav_cp_sub.GMA",xyuvav_cp_sub);
        //build up the image chip array
        imgchip_i0=(GMA_float**)malloc(sizeof(GMA_float*)*num_grid_sub);
        imgchip_i1=(GMA_float**)malloc(sizeof(GMA_float*)*num_grid_sub);
        
        //printf("Initiating CP image matching\n");
        for(cnt2=-1;cnt2<=2;cnt2++)  //index for kernel
        {
            if(cnt2<0)
            {//original image
                for(cnt1=0;cnt1<num_grid_sub;cnt1++)
                {
                    //just extract the image chips; no colvolution
                    int32_t uv_center[2];
                    uv_center[0]=(int32_t)xyuvav_cp_sub->val[cnt1][2];
                    uv_center[1]=(int32_t)xyuvav_cp_sub->val[cnt1][3];
                    imgchip_i0[cnt1]=GMA_float_create(ocw_chip*2+1,ocw_chip*2+1);
                    imgchip_i1[cnt1]=GMA_float_create(ocw_chip*2+1,ocw_chip*2+1);
                    for(cntv=-ocw_chip;cntv<=ocw_chip;cntv++)
                    {
                        for(cntu=-ocw_chip;cntu<=ocw_chip;cntu++)
                        {
                            imgchip_i0[cnt1]->val[cntv+ocw_chip][cntu+ocw_chip]=i0->val[uv_center[1]+cntv][uv_center[0]+cntu];
                            imgchip_i1[cnt1]->val[cntv+ocw_chip][cntu+ocw_chip]=i1->val[uv_center[1]+cntv][uv_center[0]+cntu];
                        }
                    }
                }
            
            }
            else
            {//convoluted image chips
                GMA_float *imgchip_i0_temp=GMA_float_create(ocw_chip*2+3,ocw_chip*2+3);
                GMA_float *imgchip_i1_temp=GMA_float_create(ocw_chip*2+3,ocw_chip*2+3);
                GMA_float *imgchip_i0_temp_out=GMA_float_create(ocw_chip*2+3,ocw_chip*2+3);
                GMA_float *imgchip_i1_temp_out=GMA_float_create(ocw_chip*2+3,ocw_chip*2+3);
                for(cnt1=0;cnt1<num_grid_sub;cnt1++)
                {
                    int32_t uv_center[2];
                    uv_center[0]=(int32_t)xyuvav_cp_sub->val[cnt1][2];
                    uv_center[1]=(int32_t)xyuvav_cp_sub->val[cnt1][3];
                    imgchip_i0[cnt1]=GMA_float_create(ocw_chip*2+1,ocw_chip*2+1);
                    imgchip_i1[cnt1]=GMA_float_create(ocw_chip*2+1,ocw_chip*2+1);
                    for(cntv=-ocw_chip-1;cntv<=ocw_chip+1;cntv++)
                    {
                        for(cntu=-ocw_chip-1;cntu<=ocw_chip+1;cntu++)
                        {
                            imgchip_i0_temp->val[ocw_chip+1+cntv][ocw_chip+1+cntu]=i0->val[uv_center[1]+cntv][uv_center[0]+cntu];
                            imgchip_i1_temp->val[ocw_chip+1+cntv][ocw_chip+1+cntu]=i1->val[uv_center[1]+cntv][uv_center[0]+cntu];
                        }
                    }
                    GMA_float_conv2(imgchip_i0_temp,kernel[cnt2],imgchip_i0_temp_out);
                    GMA_float_conv2(imgchip_i1_temp,kernel[cnt2],imgchip_i1_temp_out);
                    
                    //copy the convoluted image chip to array
                    for(cntv=-ocw_chip;cntv<=ocw_chip;cntv++)
                    {
                        for(cntu=-ocw_chip;cntu<=ocw_chip;cntu++)
                        {
                            imgchip_i0[cnt1]->val[cntv+ocw_chip][cntu+ocw_chip]=imgchip_i0_temp_out->val[cntv+ocw_chip+1][cntu+ocw_chip+1];
                            imgchip_i1[cnt1]->val[cntv+ocw_chip][cntu+ocw_chip]=imgchip_i1_temp_out->val[cntv+ocw_chip+1][cntu+ocw_chip+1];
                        }
                    }
                }
                GMA_float_destroy(imgchip_i0_temp);
                GMA_float_destroy(imgchip_i1_temp);
                GMA_float_destroy(imgchip_i0_temp_out);
                GMA_float_destroy(imgchip_i1_temp_out);
            }
            /*
            if(flag_is362273>=0)    //for debugging
            {
                sprintf(filename_xyuvav_sub,"../MIMC2_C_GMA/i0chip_%d.GMA",cnt2+1);
                GMA_float_save(filename_xyuvav_sub,imgchip_i0[flag_is362273]);
                sprintf(filename_xyuvav_sub,"../MIMC2_C_GMA/i1chip_%d.GMA",cnt2+1);
                GMA_float_save(filename_xyuvav_sub,imgchip_i1[flag_is362273]);
                
            }*/

            // image chip build complete. Proceeding to the image matching
            //TODO: parallelize this part
            float uvncc_cp[3];
            int32_t ocw_refchip;
            GMA_float *refchip;
            for(cnt3=1;cnt3<3;cnt3++)//counter for ocw
            {
                ocw_refchip=param_mimc2.vec_ocw[cnt3];
                
                //printf("ocw=%d, kernel=%d\n",cnt3,cnt2);
                //printf("%d...",cnt3+(cnt2-1)*2);
                //printf("...%d",(cnt3-1)*4+cnt2+1);
                #pragma omp parallel private(refchip,cntv,cntu,uvncc_cp) shared(dp_cp_sub)
                {
                    refchip=GMA_float_create(ocw_refchip*2+1,ocw_refchip*2+1);
                    #pragma omp for schedule(dynamic)   //TODO: consider parallelizing at upper level of the loop
                    for(cnt1=0;cnt1<num_grid_sub;cnt1++)
                    {
                        //original forward matching
                        //extract refchip from i0
                        for(cntv=-ocw_refchip;cntv<=ocw_refchip;cntv++)
                        {
                            for(cntu=-ocw_refchip;cntu<=ocw_refchip;cntu++)
                            {
                                refchip->val[ocw_refchip+cntv][ocw_refchip+cntu]=imgchip_i0[cnt1]->val[ocw_chip+cntv][ocw_chip+cntu];
                            }
                        }
                
                        //prepare for the sarea from i1 -> imgchip_i1[cnt1]
                    
                        //perform feature tracking
                        find_ncc_peak(refchip, imgchip_i1[cnt1], uv_pivot_cp, uvncc_cp);
                        dp_cp_sub[(cnt3-1)*8+(cnt2+1)*2]->val[cnt1][0]=uvncc_cp[0];
                        dp_cp_sub[(cnt3-1)*8+(cnt2+1)*2]->val[cnt1][1]=uvncc_cp[1];
                        dp_cp_sub[(cnt3-1)*8+(cnt2+1)*2]->val[cnt1][2]=uvncc_cp[2];

                        //swapped forward matching
                        //prepare for the refchip from i1
                        for(cntv=-ocw_refchip;cntv<=ocw_refchip;cntv++)
                        {
                            for(cntu=-ocw_refchip;cntu<=ocw_refchip;cntu++)
                            {
                                refchip->val[ocw_refchip+cntv][ocw_refchip+cntu]=imgchip_i1[cnt1]->val[ocw_chip+cntv][ocw_chip+cntu];
                            }
                        }
                        //prepare for the sarea from i0 -> imgchip_i0[cnt1]
                    
                        //perform feature tracking
                        //NOTE: No need to reverse the nc_pivot_cp (rectangular-shaped distribution)
                        find_ncc_peak(refchip, imgchip_i0[cnt1], uv_pivot_cp, uvncc_cp);
                        dp_cp_sub[(cnt3-1)*8+(cnt2+1)*2+1]->val[cnt1][0]=-uvncc_cp[0];
                        dp_cp_sub[(cnt3-1)*8+(cnt2+1)*2+1]->val[cnt1][1]=-uvncc_cp[1];
                        dp_cp_sub[(cnt3-1)*8+(cnt2+1)*2+1]->val[cnt1][2]=uvncc_cp[2];

                        
                    }//for(cnt1=0;cnt1<num_grid_sub;cnt1++)
                    GMA_float_destroy(refchip);
                }//#pragma omp parallel private(refchip,cntv,cntu,uvncc_cp) shared(dp_cp_sub)
            }//for(cnt3=1;cnt3<3;cnt3++)//counter for ocw
            //de-allocate the used image chips
            for(cnt1=0;cnt1<num_grid_sub;cnt1++)
            {
                GMA_float_destroy(imgchip_i0[cnt1]);
                GMA_float_destroy(imgchip_i1[cnt1]);
            }
            
        }//for(cnt2=-1;cnt2<=2;cnt2++)//counter for the kernel
        //printf("\n");
        //cluster the dp_cp_sub
        mvn_dp_cp_sub=calc_mean_var_num_dp_cluster(dp_cp_sub,16);
        free(imgchip_i0);
        free(imgchip_i1);

        //find the prominent displacements
        int32_t id_xyuvav;
        for(cnt1=0;cnt1<num_grid_sub;cnt1++)
        {
            int32_t cnt_cluster;
            for(cnt_cluster=0;cnt_cluster<mvn_dp_cp_sub[cnt1]->nrows;cnt_cluster++)
            {
                if(mvn_dp_cp_sub[cnt1]->val[cnt_cluster][4]>=0.6)
                {
                    id_xyuvav=(int32_t)(xyuvav_cp_sub->val[cnt1][6]);
                    sduv[0]+=mvn_dp_cp_sub[cnt1]->val[cnt_cluster][0];
                    sduv[1]+=mvn_dp_cp_sub[cnt1]->val[cnt_cluster][1];
                    flag_cp->val[id_xyuvav][0]=1;
                    num_cp_current++;
                }
            }
        }

        //de-allocate mvn_dp_cp_sub
        for(cnt1=0;cnt1<num_grid_sub;cnt1++)
        {
            GMA_float_destroy(mvn_dp_cp_sub[cnt1]);
        }
        free(mvn_dp_cp_sub);


        if(num_cp<=num_cp_current)
        {
            flag_sucessful_ft_cp=1;
            printf("Sufficient # of CP found: %d>=%d\n",num_cp_current,num_cp);
            GMA_double_destroy(xyuvav_cp_sub);
            break;
        }
        else
        {
            GMA_double_destroy(xyuvav_cp_sub);
            printf("#CP found so far: %d\n",num_cp_current);
        }
        
        //TODO: dp_cp_sub might be leak

        for(cnt1=0;cnt1<16;cnt1++)  //NOTE: the "16" means the # of matching attempts per grid in CP measurement
        {
            //dp_cp_sub[cnt1]=GMA_float_create(num_grid_sub,3);
            GMA_float_destroy(dp_cp_sub[cnt1]);
        }


    }//for(cnt=0;cnt<20;cnt++)
    free(array_index_segment);
    

    //Consider the image pair to be valid if if has more CPs than (or equal to) the least number of thresholds
    if(num_cp_current<num_cp && num_cp_current>=param_mimc2.num_cp_min)
    {
        flag_sucessful_ft_cp=1;
        printf("Least number of CPs found(%d>=%d). Considering this pair has enough CPs\n",num_cp_current,param_mimc2.num_cp_min);
    }

    //calculat the offset if the vmap pair is a valid one
    if(flag_sucessful_ft_cp)
    {
        float du_cp=sduv[0]/(float)num_cp_current;
        float dv_cp=sduv[1]/(float)num_cp_current;
        if(du_cp>0)
        {
            offset[0]=(int32_t)(du_cp+0.5);
        }
        else
        {
            offset[0]=(int32_t)(du_cp-0.5);
        }

        if(dv_cp>0)
        {
            offset[1]=(int32_t)(dv_cp+0.5);
        }
        else
        {
            offset[1]=(int32_t)(dv_cp-0.5);
        }
    }
    else
    {
        printf("Not enough # of successful CP measurement (%d<%d)\n",num_cp_current,param_mimc2.num_cp_min);
        GMA_double_destroy(xyuvav_cp);
        GMA_int32_destroy(uv_pivot_cp);
        return -1;
    }

    //deallocate the variables
    GMA_double_destroy(xyuvav_cp);
    GMA_int32_destroy(uv_pivot_cp);
    
    //free(flag_cp);

    //revert the original value of num_dp;
    num_dp=num_dp_original;
    return 1;
}

void GMA_double_randperm_row(GMA_double *var)
{   
    GMA_double *var_temp=GMA_double_create(var->nrows,var->ncols);

    int32_t cntrow,cntcol;
    int32_t limit_idx_row=var->nrows-1;
    int32_t idx_row;
    double N_A_N=sqrt(-1.0);

    var_temp=GMA_double_create(var->nrows,var->ncols);
    
    //duplicate the input GMA array. Also empty rhe
    for(cntrow=0;cntrow<var->nrows;cntrow++)
    {
        for(cntcol=0;cntcol<var->ncols;cntcol++)
        {
            var_temp->val[cntrow][cntcol]=var->val[cntrow][cntcol];
            var->val[cntrow][cntcol]=N_A_N;
        }
    }
    
    //deterimne the vector of the row order
    srand(time(NULL));
    
    for(limit_idx_row=var_temp->nrows-1;limit_idx_row>=0;limit_idx_row--)
    {
        //determine the ID to put
        if(limit_idx_row!=0)
        {
            idx_row=(int32_t)(rand()%limit_idx_row);
        }
        else
        {
            idx_row=0;
        }
        
        for(cntcol=0;cntcol<var->ncols;cntcol++)
        {
            var->val[limit_idx_row][cntcol]=var_temp->val[idx_row][cntcol];
            //move the head row value to the idx_row
            var_temp->val[idx_row][cntcol]=var_temp->val[0][cntcol];
            //move the tail row in var_temp to the head
            var_temp->val[0][cntcol]=var_temp->val[limit_idx_row][cntcol];
        }
    }
    GMA_double_destroy(var_temp);
}


GMA_int32** get_uv_pivot(GMA_double *xyuvav, float dt, param param_mimc2, int32_t ocw, GMA_float *i1)
{   //TODO make this module more efficient
    
    //unsigned int numgrid=xyuvav->nrows;
    int32_t cnt,cnt1;
    int32_t num_pivot;
    double length_pivot;
    float u,v,incr_u,incr_v,norm_incr;
    float theta;
    GMA_int32 **out=malloc(sizeof(GMA_int32*)*(xyuvav->nrows));
    //printf("Calculating uv_pivot\n");
    for(cnt=0;cnt<xyuvav->nrows;cnt++)
    {
        u=0.0;
        v=0.0;
        //TODO: Avoid using the triangular function. Replace them with the utilization of unit vector
        theta=atan2(xyuvav->val[cnt][5],xyuvav->val[cnt][4]);
        incr_u=cos(theta);
        incr_v=sin(theta);
        if(fabs(incr_u)>fabs(incr_v))   //normalize the increment (to avoid duplicated pivots)
        {
            incr_u=incr_u/fabs(incr_u);
            incr_v=incr_v/fabs(incr_u);
        }
        else
        {
            incr_u=incr_u/fabs(incr_v);
            incr_v=incr_v/fabs(incr_v);
        }
        norm_incr=sqrt(incr_u*incr_u+incr_v*incr_v);
        length_pivot=sqrt(xyuvav->val[cnt][4]*xyuvav->val[cnt][4]+xyuvav->val[cnt][5]*xyuvav->val[cnt][5])/param_mimc2.mpp/365*dt*param_mimc2.AW_SF+param_mimc2.AW_CRE+1;
        //count up the number of pivots
        num_pivot=0;
        while(u+(float)(xyuvav->val[cnt][2])-(float)ocw>0 &&
              u+(float)(xyuvav->val[cnt][2])+(float)ocw<(float)(i1->ncols-1) &&
              v+(float)(xyuvav->val[cnt][3])-(float)ocw>0 &&
              v+(float)(xyuvav->val[cnt][3])+(float)ocw<(float)(i1->nrows-1) &&
              length_pivot>(double)(norm_incr*(double)num_pivot))
        {
            num_pivot++;
            u+=incr_u;
            v+=incr_v;
        }
        
        u=0.0;
        v=0.0;
        out[cnt]=GMA_int32_create(num_pivot,2);
        out[cnt]->val[0][0]=0;
        out[cnt]->val[0][1]=0;
        for(cnt1=1;cnt1<num_pivot;cnt1++)
        {
            u+=incr_u;
            v+=incr_v;
            out[cnt]->val[cnt1][0]=(int32_t)(u+0.5);
            out[cnt]->val[cnt1][1]=-(int32_t)(v+0.5);
        }
    }
    printf("\n");
    return out;
}


char investigate_valid_grid(GMA_float *refchip,GMA_float *sarea)
{
    char out;
    int32_t num_invalid_refchip=0,num_invalid_sarea=0;
    float std_refchip,std_sarea;
    int32_t cntu,cntv;
    float numpx_refchip,numpx_sarea;
    float max_ratio_invalid=0.8;
    
    numpx_refchip=(float)(refchip->nrows*refchip->ncols);
    numpx_sarea=(float)(sarea->nrows*sarea->ncols);

    //check refchip
    for(cntu=0;cntu<refchip->ncols;cntu++)
    {
        for(cntv=0;cntv<refchip->nrows;cntv++)
        {
            if(refchip->val[cntv][cntu]<MIN_DN) num_invalid_refchip++;
        }
    }

    //check sarea
    for(cntu=0;cntu<sarea->ncols;cntu++)
    {
        for(cntv=0;cntv<sarea->nrows;cntv++)
        {
            if(sarea->val[cntv][cntu]<MIN_DN) num_invalid_sarea++;
        }
    }
    
    if((float)num_invalid_refchip/numpx_refchip>max_ratio_invalid || (float)num_invalid_sarea/numpx_sarea>max_ratio_invalid)
    {
        out=0;
    }
    else
    {
        out=1;
    }
    return out;
}


void find_ncc_peak(GMA_float *refchip, GMA_float *sarea, GMA_int32 *uv_pivot, float *uvncc)
{
    //TODO: A special mode of this function is necessary, especially when there is no void area line
    GMA_float *cmap;
    int32_t cnt_grid,cnt_pivot;
    int dx2,dy2,Dx2,Dy2,ocw;
    int duv[2],pivot[2];
    int uv_peak[2];
    float nccmax;
    int32_t cnt1,cnt2,cnt3,cnt4;
    int32_t uv0[2];
    const float N_A_N=sqrt(-1.0);
    int32_t nsample, flag_newncc;
    double sx,sy,sxx,syy,sxy;
    float ncc9[9];
    double coeff_poly[6];

    Dx2=(int)(sarea->ncols);
    Dy2=(int)(sarea->nrows);
    dx2=Dx2/2;
    dy2=Dy2/2;
    ocw=(int)(refchip->nrows)/2;
    uv_peak[0]=dx2;
    uv_peak[1]=dy2;
        
    uvncc[0]=0.0;
    uvncc[1]=0.0;
    uvncc[2]=-2.0;
    
    //create cmap and sarea
    cmap=GMA_float_create(Dy2,Dx2);
    for(cnt1=-dx2;cnt1<dx2;cnt1++)  //TODO: re-think about the loop boundary
    {
        for(cnt2=-dy2;cnt2<dy2;cnt2++) cmap->val[cnt2+dy2][cnt1+dx2]=-2.0;
    }

    if(!investigate_valid_grid(refchip,sarea))
    {
        uvncc[0]=N_A_N;
        uvncc[1]=N_A_N;
        uvncc[2]=-3;
    }
    else
    {
        for(cnt_pivot=0;cnt_pivot<uv_pivot->nrows;cnt_pivot++)    //loop through the pivot point
        {
            pivot[0]=uv_pivot->val[cnt_pivot][0]+dx2; //initial pivot
            pivot[1]=uv_pivot->val[cnt_pivot][1]+dy2;
            nccmax=-2;
            duv[0]=-1;
            duv[1]=-1;
            flag_newncc=1;  //fake values to get into the while loop
            while((duv[0]!=0||duv[1]!=0)&&flag_newncc!=0)
            {   
                duv[0]=0;
                duv[1]=0;
                if(pivot[0]-ocw<=1 || pivot[0]+ocw>=Dx2-1 || pivot[1]-ocw<=1 || pivot[1]+ocw>=Dy2-1)  //boundary check
                {

                    break;
                }
                flag_newncc=0;
                for(cnt1=-1;cnt1<=1;cnt1++)
                {
                    for(cnt2=-1;cnt2<=1;cnt2++)
                    {
                        if(cmap->val[pivot[1]+cnt2][pivot[0]+cnt1]<-1.0) //when NCC is not calculated
                        {
                            //calculate NCC for the location [pivot[0]+cnt1][pivot[1]+cnt2]
                            flag_newncc++;
                            nsample=0;
                            sy=0; sx=0; sxx=0; sxy=0; syy=0;
                            for(cnt3=-ocw;cnt3<=ocw;cnt3++)
                            {
                                for(cnt4=-ocw;cnt4<=ocw;cnt4++)
                                {
                                    if(refchip->val[cnt4+ocw][cnt3+ocw]>=MIN_DN && sarea->val[pivot[1]+cnt2+cnt4][pivot[0]+cnt1+cnt3] >=MIN_DN)
								    {	//Concept of null exclusion. Refer to Ahn and Howat [2011]
                                        nsample++;
                                        sy+=sarea->val[pivot[1]+cnt2+cnt4][pivot[0]+cnt1+cnt3];
                                        sx+=refchip->val[cnt4+ocw][cnt3+ocw];
                                        sxx+=refchip->val[cnt4+ocw][cnt3+ocw]*refchip->val[cnt4+ocw][cnt3+ocw];
                                        syy+=sarea->val[pivot[1]+cnt2+cnt4][pivot[0]+cnt1+cnt3]*sarea->val[pivot[1]+cnt2+cnt4][pivot[0]+cnt1+cnt3];
                                        sxy+=refchip->val[cnt4+ocw][cnt3+ocw]*sarea->val[pivot[1]+cnt2+cnt4][pivot[0]+cnt1+cnt3];
                                    }
                                }
                            }
						    cmap->val[pivot[1]+cnt2][pivot[0]+cnt1]=(float)((nsample*sxy-sx*sy)/sqrt((nsample*sxx-sx*sx)*(nsample*syy-sy*sy))); //ncc calculation in spatial domain
                        }
                        if(cmap->val[pivot[1]+cnt2][pivot[0]+cnt1] > nccmax)
                        { 
                            nccmax=cmap->val[pivot[1]+cnt2][pivot[0]+cnt1];
                            duv[0]=cnt1;
                            duv[1]=cnt2;
                        }
                    }
                }
                pivot[0]+=duv[0];
                pivot[1]+=duv[1];
            }
            if(nccmax>uvncc[2])
            {
                uv_peak[0]=pivot[0];
                uv_peak[1]=pivot[1];
                uvncc[2]=nccmax;
            }
        }   //peak location along with the NCC is now at the array uvncc

        // quadratic fitting of the NCC peak
        
        ncc9[0]=cmap->val[uv_peak[1]-1][uv_peak[0]-1];
        ncc9[1]=cmap->val[uv_peak[1]-1][uv_peak[0]];
        ncc9[2]=cmap->val[uv_peak[1]-1][uv_peak[0]+1];
        ncc9[3]=cmap->val[uv_peak[1]][uv_peak[0]-1];
        ncc9[4]=cmap->val[uv_peak[1]][uv_peak[0]];
        ncc9[5]=cmap->val[uv_peak[1]][uv_peak[0]+1];
        ncc9[6]=cmap->val[uv_peak[1]+1][uv_peak[0]-1];
        ncc9[7]=cmap->val[uv_peak[1]+1][uv_peak[0]];
        ncc9[8]=cmap->val[uv_peak[1]+1][uv_peak[0]+1];
        //TODO: deal with the case that the peak is at the edges of the cmap
        
        //printf("Calculating Coefficients...");
        coeff_poly[0]=6*ncc9[0] -12*ncc9[1] +6*ncc9[2] +6*ncc9[3] -12*ncc9[4] +6*ncc9[5] +6*ncc9[6] -12*ncc9[7] +6*ncc9[8];
        coeff_poly[1]=9*ncc9[0] -9*ncc9[2] -9*ncc9[6] +9*ncc9[8];
        coeff_poly[2]=6*ncc9[0] +6*ncc9[1] +6*ncc9[2] -12*ncc9[3] -12*ncc9[4] -12*ncc9[5] +6*ncc9[6] +6*ncc9[7] +6*ncc9[8];
        coeff_poly[3]=-6*ncc9[0] +6*ncc9[2] -6*ncc9[3] +6*ncc9[5] -6*ncc9[6] +6*ncc9[8];
        coeff_poly[4]=-6*ncc9[0] -6*ncc9[1] -6*ncc9[2] +6*ncc9[6] +6*ncc9[7] +6*ncc9[8];
        coeff_poly[5]=-4*ncc9[0] +8*ncc9[1] -4*ncc9[2] +8*ncc9[3] +20*ncc9[4] +8*ncc9[5] -4*ncc9[6] +8*ncc9[7] -4*ncc9[8];
        coeff_poly[0]/=36;
        coeff_poly[1]/=36;
        coeff_poly[2]/=36;
        coeff_poly[3]/=36;
        coeff_poly[4]/=36;
        coeff_poly[5]/=36;

        //printf("Calculating fitting results...");
        uvncc[0]=-2*coeff_poly[2]*coeff_poly[3]+coeff_poly[1]*coeff_poly[4];
        uvncc[1]=-2*coeff_poly[0]*coeff_poly[4]+coeff_poly[1]*coeff_poly[3];
        uvncc[0]/=4*coeff_poly[0]*coeff_poly[2]-coeff_poly[1]*coeff_poly[1];
        uvncc[1]/=4*coeff_poly[0]*coeff_poly[2]-coeff_poly[1]*coeff_poly[1];
        uvncc[0]+=(float)(uv_peak[0]-dx2);
        uvncc[1]+=(float)(uv_peak[1]-dy2);
        
        //TODO: also calculate fitted NCC peak value if there is no significant computation resource usage
        //TODO: Cancel the sub-pixel interpolation if the fitting result goes too far from the original peak result
        //printf("Successful! - [%f, %f]\n\n",uvncc[0],uvncc[1]);
        
        
    }
    GMA_float_destroy(cmap);
    //printf("\n");
    


}


//GMA_float* matching_ncc_dlc_2(GMA_float *i0,GMA_float *i1,GMA_double *xyuvav,GMA_int32 **uv_pivot,int32_t ocw,float AW_CRE,float AW_SF)
GMA_float* matching_ncc_dlc_2(GMA_float *i0,GMA_float *i1,GMA_double *xyuvav,int32_t *offset,GMA_int32 **uv_pivot,int32_t ocw,float AW_CRE,float AW_SF)
{
    GMA_float *out=GMA_float_create(xyuvav->nrows,3); //n*3 matrix
    GMA_float *refchip,*sarea;
    float uvncc[3];
    int32_t cnt_grid;
    //unsigned int uv0[2];
    int32_t uv0[2];
    int32_t cnt1,cnt2,cnt3,cnt4;

    //refchip=GMA_float_create((unsigned int)ocw*2+1,(unsigned int)ocw*2+1);
    #pragma omp parallel private(uv0,uvncc,refchip,sarea) shared(out) 
    {
        refchip=GMA_float_create(ocw*2+1,ocw*2+1);
        #pragma omp for schedule(dynamic)
        for(cnt_grid=0;cnt_grid<xyuvav->nrows;cnt_grid++)
        {
            uv0[0]=(int32_t)xyuvav->val[cnt_grid][2];
            uv0[1]=(int32_t)xyuvav->val[cnt_grid][3];

            extract_refchip(i0,uv0,ocw,refchip);

            uv0[0]+=offset[0];
            uv0[1]+=offset[1];
            
            sarea=extract_sarea(i1, uv0, ocw, uv_pivot[cnt_grid]);
            find_ncc_peak(refchip, sarea, uv_pivot[cnt_grid], uvncc);
            out->val[cnt_grid][0]=uvncc[0];
            out->val[cnt_grid][1]=uvncc[1];
            out->val[cnt_grid][2]=uvncc[2];
            
            GMA_float_destroy(sarea);

        }
        GMA_float_destroy(refchip);
    }
    return out;
}


void extract_refchip(GMA_float *i0, int32_t *uv0, int32_t ocw,GMA_float *refchip)
{
    int32_t cnt1,cnt2;
    for(cnt1=-ocw;cnt1<=ocw;cnt1++)
    {
        for(cnt2=-ocw;cnt2<=ocw;cnt2++)
        {
            refchip->val[cnt2+ocw][cnt1+ocw]=i0->val[uv0[1]+cnt2][uv0[0]+cnt1];
        }
    }
}

GMA_float* extract_sarea(GMA_float *i1, int32_t *uv0, int32_t ocw, GMA_int32 *uv_pivot)
{
    GMA_float *sarea;
    int32_t cnt1,cnt2;
    int dx2, Dx2, dy2, Dy2;

    dx2=abs(uv_pivot->val[uv_pivot->nrows-1][0])+ocw+2;
    dy2=abs(uv_pivot->val[uv_pivot->nrows-1][1])+ocw+2;
    Dx2=dx2*2+1;
    Dy2=dy2*2+1;
    sarea=GMA_float_create(Dy2,Dx2);

    for(cnt1=-dx2;cnt1<dx2;cnt1++)
    {
        for(cnt2=-dy2;cnt2<dy2;cnt2++)
        {
            int32_t coord_u,coord_v;
            coord_u=uv0[0]+cnt1;
            coord_v=uv0[1]+cnt2;
            //boundary check
            if(coord_u>=0 && coord_u<i1->ncols && coord_v>=0 && coord_v<i1->nrows)
            {
                sarea->val[cnt2+dy2][cnt1+dx2]=i1->val[coord_v][coord_u];
            }
            else
            {
                sarea->val[cnt2+dy2][cnt1+dx2]=0.0;
            }
            
        } 
    }

    return sarea;
}


GMA_float** mimc2_postprocess(GMA_float **dp, GMA_double *xyuvav, float dt)
{
    int32_t cnt;
    GMA_float **vxyexyqual;
    GMA_float **mvn_dp;
    GMA_int32 *dpf0;    //TODO: use GMA_int8 instead of GMA_int32 (for memory efficiency)
    GMA_int32 *ruv_neighbor;
    GMA_uint8 *mask_nomatching;
    GMA_float *dpf_dx, *dpf_dy;
    float N_A_N=sqrt(-1.0);

    printf("Calculating mean/std/qual (MSQ):\n");
    mvn_dp=calc_mean_var_num_dp_cluster(dp,num_dp);
    printf("MSQ calculation successful\n\n");

    printf("Finding grids without any matching results");
    mask_nomatching=get_mask_nomatching(mvn_dp);

    //determine dpf0
    printf("Determining prominent displacements\n");
    dpf0=get_dpf0(mvn_dp,0.6);      //TODO: soft-code the threshold
    printf("Successful\n\n");

    //GMA_int32_save("../MIMC2_C_GMA/dpf0.GMA",dpf0);    //for checking the result 
    //Debug code
    
    printf("Calculating meighboring grid mask for dpf2\n");
    ruv_neighbor=get_ruv_neighbor(xyuvav,param_mimc2.radius_neighbor_dpf1);
    printf("Successful: %d neignboring grids\n\n",ruv_neighbor->nrows);

    printf("Calculating initial most probable displacements (dpf2)\n");
    dpf_dx=GMA_float_create(dimy_vmap,dimx_vmap);
    dpf_dy=GMA_float_create(dimy_vmap,dimx_vmap);
    
    get_dpf1(dpf0,dpf_dx,dpf_dy,ruv_neighbor,mvn_dp,xyuvav);
    printf("Successful\n\n");

    //Perform pseudosmoothing
    printf("Performing pseudosmoothing\n");
    GMA_int32_destroy(ruv_neighbor);
    ruv_neighbor=get_ruv_neighbor(xyuvav,param_mimc2.radius_neighbor_ps);
    get_dpf_pseudosmoothing(dpf0,dpf_dx,dpf_dy,ruv_neighbor,mvn_dp,xyuvav);
    printf("pseudosmoothing successful!\n\n");

    //convert the cluster map (i.e. dpf0) into vmap (i.e. vxyexyqual)
    vxyexyqual=(GMA_float**)malloc(sizeof(GMA_float*)*5);
    vxyexyqual[0]=GMA_float_create(dimy_vmap,dimx_vmap);    //vx
    vxyexyqual[1]=GMA_float_create(dimy_vmap,dimx_vmap);    //vy
    vxyexyqual[2]=GMA_float_create(dimy_vmap,dimx_vmap);    //ex
    vxyexyqual[3]=GMA_float_create(dimy_vmap,dimx_vmap);    //ey
    vxyexyqual[4]=GMA_float_create(dimy_vmap,dimx_vmap);    //qual
    
    printf("checkpoint 1\n");
    
    int32_t cntv,cntu;
    int32_t id_vec;
    int32_t id_cluster;
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            id_vec=cntv*dimx_vmap+cntu;
            id_cluster=dpf0->val[cntv][cntu];
            if(id_cluster>=0)
            {
                vxyexyqual[0]->val[cntv][cntu]=mvn_dp[id_vec]->val[id_cluster][0];
                vxyexyqual[1]->val[cntv][cntu]=mvn_dp[id_vec]->val[id_cluster][1];
                vxyexyqual[2]->val[cntv][cntu]=mvn_dp[id_vec]->val[id_cluster][2];
                vxyexyqual[3]->val[cntv][cntu]=mvn_dp[id_vec]->val[id_cluster][3];
                vxyexyqual[4]->val[cntv][cntu]=mvn_dp[id_vec]->val[id_cluster][4];
            }
            else
            {
                vxyexyqual[0]->val[cntv][cntu]=N_A_N;
                vxyexyqual[1]->val[cntv][cntu]=N_A_N;
                vxyexyqual[2]->val[cntv][cntu]=N_A_N;
                vxyexyqual[3]->val[cntv][cntu]=N_A_N;
                vxyexyqual[4]->val[cntv][cntu]=N_A_N;
            }
        }
    }
    printf("checkpoint 2\n");
    //deallocate the allocated variables

    GMA_float_destroy(dpf_dx);
    GMA_float_destroy(dpf_dy);
    GMA_int32_destroy(ruv_neighbor);
    GMA_int32_destroy(dpf0);
    GMA_uint8_destroy(mask_nomatching);
    printf("checkpoint 3\n");
    for(cnt=0;cnt<num_grid;cnt++)
    {
        GMA_float_destroy(mvn_dp[cnt]);
    }
    free(mvn_dp);
    printf("checkpoint 4\n");
    return vxyexyqual;

}


GMA_float** calc_mean_var_num_dp_cluster(GMA_float** dp,int32_t num_dpoi)
{
    int32_t cnt_grid,cnt_dp,num_valid_dp;
    uint8_t cnt_id,max_id;
    int32_t num_grid=dp[0]->nrows;
    GMA_float **out;
    GMA_uint8 id_cluster;
    GMA_float dp_stack;
    float min_ncc=0.1;
    float sx[num_dpoi],sxx[num_dpoi],sy[num_dpoi],syy[num_dpoi];
    int nsample[num_dpoi];
    //float *sx, *sxx, *sy, *syy;
    //int *nsample;
    //sx=malloc(sizeof(float)*num_dpoi);
    //sxx=malloc(sizeof(float)*num_dpoi);
    //sy=malloc(sizeof(float)*num_dpoi);
    //syy=malloc(sizeof(float)*num_dpoi);
    //nsample=malloc(sizeof(int)*num_dpoi);

    
    float min_dist=0.5; //TODO: soft-code this

    //printf("Allocating output values\n");
    out=(GMA_float**)malloc(sizeof(GMA_float*)*dp[0]->nrows);

    //printf("Initializing dp_stack\n");
    //initialze dp_stack
    dp_stack.nrows=num_dpoi;
    dp_stack.ncols=2;
    dp_stack.val=(float**)malloc(sizeof(float*)*dp_stack.nrows);
    for(cnt_dp=0;cnt_dp<num_dpoi;cnt_dp++)
    {
        dp_stack.val[cnt_dp]=(float*)malloc(sizeof(float)*dp_stack.ncols);
    }
    //TODO: modify the part above to be in compliance with the struct definition (i.e. allocate space at dp_stack->data)

    //initialize id_cluster
    //printf("Initializing id_cluster\n");
    id_cluster.nrows=num_dpoi;
    id_cluster.ncols=1;
    id_cluster.val=(uint8_t**)malloc(sizeof(uint8_t*)*dp_stack.nrows);
    for(cnt_dp=0;cnt_dp<num_dpoi;cnt_dp++)
    {
        id_cluster.val[cnt_dp]=(uint8_t*)malloc(sizeof(float)*id_cluster.ncols);
    }

    //TODO: parallelize this part
    //printf("Looping through grid\n");
    for(cnt_grid=0;cnt_grid<num_grid;cnt_grid++)
    {
        //build a stack
        num_valid_dp=0;
        for(cnt_dp=0;cnt_dp<num_dpoi;cnt_dp++)
        {
            if(dp[cnt_dp]->val[cnt_grid][2]>min_ncc)
            {
                dp_stack.val[num_valid_dp][0]=dp[cnt_dp]->val[cnt_grid][0];
                dp_stack.val[num_valid_dp][1]=dp[cnt_dp]->val[cnt_grid][1];
                num_valid_dp++;
            }
            //trim out the tail
            dp_stack.nrows=num_valid_dp;
            id_cluster.nrows=num_valid_dp;
        }

        //perform clustering
        cluster_euclidian(&dp_stack, min_dist, &id_cluster);

        //find the maximum value in id_cluster
        max_id=0;
        for(cnt_dp=0;cnt_dp<id_cluster.nrows;cnt_dp++)
        {
            if(max_id<id_cluster.val[cnt_dp][0])
            {
                max_id=id_cluster.val[cnt_dp][0];
            }
        }
        
        out[cnt_grid]=GMA_float_create(max_id,5);
        for(cnt_dp=0;cnt_dp<max_id;cnt_dp++)
        {
            sx[cnt_dp]=0;
            sy[cnt_dp]=0;
            sxx[cnt_dp]=0;
            syy[cnt_dp]=0;
            nsample[cnt_dp]=0;
        }
        
        //[1111111, 2222222, 333333, 444444, 5555555555]
        //[mean_vx, mean_vy, var_vx, var_vy, num_sample]
        for(cnt_dp=0;cnt_dp<id_cluster.nrows;cnt_dp++)
        {
            sx[id_cluster.val[cnt_dp][0]-1]+=dp_stack.val[cnt_dp][0];
            sy[id_cluster.val[cnt_dp][0]-1]+=dp_stack.val[cnt_dp][1];
            sxx[id_cluster.val[cnt_dp][0]-1]+=dp_stack.val[cnt_dp][0]*dp_stack.val[cnt_dp][0];
            syy[id_cluster.val[cnt_dp][0]-1]+=dp_stack.val[cnt_dp][1]*dp_stack.val[cnt_dp][1];
            nsample[id_cluster.val[cnt_dp][0]-1]+=1;
        }

        //each row contains the statistics of each cluster
        for(cnt_dp=0;cnt_dp<max_id;cnt_dp++)
        {
            out[cnt_grid]->val[cnt_dp][0]=sx[cnt_dp]/(float)nsample[cnt_dp];  //mean vx
            out[cnt_grid]->val[cnt_dp][1]=sy[cnt_dp]/(float)nsample[cnt_dp];  //mean vy
            out[cnt_grid]->val[cnt_dp][2]=sxx[cnt_dp]/(float)nsample[cnt_dp]
                                          -out[cnt_grid]->val[cnt_dp][0]*out[cnt_grid]->val[cnt_dp][0];  //var vx
            out[cnt_grid]->val[cnt_dp][3]=syy[cnt_dp]/(float)nsample[cnt_dp]
                                          -out[cnt_grid]->val[cnt_dp][1]*out[cnt_grid]->val[cnt_dp][1];  //var vy
            out[cnt_grid]->val[cnt_dp][4]=(float)nsample[cnt_dp]/(float)num_dpoi;  //# of samples
        }
        
        //TODO: Consider the way about reducing the overhead from memory allocation/deallocation
    }

    //deallocate dp_stack
    dp_stack.nrows=num_dpoi;
    id_cluster.nrows=num_dpoi;
    for(cnt_dp=0;cnt_dp<num_dpoi;cnt_dp++)
    {
        free(dp_stack.val[cnt_dp]);
    }
    free(dp_stack.val);
    
    //deallocate id_cluster
    for(cnt_dp=0;cnt_dp<num_dpoi;cnt_dp++)
    {
        free(id_cluster.val[cnt_dp]);
    }
    free(id_cluster.val);
    return out;

    //free(sx);
    //free(sxx);
    //free(sy);
    //free(syy);
    //free(nsample);
}


void cluster_euclidian(GMA_float *dp_stack, float min_dist, GMA_uint8 *id_cluster)
{   //TODO : improve the efficiency
    //GMA_uint8 *out;
    GMA_uint8 *mtrx_dist;
    float dx,dy;
    float min_dist_sq=min_dist*min_dist;
    uint8_t id_curr;
    int32_t cnt1,cnt2;

    //initialize the id_cluster
    for(cnt1=0;cnt1<id_cluster->nrows;cnt1++)
    {
        id_cluster->val[cnt1][0]=0;
    }

    //calculate the distance matrix
    mtrx_dist=GMA_uint8_create(dp_stack->nrows,dp_stack->nrows);
    for(cnt1=0;cnt1<dp_stack->nrows;cnt1++)
    {
        for(cnt2=cnt1;cnt2<dp_stack->nrows;cnt2++)
        {
            dx=dp_stack->val[cnt2][0]-dp_stack->val[cnt1][0];
            dy=dp_stack->val[cnt2][1]-dp_stack->val[cnt1][1];
            if(dx*dx+dy*dy<min_dist_sq)
            {
                mtrx_dist->val[cnt1][cnt2]=1;
                mtrx_dist->val[cnt2][cnt1]=1;
            }
            else
            {
                mtrx_dist->val[cnt1][cnt2]=0;
                mtrx_dist->val[cnt2][cnt1]=0;
            }
        }
    }

    id_curr=0;
    for(cnt1=0;cnt1<dp_stack->nrows;cnt1++)
    {
        if(id_cluster->val[cnt1][0]==0)
        {
            id_curr++;
            mark_row(cnt1, id_curr, id_cluster, mtrx_dist);
        }
    }

    GMA_uint8_destroy(mtrx_dist);
}


void mark_row(int32_t nrow, uint8_t id_dp, GMA_uint8 *id_cluster, GMA_uint8 *mtrx_dist)
{
    int32_t cnt;
    for(cnt=0;cnt<id_cluster->nrows;cnt++)
    {
        if(mtrx_dist->val[nrow][cnt] && id_cluster->val[cnt][0]==0)
        {
            id_cluster->val[cnt][0]=id_dp;
            mark_row(cnt, id_dp, id_cluster, mtrx_dist);
        }
    }
}


GMA_uint8* get_mask_nomatching(GMA_float **mvn_dp)
{
    int32_t cntu,cntv,cntvec;
    GMA_uint8 *out;
    out=GMA_uint8_create(dimy_vmap,dimx_vmap);

    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            cntvec=cntv*dimx_vmap+cntu;
            if(mvn_dp[cntvec]->nrows==0)
            {
                out->val[cntv][cntu]=1;
                
            }
            else
            {
                out->val[cntv][cntu]=0;
            }
        }
    }
    //GMA_uint8_save("../MIMC2_C_GMA/mask_nomatching.GMA",out);
    return out;
}


GMA_int32* get_dpf0(GMA_float **mvn_dp, float min_matching_ratio)   //to get the most probable dps
{   //input parameter description  [000][111][222][333][4444]
    //mvn_dp: array of GMA_float*. [mvx][mvy][evx][evy][qual]
    //min_matching_ratio: Minimum value of matching ratio to be considered as a prominent displacement (dpf0).
    //                    This value is usually 0.6
    //output value description
    //out: [dimy_vmap][dimx_vmap] GMA_int32. Each cell contains the cluster ID of the DPF0
    //     The cell value is -1 (or less than zero) if there isn`t any corresponding cluster
    
    int32_t cnt_grid,cnt_cluster;
    int32_t cntx,cnty;
    uint8_t flag_assigned;
    GMA_int32 *out=GMA_int32_create(dimy_vmap,dimx_vmap);

    //initialize the output
    for(cnty=0;cnty<dimy_vmap;cnty++)
    {
        for(cntx=0;cntx<dimx_vmap;cntx++)
        {
            cnt_grid=cnty*dimx_vmap+cntx;
            flag_assigned=0;
            for(cnt_cluster=0;cnt_cluster<mvn_dp[cnt_grid]->nrows;cnt_cluster++)
            {
                if(mvn_dp[cnt_grid]->val[cnt_cluster][4]>min_matching_ratio)
                {
                    out->val[cnty][cntx]=cnt_cluster;
                    flag_assigned=1;
                    break;
                }
            }
            if(!flag_assigned)
            {
                out->val[cnty][cntx]=-1;
            }
        }
    }
    //GMA_int32_save("../MIMC2_C_GMA/dpf0.GMA",out);
    return out;
    
}


GMA_int32* get_ruv_neighbor(GMA_double *xyuvav,float radius_neighbor)
{
    GMA_int32 *out, *vec_buffer;
    //unsigned int dimu,dimv,cu,cv;   //cu, cv: center of the 2d grid
    int32_t cu,cv;   //cu, cv: center of the 2d grid
    int32_t cntu,cntv;
    int32_t num_grid_neighbor;
    float dxy[2];
    GMA_float *field_x, *field_y;
    float sq_dist,cx,cy;        //xy coordinate of the 2d grid

    cu=dimx_vmap/2;
    cv=dimy_vmap/2;

    // create the mesh grid of x and y
    field_x=GMA_float_create(dimy_vmap,dimx_vmap);
    field_y=GMA_float_create(dimy_vmap,dimx_vmap);
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            field_x->val[cntv][cntu]=xyuvav->val[cntu][0];
            field_y->val[cntv][cntu]=xyuvav->val[cntv*dimx_vmap][1];
        }
    }
    
    cx=field_x->val[cv][cu];
    cy=field_y->val[cv][cu];

    vec_buffer=GMA_int32_create(num_grid,2);

    num_grid_neighbor=0;
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            dxy[0]=field_x->val[cntv][cntu]-cx;
            dxy[1]=field_y->val[cntv][cntu]-cy;
            sq_dist=dxy[0]*dxy[0]+dxy[1]*dxy[1];

            if(sq_dist<=(radius_neighbor*param_mimc2.meter_per_spacing)*(radius_neighbor*param_mimc2.meter_per_spacing))
            {
                vec_buffer->val[num_grid_neighbor][0]=cntu-cu;
                vec_buffer->val[num_grid_neighbor][1]=cntv-cv;
                num_grid_neighbor++;
            }
        }
    }

    out=GMA_int32_create(num_grid_neighbor,2);
    for(cntv=0;cntv<num_grid_neighbor;cntv++)
    {
        out->val[cntv][0]=vec_buffer->val[cntv][0];
        out->val[cntv][1]=vec_buffer->val[cntv][1];
    }

    GMA_int32_destroy(vec_buffer);
    GMA_float_destroy(field_x);
    GMA_float_destroy(field_y);

    return out;
}


void get_dpf1(GMA_int32 *dpf0, GMA_float *dpf_dx, GMA_float *dpf_dy, GMA_int32 *ruv_neighbor, GMA_float **mvn_dp, GMA_double *xyuvav)
{
    int32_t cntu,cntv,cnt_grid,cnt_neighbor,cnt_cluster,id_cluster;
    int32_t NOI,num_processed,num_grid_neighbor,num_grid_unprocessed;
    int32_t uvoi[2];
    GMA_float *dx_buffer, *dy_buffer;
    GMA_float *vec_ruv_w_scale_mag;
    GMA_float *mtrx_noi;
    float dpe[2],dxy_neighbor[2],dxy[2],apv_neighbor[2];
    float mag_dpe,mag_dxy_neighbor;
    float N_A_N=sqrt(-1.0);
    float sum_mag,sq_dist,sq_dist_closest;
    char buffer_filename[2048];

    
    dx_buffer=GMA_float_create(dimy_vmap,dimx_vmap);
    dy_buffer=GMA_float_create(dimy_vmap,dimx_vmap);
    
    mtrx_noi=GMA_float_create(dimy_vmap,dimx_vmap);
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            mtrx_noi->val[cntv][cntu]=1.0;
        }
    }

    

    //vec_ruv_w_scale_mag=GMA_float_create(ruv_neighbor->nrows,4);
    
    //reconstruct dpf_dx and dpf_dy. Also initialize the dx_buffer and dy_buffer
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            cnt_grid=cntv*dimx_vmap+cntu;
            if(dpf0->val[cntv][cntu]>=0)
            {
                dpf_dx->val[cntv][cntu]=mvn_dp[cnt_grid]->val[dpf0->val[cntv][cntu]][0];
                dpf_dy->val[cntv][cntu]=mvn_dp[cnt_grid]->val[dpf0->val[cntv][cntu]][1];
            }
            else
            {
                dpf_dx->val[cntv][cntu]=N_A_N;
                dpf_dy->val[cntv][cntu]=N_A_N;
            }
            dx_buffer->val[cntv][cntu]=N_A_N;
            dy_buffer->val[cntv][cntu]=N_A_N;
        }
    }
    
    //debugging code
    //GMA_float_save("/Users/seongsu/Desktop/bad_vmaptar_example/dpf0_gx.gma",dpf_dx);
    //GMA_float_save("/Users/seongsu/Desktop/bad_vmaptar_example/dpf0_gy.gma",dpf_dy);


    NOI=0;
    num_grid_unprocessed=1;
    int32_t thres_num_grid_neighbor;
    for(thres_num_grid_neighbor=ruv_neighbor->nrows-1;thres_num_grid_neighbor>=3;thres_num_grid_neighbor--)
    {
        float thres_weight=0.5;
        float factor_mpy_to_px=1.0/365.0*dt/param_mimc2.mpp;
        while(num_grid_unprocessed!=0 && thres_weight>=0.5)
        {
            thres_weight-=0.02;
            num_processed=1;
            while(num_processed!=0)
            {   
                NOI++;
                num_processed=0;
                //#pragma omp parallel private(cnt_grid,id_cluster,num_grid_neighbor,uvoi, vec_ruv_w_scale_mag, dpe, dxy_neighbor, dxy, sq_dist, sq_dist_closest, sum_mag) \
                                       shared(dpf0, dpf_dx, dpf_dy, ruv_neighbor, mvn_dp, xyuvav, dx_buffer, dy_buffer, N_A_N)
                // {
                vec_ruv_w_scale_mag=GMA_float_create(ruv_neighbor->nrows,7);    
                //#pragma omp for schedule(dynamic) collapse(2) 
                for(cntv=0;cntv<dimy_vmap;cntv++)
                {   //loop around the array and process if the visited grid needs to
                    for(cntu=0;cntu<dimx_vmap;cntu++)
                    {
                        cnt_grid=cntv*dimx_vmap+cntu;
                        if(isnan(dpf_dx->val[cntv][cntu]+dpf_dy->val[cntv][cntu]) && mvn_dp[cnt_grid]->nrows!=0)
                        {      //dx validity check             dy validity check                    nomatching check
                            sum_mag=0;
                            num_grid_neighbor=0;
                            dpe[0]=xyuvav->val[cnt_grid][4]*factor_mpy_to_px;
                            dpe[1]=-xyuvav->val[cnt_grid][5]*factor_mpy_to_px;
                            mag_dpe=sqrt(dpe[0]*dpe[0]+dpe[1]*dpe[1]);
                            //printf("factor=%f,dpe=[%e,%e]\n",factor_mpy_to_px,dpe[0],dpe[1]);
                            for(cnt_neighbor=0;cnt_neighbor<ruv_neighbor->nrows;cnt_neighbor++)
                            {
                                uvoi[0]=cntu+ruv_neighbor->val[cnt_neighbor][0];
                                uvoi[1]=cntv+ruv_neighbor->val[cnt_neighbor][1];
                                if(uvoi[0]>=0 && uvoi[0]<dimx_vmap && uvoi[1]>=0 && uvoi[1]<dimy_vmap) //boundary check
                                {
                                    dxy_neighbor[0]=dpf_dx->val[uvoi[1]][uvoi[0]];
                                    dxy_neighbor[1]=dpf_dy->val[uvoi[1]][uvoi[0]];
                                    mag_dxy_neighbor=sqrt(dxy_neighbor[0]*dxy_neighbor[0]+dxy_neighbor[1]*dxy_neighbor[1]);
                                    if(!isnan(dxy_neighbor[0]+dxy_neighbor[1]))
                                    //TODO: change: if(!isnan(dxy_neighbor[0]+dxy_neighbor[1]))
                                    {
                                        apv_neighbor[0]=(float)(xyuvav->val[uvoi[1]*dimx_vmap+uvoi[0]][4])*factor_mpy_to_px;
                                        apv_neighbor[1]=-(float)(xyuvav->val[uvoi[1]*dimx_vmap+uvoi[0]][5])*factor_mpy_to_px;

                                        vec_ruv_w_scale_mag->val[num_grid_neighbor][0]=(float)ruv_neighbor->val[cnt_neighbor][0];//relative u coordinate
                                        vec_ruv_w_scale_mag->val[num_grid_neighbor][1]=(float)ruv_neighbor->val[cnt_neighbor][1];//relative v coordinate
                                        vec_ruv_w_scale_mag->val[num_grid_neighbor][4]=sqrt(dxy_neighbor[0]*dxy_neighbor[0]+dxy_neighbor[1]*dxy_neighbor[1]);//magnitude
                                        vec_ruv_w_scale_mag->val[num_grid_neighbor][5]=sqrt(apv_neighbor[0]*apv_neighbor[0]+apv_neighbor[1]*apv_neighbor[1]);//magnitude
                                        vec_ruv_w_scale_mag->val[num_grid_neighbor][6]=mtrx_noi->val[uvoi[1]][uvoi[0]];//magnitude
                                        //vec_ruv_w_scale_mag->val[num_grid_neighbor][3]=sqrt((dxy_neighbor[0]*dxy_neighbor[0]+dxy_neighbor[1]*dxy_neighbor[1])
                                        //                                             /(apv_neighbor[0]*apv_neighbor[0]+apv_neighbor[1]*apv_neighbor[1]));//scale factor
                                        vec_ruv_w_scale_mag->val[num_grid_neighbor][3]=vec_ruv_w_scale_mag->val[num_grid_neighbor][4]
                                                                                        /sqrt(apv_neighbor[0]*apv_neighbor[0]+apv_neighbor[1]*apv_neighbor[1]);//scale factor
                                        
                                        //vec_ruv_w_scale_mag->val[num_grid_neighbor][3]=sqrt(dxy_neighbor[0]*dxy_neighbor[0]+dxy_neighbor[1]*dxy_neighbor[1])-
                                        //                                         sqrt(apv_neighbor[0]*apv_neighbor[0]+apv_neighbor[1]*apv_neighbor[1]);//Residuals
                                        //vec_ruv_w_scale_mag->val[num_grid_neighbor][3]=sqrt((dxy_neighbor[0]-apv_neighbor[0])*(dxy_neighbor[0]-apv_neighbor[0])+(dxy_neighbor[1]-apv_neighbor[1])*(dxy_neighbor[1]-apv_neighbor[1]));
                                        num_grid_neighbor++;
                                    }
                                }
                            }
                    
                            if(num_grid_neighbor>=thres_num_grid_neighbor)
                            {//Enough neighboring samples. perform the interpolation
                                
                                float w_min=1E+37;
                                float w_max=-1E+37;
                                float max_noi=1.0;
                                int32_t id_w_max=0;
                                int32_t id_w_min=0;
                                for(cnt_neighbor=0;cnt_neighbor<num_grid_neighbor;cnt_neighbor++)
                                {
                                    //Calculate the weight using dot product
                                    dxy[0]=vec_ruv_w_scale_mag->val[cnt_neighbor][0];
                                    dxy[1]=vec_ruv_w_scale_mag->val[cnt_neighbor][1];
                                    float mag_dxy=sqrt(dxy[0]*dxy[0]+dxy[1]*dxy[1]);
                                    float w_candidate=(dpe[0]*dxy[0]+dpe[1]*dxy[1])/(mag_dpe*mag_dxy);
                                    w_candidate=w_candidate>0?w_candidate:-w_candidate;
                                    
                                    //cut off if the inner product is less than the threshold
                                    //TODO: consider distance to determine the weight
                                    //Might be better to make LUT when calculating the ruv
                                    if(w_candidate >= thres_weight ) //dot product
                                    {
                                        //vec_ruv_w_scale_mag->val[cnt_neighbor][2]=w_candidate*w_candidate;//attenuate the weight by squaring it (trial code)
                                        vec_ruv_w_scale_mag->val[cnt_neighbor][2]=w_candidate;
                                        //check the min/max
                                        if(vec_ruv_w_scale_mag->val[cnt_neighbor][3]>w_max)
                                        {
                                            w_max=vec_ruv_w_scale_mag->val[cnt_neighbor][3];
                                            id_w_max=cnt_neighbor;
                                        }

                                        if(vec_ruv_w_scale_mag->val[cnt_neighbor][3]<w_min)
                                        {
                                            w_min=vec_ruv_w_scale_mag->val[cnt_neighbor][3];
                                            id_w_min=cnt_neighbor;
                                        }
                                    }
                                    else
                                    {
                                        vec_ruv_w_scale_mag->val[cnt_neighbor][2]=0.0;
                                    }
                                }
                                //cut out the maximum and minimum values
                                vec_ruv_w_scale_mag->val[id_w_max][2]=0.0;
                                vec_ruv_w_scale_mag->val[id_w_min][2]=0.0;


                                //debug code
                                if(cntu==15 && cntv==30)
                                {
                                    printf("dpe=[%f,%f]\n",dpe[0],dpe[1]);
                                    printf("vec_ruv_w_scale_mag:\n");
                                    GMA_float_print(vec_ruv_w_scale_mag);
                                }
                                //end of the debugging code

                                float sum_mag_times_w=0.0,sum_w=0.0,sum_w2=0.0,w2;
                                float sum_w_dp=0.0,sum_w_dpe=0.0,sum_noi=0.0;
                                for(cnt_neighbor=0;cnt_neighbor<num_grid_neighbor;cnt_neighbor++)
                                {
                                    //w2=1/(1+expf(-mag_dpe+5));
                                    w2=1/(1+expf(-vec_ruv_w_scale_mag->val[cnt_neighbor][5]+5))/max_noi;
                                    //printf("mag_dpe=%f, w2=%f\n",mag_dpe,w2);
                                    //sum_mag_times_w+=vec_ruv_w_scale_mag->val[cnt_neighbor][2]*vec_ruv_w_scale_mag->val[cnt_neighbor][3];
                                    sum_mag_times_w+=vec_ruv_w_scale_mag->val[cnt_neighbor][2]*vec_ruv_w_scale_mag->val[cnt_neighbor][3]*w2;
                                    sum_w+=vec_ruv_w_scale_mag->val[cnt_neighbor][2];
                                    //sum_w2+=vec_ruv_w_scale_mag->val[cnt_neighbor][2]*w2;
                                    

                                    //debug code
                                    if(cntu==15 && cntv==30)
                                    {
                                        printf("cnt_neighbor=%d, mag_dpe=%f, w2=%f\n",cnt_neighbor,mag_dpe,w2);
                                    }
                                    //end of the debugging code

                                    //test code - 3/12/2018 trial 1
                                    sum_w_dp+=vec_ruv_w_scale_mag->val[cnt_neighbor][2]*w2*vec_ruv_w_scale_mag->val[cnt_neighbor][4]/vec_ruv_w_scale_mag->val[cnt_neighbor][6];
                                    sum_w_dpe+=vec_ruv_w_scale_mag->val[cnt_neighbor][2]*w2*vec_ruv_w_scale_mag->val[cnt_neighbor][5]/vec_ruv_w_scale_mag->val[cnt_neighbor][6];
                                    sum_noi+=vec_ruv_w_scale_mag->val[cnt_neighbor][6];
                                    sum_w2+=vec_ruv_w_scale_mag->val[cnt_neighbor][2]*w2/vec_ruv_w_scale_mag->val[cnt_neighbor][6];


                                }

                                //debug code
                                if(cntu==15 && cntv==30)
                                {
                                    printf("sum_mag_times_w=%f\n",sum_mag_times_w);
                                    printf("sum_w=%f\n",sum_w);
                                    printf("sum_w2=%f\n",sum_w2);
                                }
                                    //end of the debugging code
                                
                                if(sum_w>=1.0)
                                {
                                    //float factor_mag=sum_mag_times_w/sum_w;
                                    //float factor_mag=sum_mag_times_w/sum_w2;
                                    //dx_buffer->val[cntv][cntu]=dpe[0]*factor_mag;
                                    //dy_buffer->val[cntv][cntu]=dpe[1]*factor_mag;
                                    
                                    //float factor_mag=1.0+sum_mag_times_w/sum_w/mag_dpe;
                                    //dx_buffer->val[cntv][cntu]=dpe[0]*factor_mag;
                                    //dy_buffer->val[cntv][cntu]=dpe[1]*factor_mag;

                                    //test code - 3/12/2018 trial 1
                                    float factor_mag=sum_w_dp/sum_w_dpe;
                                    dx_buffer->val[cntv][cntu]=dpe[0]*factor_mag;
                                    dy_buffer->val[cntv][cntu]=dpe[1]*factor_mag;
                                    mtrx_noi->val[cntv][cntu]=sum_noi/num_grid_neighbor+1;
                                    num_processed++;
                                }

                                


                            }//if(num_grid_neighbor>=5)
                        }//if(isnan(dpf_dx->val[cntv][cntu]) && isnan(dpf_dx->val[cntv][cntu]) && mvn_dp[cnt_grid]->nrows!=0)
                    }//for(cntu=0;cntu<dimx_vmap;cnty++)
                }//for(cntv=0;cntv<dimy_vmap;cntv++)
                GMA_float_destroy(vec_ruv_w_scale_mag);
                //}//#pragma omp parallel private....
                
                //put the interpolated values in the buffer to dpf_dx and dpf_dy
                for(cntv=0;cntv<dimy_vmap;cntv++)
                {
                    for(cntu=0;cntu<dimx_vmap;cntu++)
                    {
                        if(!isnan(dx_buffer->val[cntv][cntu]) && !isnan(dy_buffer->val[cntv][cntu]))
                        {
                            dpf_dx->val[cntv][cntu]=dx_buffer->val[cntv][cntu];
                            dpf_dy->val[cntv][cntu]=dy_buffer->val[cntv][cntu];
                            dx_buffer->val[cntv][cntu]=N_A_N;
                            dy_buffer->val[cntv][cntu]=N_A_N;
                        }
                    }
                }
                //sprintf(buffer_filename,"/Users/seongsu/Desktop/bad_vmaptar_example/dpf2_dx_%02d.GMA",NOI);
                //GMA_float_save(buffer_filename,dpf_dx);
                //sprintf(buffer_filename,"/Users/seongsu/Desktop/bad_vmaptar_example/dpf2_dy_%02d.GMA",NOI);
                //GMA_float_save(buffer_filename,dpf_dy);
            }//while(num_processed!=0)
    
            //scan through the dpf0 and see if there is any grid that does not have corresponding clusters
            num_grid_unprocessed=0;
            for(cntv=0;cntv<dimy_vmap;cntv++)
            {
                for(cntu=0;cntu<dimx_vmap;cntu++)
                {
                    cnt_grid=cntv*dimx_vmap+cntu;
                    if((isnan(dpf_dx->val[cntv][cntu]) || isnan(dpf_dy->val[cntv][cntu])) && mvn_dp[cnt_grid]->nrows!=0 )
                    {
                        num_grid_unprocessed++;
                    }
                }
            }
            printf("NOI=%d, weight threshold=%.2f, #remaining=%d\n",NOI,thres_weight,num_grid_unprocessed);

        }//while(num_grid_unprocessed!=0)
    }//for(thres_num_grid_neighbor=ruv_neighbor->nrows-1;thres_num_grid_neighbor>=3;thres_num_grid_neighbor--)

    //GMA_float_save("../MIMC2_C_GMA/dpf2_dx_interp.GMA",dpf_dx);
    //GMA_float_save("../MIMC2_C_GMA/dpf2_dy_interp.GMA",dpf_dy);


    //Smoothing the interpolated grids (dpf_dx and dpf_dy)
    printf("Smoothing the inteprolation results before finding the corresponding clusters\n");
    float sum_dx,sum_dy;
    float num_disp;
    int32_t du,dv;
    for(cntv=1;cntv<dimy_vmap-1;cntv++)
    {
        for(cntu=1;cntu<dimx_vmap-1;cntu++)
        {
            //if(dpf0->val[cntv][cntu]<0 && mvn_dp[cnt_grid]->nrows!=0)
            if(dpf0->val[cntv][cntu]<0 && !isnan(dpf_dx->val[cntv][cntu]+dpf_dy->val[cntv][cntu]))
            {
                //printf("Smoothing comes in! [cntu,cntv]=[%d,%d]\n",cntu,cntv);
                num_disp=0.0;
                sum_dx=0.0;
                sum_dy=0.0;
                for(dv=-1;dv<=1;dv++)
                {
                    for(du=-1;du<=1;du++)
                    {
                        if(!isnan(dpf_dx->val[cntv+dv][cntu+du]+dpf_dy->val[cntv+dv][cntu+du]))
                        {
                            sum_dx+=dpf_dx->val[cntv+dv][cntu+du];
                            sum_dy+=dpf_dy->val[cntv+dv][cntu+du];
                            num_disp=num_disp+1;
                        }
                        
                    }
                }
                dx_buffer->val[cntv][cntu]=sum_dx/num_disp;
                dy_buffer->val[cntv][cntu]=sum_dy/num_disp;
            }
            else
            {
                dx_buffer->val[cntv][cntu]=dpf_dx->val[cntv][cntu];
                dy_buffer->val[cntv][cntu]=dpf_dy->val[cntv][cntu];
            }

        }
    }

    for(cntv=1;cntv<dimy_vmap-1;cntv++)
    {
        for(cntu=1;cntu<dimx_vmap-1;cntu++)
        {
            dpf_dx->val[cntv][cntu]=dx_buffer->val[cntv][cntu];
            dpf_dy->val[cntv][cntu]=dy_buffer->val[cntv][cntu];
        }
    }

    printf("Smoothing complete\n");
    //sprintf(buffer_filename,"/Users/seongsu/Desktop/bad_vmaptar_example/dpf2_dx_%02d.GMA",NOI+1);
    //GMA_float_save(buffer_filename,dpf_dx);
    //sprintf(buffer_filename,"/Users/seongsu/Desktop/bad_vmaptar_example/dpf2_dy_%02d.GMA",NOI+1);
    //GMA_float_save(buffer_filename,dpf_dy);



    
    printf("dpf2 approximation successful. proceeding to final adjustment\n");

    //put mean cluster dps and IDs closest to the interpolated displacements
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            cnt_grid=cntv*dimx_vmap+cntu;
            if(dpf0->val[cntv][cntu]<0 && mvn_dp[cnt_grid]->nrows!=0)
            {
                sq_dist_closest=1E+37;
                id_cluster=0;
                for(cnt_cluster=0;cnt_cluster<mvn_dp[cnt_grid]->nrows;cnt_cluster++)
                {
                    dxy[0]=dpf_dx->val[cntv][cntu]-mvn_dp[cnt_grid]->val[cnt_cluster][0];
                    dxy[1]=dpf_dy->val[cntv][cntu]-mvn_dp[cnt_grid]->val[cnt_cluster][1];
                    sq_dist=dxy[0]*dxy[0]+dxy[1]*dxy[1];

                    if(sq_dist<sq_dist_closest)
                    {
                        sq_dist_closest=sq_dist;
                        id_cluster=cnt_cluster;
                    }
                }
                dpf0->val[cntv][cntu]=id_cluster;
                dpf_dx->val[cntv][cntu]=mvn_dp[cnt_grid]->val[id_cluster][0];
                dpf_dy->val[cntv][cntu]=mvn_dp[cnt_grid]->val[id_cluster][1];
            }
        }
    }

    //GMA_float_save("/Users/seongsu/Desktop/bad_vmaptar_example/dpf2_dx_final.GMA",dpf_dx);
    //GMA_float_save("/Users/seongsu/Desktop/bad_vmaptar_example/dpf2_dy_final.GMA",dpf_dy);
    //GMA_int32_save("../MIMC2_C_GMA/dpf2_cluster.GMA",dpf0);

    printf("dpf1 adjustment successful. finalizing the process\n");

    //de-allocate array
    GMA_float_destroy(dx_buffer);
    GMA_float_destroy(dy_buffer);
    //GMA_float_destroy(vec_ruv_w_scale_mag);
}

void get_dpf1_original(GMA_int32 *dpf0, GMA_float *dpf_dx, GMA_float *dpf_dy, GMA_int32 *ruv_neighbor, GMA_float **mvn_dp, GMA_double *xyuvav)  //get the initial dp for pseudosmoothing
{
    int32_t cntu,cntv,cnt_grid,cnt_neighbor,cnt_cluster,id_cluster;
    int32_t NOI,num_processed,num_grid_neighbor,num_grid_unprocessed;
    int32_t uvoi[2];
    GMA_float *dx_buffer, *dy_buffer;
    GMA_float *vec_ruv_w_scale_mag;
    float dpe[2],dxy_neighbor[2],dxy[2],apv_neighbor[2];
    float N_A_N=sqrt(-1.0);
    float sum_mag,sq_dist,sq_dist_closest;
    char buffer_filename[2048];

    
    dx_buffer=GMA_float_create(dimy_vmap,dimx_vmap);
    dy_buffer=GMA_float_create(dimy_vmap,dimx_vmap);
    //vec_ruv_w_scale_mag=GMA_float_create(ruv_neighbor->nrows,4);
    
    //reconstruct dpf_dx and dpf_dy. Also initialize the dx_buffer and dy_buffer
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            cnt_grid=cntv*dimx_vmap+cntu;
            if(dpf0->val[cntv][cntu]>=0)
            {
                dpf_dx->val[cntv][cntu]=mvn_dp[cnt_grid]->val[dpf0->val[cntv][cntu]][0];
                dpf_dy->val[cntv][cntu]=mvn_dp[cnt_grid]->val[dpf0->val[cntv][cntu]][1];
            }
            else
            {
                dpf_dx->val[cntv][cntu]=N_A_N;
                dpf_dy->val[cntv][cntu]=N_A_N;
            }
            dx_buffer->val[cntv][cntu]=N_A_N;
            dy_buffer->val[cntv][cntu]=N_A_N;
        }
    }

    //debugging code
    GMA_float_save("/Users/seongsu/Desktop/bad_vmaptar_example/dpf0_gx.gma",dpf_dx);
    GMA_float_save("/Users/seongsu/Desktop/bad_vmaptar_example/dpf0_gy.gma",dpf_dy);


    NOI=0;
    num_grid_unprocessed=1;
    float thres_weight=1.0;
    while(num_grid_unprocessed!=0 && thres_weight>=0.5)
    {
        thres_weight-=0.02;
        num_processed=1;
        while(num_processed!=0)
        {   
            NOI++;
            num_processed=0;
            //#pragma omp parallel private(cnt_grid,id_cluster,num_grid_neighbor,uvoi, vec_ruv_w_scale_mag, dpe, dxy_neighbor, dxy, sq_dist, sq_dist_closest, sum_mag) \
                                            shared(dpf0, dpf_dx, dpf_dy, ruv_neighbor, mvn_dp, xyuvav, dx_buffer, dy_buffer, N_A_N)
           // {
            vec_ruv_w_scale_mag=GMA_float_create(ruv_neighbor->nrows,4);    
            //#pragma omp for schedule(dynamic) collapse(2) 
            for(cntv=0;cntv<dimy_vmap;cntv++)
            {   //loop around the array and process if the visited grid needs to
                for(cntu=0;cntu<dimx_vmap;cntu++)
                {
                    cnt_grid=cntv*dimx_vmap+cntu;
                    if(isnan(dpf_dx->val[cntv][cntu]) && isnan(dpf_dy->val[cntv][cntu]) && mvn_dp[cnt_grid]->nrows!=0)
                    {      //dx validity check             dy validity check                    nomatching check
                        sum_mag=0;
                        num_grid_neighbor=0;
                        dpe[0]=xyuvav->val[cnt_grid][4]/365*dt/param_mimc2.mpp;
                        dpe[1]=-xyuvav->val[cnt_grid][5]/365*dt/param_mimc2.mpp;
                        for(cnt_neighbor=0;cnt_neighbor<ruv_neighbor->nrows;cnt_neighbor++)
                        {
                            uvoi[0]=cntu+ruv_neighbor->val[cnt_neighbor][0];
                            uvoi[1]=cntv+ruv_neighbor->val[cnt_neighbor][1];
                            if(uvoi[0]>=0 && uvoi[0]<dimx_vmap && uvoi[1]>=0 && uvoi[1]<dimy_vmap)//boundary check
                            {
                                dxy_neighbor[0]=dpf_dx->val[uvoi[1]][uvoi[0]];
                                dxy_neighbor[1]=dpf_dy->val[uvoi[1]][uvoi[0]];
                                
                                
                                
                                if(!isnan(dxy_neighbor[0]) && !isnan(dxy_neighbor[1]))
                                {
                                    apv_neighbor[0]=(float)(xyuvav->val[uvoi[1]*dimx_vmap+uvoi[0]][4])/365*dt/param_mimc2.mpp;
                                    apv_neighbor[1]=-(float)(xyuvav->val[uvoi[1]*dimx_vmap+uvoi[0]][5])/365*dt/param_mimc2.mpp;
                                
                                    vec_ruv_w_scale_mag->val[num_grid_neighbor][0]=(float)ruv_neighbor->val[cnt_neighbor][0];//relative u coordinate
                                    vec_ruv_w_scale_mag->val[num_grid_neighbor][1]=(float)ruv_neighbor->val[cnt_neighbor][1];//relative v coordinate
                                    //vec_ruv_w_scale_mag->val[num_grid_neighbor][3]=sqrt(dxy_neighbor[0]*dxy_neighbor[0]+dxy_neighbor[1]*dxy_neighbor[1]);//magnitude
                                    vec_ruv_w_scale_mag->val[num_grid_neighbor][3]=sqrt((dxy_neighbor[0]*dxy_neighbor[0]+dxy_neighbor[1]*dxy_neighbor[1])/(apv_neighbor[0]*apv_neighbor[0]+apv_neighbor[1]*apv_neighbor[1]));//scale factor
                                    //49,57
                                    if(uvoi[0]==48 && uvoi[0]==56)
                                    {
                                        printf("%f\t",sqrt((dxy_neighbor[0]*dxy_neighbor[0]+dxy_neighbor[1]*dxy_neighbor[1])/(apv_neighbor[0]*apv_neighbor[0]+apv_neighbor[1]*apv_neighbor[1])));


                                    }
                                    num_grid_neighbor++;
                                }
                            }
                        }
                    
                        if(num_grid_neighbor>=5)
                        {//Enough neighboring samples. perform the interpolation
                            //normalize dpe
                            float mag_dpe=sqrt(dpe[0]*dpe[0]+dpe[1]*dpe[1]);
                            dpe[0]/=mag_dpe;
                            dpe[1]/=mag_dpe;
                            float w_min=1E+37;
                            float w_max=-1E+37;
                            int32_t id_w_max=0;
                            int32_t id_w_min=0;
                            for(cnt_neighbor=0;cnt_neighbor<num_grid_neighbor;cnt_neighbor++)
                            {
                                //calculate inner product
                                dxy[0]=vec_ruv_w_scale_mag->val[cnt_neighbor][0];
                                dxy[1]=vec_ruv_w_scale_mag->val[cnt_neighbor][1];
                                float mag_dxy=sqrt(dxy[0]*dxy[0]+dxy[1]*dxy[1]);
                                dxy[0]/=mag_dxy;
                                dxy[1]/=mag_dxy;

                                //Calculate the weight using tot product
                                //cut off if the inner product is less than the threshold (0.7)
                                //TODO :make the threshold  (0.7) to be adjustable
                                float w_candidate=dpe[0]*dxy[0]+dpe[1]*dxy[1];
                                if(w_candidate<0)
                                {
                                    w_candidate=-w_candidate;
                                }
                            
                                if(w_candidate >= thres_weight ) //dot product
                                {
                                    vec_ruv_w_scale_mag->val[cnt_neighbor][2]=w_candidate;
                                    
                                    //check the min/max
                                    if(vec_ruv_w_scale_mag->val[cnt_neighbor][3]>w_max)
                                    {
                                        w_max=vec_ruv_w_scale_mag->val[cnt_neighbor][3];
                                        id_w_max=cnt_neighbor;
                                    }

                                    if(vec_ruv_w_scale_mag->val[cnt_neighbor][3]<w_min)
                                    {
                                        w_min=vec_ruv_w_scale_mag->val[cnt_neighbor][3];
                                        id_w_min=cnt_neighbor;
                                    }

                                }
                                else
                                {
                                    vec_ruv_w_scale_mag->val[cnt_neighbor][2]=0.0;
                                }
                            }
                            //cut out the maximum and minimum values
                            vec_ruv_w_scale_mag->val[id_w_max][2]=0.0;
                            vec_ruv_w_scale_mag->val[id_w_min][2]=0.0;

                            float sum_mag_times_w=0,sum_w=0;
                            for(cnt_neighbor=0;cnt_neighbor<num_grid_neighbor;cnt_neighbor++)
                            {
                                sum_mag_times_w+=vec_ruv_w_scale_mag->val[cnt_neighbor][2]*vec_ruv_w_scale_mag->val[cnt_neighbor][3];
                                sum_w+=vec_ruv_w_scale_mag->val[cnt_neighbor][2];
                            }

                            if(sum_w>=1.0)
                            {
                                //float mag_dpe=
                                float factor_mag=sum_mag_times_w/sum_w;
                                dx_buffer->val[cntv][cntu]=dpe[0]*factor_mag*mag_dpe;
                                dy_buffer->val[cntv][cntu]=dpe[1]*factor_mag*mag_dpe;
                                num_processed++;
                            }
                        }//if(num_grid_neighbor>=5)
                    }//if(isnan(dpf_dx->val[cntv][cntu]) && isnan(dpf_dx->val[cntv][cntu]) && mvn_dp[cnt_grid]->nrows!=0)
                }//for(cntu=0;cntu<dimx_vmap;cnty++)
            }//for(cntv=0;cntv<dimy_vmap;cntv++)
            GMA_float_destroy(vec_ruv_w_scale_mag);
            //}//#pragma omp parallel private....
            //put the interpolated values in the buffer to dpf_dx and dpf_dy
            for(cntv=0;cntv<dimy_vmap;cntv++)
            {
                for(cntu=0;cntu<dimx_vmap;cntu++)
                {
                    if(!isnan(dx_buffer->val[cntv][cntu]) && !isnan(dy_buffer->val[cntv][cntu]))
                    {
                        dpf_dx->val[cntv][cntu]=dx_buffer->val[cntv][cntu];
                        dpf_dy->val[cntv][cntu]=dy_buffer->val[cntv][cntu];
                        dx_buffer->val[cntv][cntu]=N_A_N;
                        dy_buffer->val[cntv][cntu]=N_A_N;
                    }
                }
            }
            sprintf(buffer_filename,"/Users/seongsu/Desktop/bad_vmaptar_example/dpf2_dx_%02d.GMA",NOI);
            //printf("%s\n",buffer_filename);
            GMA_float_save(buffer_filename,dpf_dx);
            sprintf(buffer_filename,"/Users/seongsu/Desktop/bad_vmaptar_example/dpf2_dy_%02d.GMA",NOI);
            //printf("%s\n",buffer_filename);
            GMA_float_save(buffer_filename,dpf_dy);
        }//while(num_processed!=0)
    
        //scan through the dpf0 and see if there is any grid that does not have corresponding clusters
        num_grid_unprocessed=0;
        for(cntv=0;cntv<dimy_vmap;cntv++)
        {
            for(cntu=0;cntu<dimx_vmap;cntu++)
            {
                cnt_grid=cntv*dimx_vmap+cntu;
                if((isnan(dpf_dx->val[cntv][cntu]) || isnan(dpf_dy->val[cntv][cntu])) && mvn_dp[cnt_grid]->nrows!=0 )
                {
                    num_grid_unprocessed++;
                }
            }
        }
        printf("NOI=%d, weight threshold=%.2f, #remaining=%d\n",NOI,thres_weight,num_grid_unprocessed);

    }//while(num_grid_unprocessed!=0)

    //GMA_float_save("../MIMC2_C_GMA/dpf2_dx_interp.GMA",dpf_dx);
    //GMA_float_save("../MIMC2_C_GMA/dpf2_dy_interp.GMA",dpf_dy);
    

    

    printf("dpf2 approximation successful. proceeding to final adjustment\n");

    //put mean cluster dps and IDs closest to the interpolated displacements
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            cnt_grid=cntv*dimx_vmap+cntu;
            if(dpf0->val[cntv][cntu]<0 && mvn_dp[cnt_grid]->nrows!=0)
            {
                sq_dist_closest=1E+37;
                id_cluster=0;
                for(cnt_cluster=0;cnt_cluster<mvn_dp[cnt_grid]->nrows;cnt_cluster++)
                {
                    dxy[0]=dpf_dx->val[cntv][cntu]-mvn_dp[cnt_grid]->val[cnt_cluster][0];
                    dxy[1]=dpf_dy->val[cntv][cntu]-mvn_dp[cnt_grid]->val[cnt_cluster][1];
                    sq_dist=dxy[0]*dxy[0]+dxy[1]*dxy[1];

                    if(sq_dist<sq_dist_closest)
                    {
                        sq_dist_closest=sq_dist;
                        id_cluster=cnt_cluster;
                    }
                }
                dpf0->val[cntv][cntu]=id_cluster;
                dpf_dx->val[cntv][cntu]=mvn_dp[cnt_grid]->val[id_cluster][0];
                dpf_dy->val[cntv][cntu]=mvn_dp[cnt_grid]->val[id_cluster][1];
            }
        }
    }

    //GMA_float_save("../MIMC2_C_GMA/dpf2_dx_final.GMA",dpf_dx);
    //GMA_float_save("../MIMC2_C_GMA/dpf2_dy_final.GMA",dpf_dy);
    //GMA_int32_save("../MIMC2_C_GMA/dpf2_cluster.GMA",dpf0);

    printf("dpf1 adjustment successful. finalizing the process\n");

    //de-allocate array
    GMA_float_destroy(dx_buffer);
    GMA_float_destroy(dy_buffer);
    //GMA_float_destroy(vec_ruv_w_scale_mag);
}

void get_dpf_pseudosmoothing(GMA_int32 *dpf,GMA_float *dpf_dx, GMA_float *dpf_dy,GMA_int32 *ruv_neighbor,GMA_float **mvn_dp, GMA_double *xyuvav)
{
    GMA_uint8 *mask_investigate, *mask_investigate_next, **mask_investigate_layer;
    GMA_uint8 **stack_mask_investigate;
    GMA_int32 *uv_neighbor;//,*uvoi;
    GMA_double *uvoi;
    GMA_double *w, *duv_neighbor, *duv_interp;
    int32_t NOI=0, max_iter=100;
    int32_t cntu,cntv,cntn,cntc;
    int32_t idx_vec_corr, id_cluster;
    int32_t num_neighbor;
    double vxyoi[2],eigvel[2],ITM[4],duv_cluster[2],duv_grid[2],duv_candidate[2];
    char flag_any_modification, flag_fluctuation;
    double sq_dist, sq_dist_min;
    int32_t num_cluster,id_cluster_closest;
    float N_A_N=sqrt(-1.0);

    eigvel[0]=1500.0/300.0;
    eigvel[1]=eigvel[0]/3.0;
    
    //determine the initial values of mask_investigate and mask_investigate_next
    printf("Initializing the layer array\n");
    mask_investigate_layer=(GMA_uint8**)malloc(2*sizeof(GMA_uint8*));
    mask_investigate_layer[0]=GMA_uint8_create(dimy_vmap,dimx_vmap);
    mask_investigate_layer[1]=GMA_uint8_create(dimy_vmap,dimx_vmap);
    mask_investigate=mask_investigate_layer[0];
    mask_investigate_next=mask_investigate_layer[1];
    GMA_float **dxy_ps_buffer;
    GMA_int32 *dpf_ps_buffer;
    
    uv_neighbor=GMA_int32_create(ruv_neighbor->nrows,ruv_neighbor->ncols);
    duv_neighbor=GMA_double_create(ruv_neighbor->nrows,ruv_neighbor->ncols);
    uvoi=GMA_double_create(1,2);
    w=GMA_double_create(ruv_neighbor->nrows,1);
    duv_interp=GMA_double_create(1,2);
    
    dxy_ps_buffer=malloc(2*sizeof(GMA_float*));
    dxy_ps_buffer[0]=GMA_float_create(dimy_vmap,dimx_vmap);
    dxy_ps_buffer[1]=GMA_float_create(dimy_vmap,dimx_vmap);
    dpf_ps_buffer=GMA_int32_create(dimy_vmap,dimx_vmap);

    //TODO: consider parallelizing the for loop below
    printf("Determining the grids to investigate, while initializing dxy_ps_buffer\n");
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {   
            idx_vec_corr=cntv*dimx_vmap+cntu;
            id_cluster=dpf->val[cntv][cntu];
            if(id_cluster<0)
            {
                mask_investigate->val[cntv][cntu]=0;
            }
            else
            {
                if(mvn_dp[idx_vec_corr]->val[id_cluster][4]>=0.6)
                {
                    
                    mask_investigate->val[cntv][cntu]=0;
                }
                else
                {
                    mask_investigate->val[cntv][cntu]=1;
                }
            }
            dxy_ps_buffer[0]->val[cntv][cntu]=N_A_N;
            dxy_ps_buffer[1]->val[cntv][cntu]=N_A_N;
            dpf_ps_buffer->val[cntv][cntu]=-1;
        }
    }

    //GMA_uint8_save("../MIMC2_C_GMA/mask_investigate.GMA",mask_investigate);
    
    printf("Allocating mask stack\n");
    //allocate stack_mask_investigate
    stack_mask_investigate=(GMA_uint8 **)malloc(max_iter*sizeof(GMA_uint8*));

    //copy the initial mask_investigate
    stack_mask_investigate[0]=GMA_uint8_create(dimy_vmap,dimx_vmap);
    //TODO: consider parallelizing the for loop below
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            stack_mask_investigate[0]->val[cntv][cntu]=mask_investigate->val[cntv][cntu];
        }
    }

    printf("Initiating while loop\n");
    flag_any_modification=1;
    flag_fluctuation=0;
    while(NOI<=100 && flag_any_modification && ~flag_fluctuation)
    {
        mask_investigate=mask_investigate_layer[NOI%2];
        mask_investigate_next=mask_investigate_layer[(NOI+1)%2];

        flag_any_modification=0;
        NOI++;
        //clean up mask_investigate_next
        //TODO: consider parallelizing the for loop below
        for(cntv=0;cntv<dimy_vmap;cntv++)
        {
            for(cntu=0;cntu<dimx_vmap;cntu++)
            {
                mask_investigate_next->val[cntv][cntu]=0;
            }
        }

        //loop through the mask_investigate
        for(cntv=0;cntv<dimy_vmap;cntv++)
        {
            for(cntu=0;cntu<dimx_vmap;cntu++)
            {
                if(mask_investigate->val[cntv][cntu]==0)
                {
                    //The case that the visited grid does not need to be investigated
                    continue;
                }
                else //The case that the visited grid NEEDS to be investigated
                {
                    //extract the u,v, du and dv
                    num_neighbor=0;
                    for(cntn=0;cntn<ruv_neighbor->nrows;cntn++)
                    {
                        int32_t uvoi_neighbor[2];
                        uvoi_neighbor[0]=cntu+ruv_neighbor->val[cntn][0];
                        uvoi_neighbor[1]=cntv+ruv_neighbor->val[cntn][1];
                        
                        if(uvoi_neighbor[0]>=0 && uvoi_neighbor[0]<dimx_vmap && //boundary check
                           uvoi_neighbor[1]>=0 && uvoi_neighbor[1]<dimy_vmap &&
                           !isnan(dpf_dx->val[uvoi_neighbor[1]][uvoi_neighbor[0]]) &&
                           !isnan(dpf_dy->val[uvoi_neighbor[1]][uvoi_neighbor[0]])) //validity check of displacement matrix
                        {
                            uv_neighbor->val[num_neighbor][0]=ruv_neighbor->val[cntn][0];
                            uv_neighbor->val[num_neighbor][1]=ruv_neighbor->val[cntn][1];
                            id_cluster=dpf->val[cntv+ruv_neighbor->val[cntn][1]][cntu+ruv_neighbor->val[cntn][0]];
                            duv_neighbor->val[num_neighbor][0]=dpf_dx->val[uvoi_neighbor[1]][uvoi_neighbor[0]];
                            duv_neighbor->val[num_neighbor][1]=dpf_dy->val[uvoi_neighbor[1]][uvoi_neighbor[0]];
                            num_neighbor++;
                        }
                    }
                    uv_neighbor->nrows=num_neighbor;
                    duv_neighbor->nrows=num_neighbor;
                    w->nrows=num_neighbor;
                    
                    if(num_neighbor>=10) // Enough number of neighboring grids
                    {
                        //printf("Enough number of neighboring displacements found: %d>=10\n",num_neighbor);
                        uvoi->val[0][0]=0.0;
                        uvoi->val[0][1]=0.0;

                        //calculate weight
                        idx_vec_corr=cntv*dimx_vmap+cntu;

                        //retrieve a priori from xyuvav
                        vxyoi[0]=xyuvav->val[idx_vec_corr][4];
                        vxyoi[1]=xyuvav->val[idx_vec_corr][5];

                        ITM[0]=(eigvel[1]*vxyoi[0]*vxyoi[0]+eigvel[0]*vxyoi[1]*vxyoi[1])/((eigvel[0]*eigvel[1])*(vxyoi[0]*vxyoi[0]+vxyoi[1]*vxyoi[1]));
                        ITM[1]=((eigvel[0]-eigvel[1])*vxyoi[0]*vxyoi[1])/((eigvel[0]*eigvel[1])*(vxyoi[0]*vxyoi[0]+vxyoi[1]*vxyoi[1]));
                        ITM[3]=(eigvel[1]*vxyoi[1]*vxyoi[1] +eigvel[0]*vxyoi[0]*vxyoi[0])/((eigvel[0]*eigvel[1])*(vxyoi[0]*vxyoi[0]+vxyoi[1]*vxyoi[1]));

                        
                        for(cntn=0;cntn<num_neighbor;cntn++)
                        {
                            w->val[cntn][0]=exp(-(ITM[0]*uv_neighbor->val[cntn][0]*uv_neighbor->val[cntn][0]
                                                        +2*ITM[1]*uv_neighbor->val[cntn][0]*uv_neighbor->val[cntn][1]
                                                        +ITM[3]*uv_neighbor->val[cntn][1]*uv_neighbor->val[cntn][1]));
                                                        //TODO: adjust the coefficients
                                                        //TODO: Think about making use of the symmetric properties of the weight vector
                        }
                        //End calculating weight

                        //calculate dui and dvi
                        quadfit2(uv_neighbor,duv_neighbor,w,uvoi,duv_interp);
                        //compare duv_interp with the mean values of the corresponding clusters
                        id_cluster=dpf->val[cntv][cntu];
                        num_cluster=mvn_dp[idx_vec_corr]->nrows;
                        
                        sq_dist_min=1E+37;
                        int8_t flag_update_grid=0;
                        for(cntc=0;cntc<num_cluster;cntc++)
                        {
                            duv_cluster[0]=mvn_dp[idx_vec_corr]->val[cntc][0];
                            duv_cluster[1]=mvn_dp[idx_vec_corr]->val[cntc][1];
                            sq_dist=(duv_interp->val[0][0]-duv_cluster[0])*(duv_interp->val[0][0]-duv_cluster[0])+(duv_interp->val[0][1]-duv_cluster[1])*(duv_interp->val[0][1]-duv_cluster[1]);
                            //NOTE: data type casting is expected at the line above
                            
                            if(sq_dist<sq_dist_min)
                            {
                                flag_update_grid=1;
                                sq_dist_min=sq_dist;
                                id_cluster_closest=cntc;
                            }
                        }

                        if(flag_update_grid)
                        {
                            duv_grid[0]=mvn_dp[idx_vec_corr]->val[id_cluster][0];
                            duv_grid[1]=mvn_dp[idx_vec_corr]->val[id_cluster][1];
                            duv_candidate[0]=mvn_dp[idx_vec_corr]->val[id_cluster_closest][0];
                            duv_candidate[1]=mvn_dp[idx_vec_corr]->val[id_cluster_closest][1];
                        }

                        if((duv_grid[0]-duv_candidate[0])*(duv_grid[0]-duv_candidate[0]) +(duv_grid[1]-duv_candidate[1])*(duv_grid[1]-duv_candidate[1]) < 0.0001 )  //TODO: make the threshold adjustable
                        {
                            continue;
                        }
                        else
                        {
                            dxy_ps_buffer[0]->val[cntv][cntu]=duv_candidate[0];
                            dxy_ps_buffer[1]->val[cntv][cntu]=duv_candidate[1];
                            dpf_ps_buffer->val[cntv][cntu]=id_cluster_closest;

                            flag_any_modification=1;

                            //mark the neighboring location to investigate in the next iteration
                            for(cntn=0;cntn<num_neighbor;cntn++)
                            {
                                if(stack_mask_investigate[0]->val[cntv+uv_neighbor->val[cntn][1]][cntu+uv_neighbor->val[cntn][0]])
                                //NOTE: The 1st layer of this stack indicates dpf0 and non-matching results grids
                                {
                                    mask_investigate_next->val[cntv+uv_neighbor->val[cntn][1]][cntu+uv_neighbor->val[cntn][0]]=1;
                                }
                            }
                        }
                    }
                }
            }//for(cntu=0;cntu<dimx_vmap;cntu++)
        }//for(cntv=0;cntv<dimy_vmap;cntv++)

        //update and reset dpf_dx, dpf_dy and dpf
        for(cntu=0;cntu<dimx_vmap;cntu++)
        {
            for(cntv=0;cntv<dimy_vmap;cntv++)
            {
                if(dpf_ps_buffer->val[cntv][cntu]>=0)
                {
                    dpf_dx->val[cntv][cntu]=dxy_ps_buffer[0]->val[cntv][cntu];
                    dpf_dy->val[cntv][cntu]=dxy_ps_buffer[1]->val[cntv][cntu];
                    dpf->val[cntv][cntu]=dpf_ps_buffer->val[cntv][cntu];
                    
                    dxy_ps_buffer[0]->val[cntv][cntu]=N_A_N;
                    dxy_ps_buffer[1]->val[cntv][cntu]=N_A_N;
                    dpf_ps_buffer->val[cntv][cntu]=-1;
                }
            }
        }

        //check the fluctuation
        
        for(cntn=NOI-1;cntn>=0;cntn--)
        {
            char flag_difference_found=0;
            for(cntv=0;cntv<dimy_vmap;cntv++)
            {
                for(cntu=0;cntu<dimx_vmap;cntu++)
                {
                    if(stack_mask_investigate[cntn]->val[cntv][cntu]!=mask_investigate_next->val[cntv][cntu])
                    {
                        flag_difference_found=1;
                        break;
                    }
                }
                if(flag_difference_found)
                {
                    break;
                }
            }
            if(flag_difference_found==0)
            {
                printf("Fluctuation detected: NOI=%d, cntn=%d\n",NOI,cntn);
                flag_fluctuation=1;
                break;
            }
        }

        if(flag_fluctuation)
        {
            printf("Fluctuation was detected - Terminating the pseudosmoothing.\n");
            NOI--;
            break;
        }
        else
        {
            //stack up mask_investigate_next to stack
            stack_mask_investigate[NOI]=GMA_uint8_create(dimy_vmap,dimx_vmap);
            int32_t num_grid_process=0;
            for(cntv=0;cntv<dimy_vmap;cntv++)
            {
                for(cntu=0;cntu<dimx_vmap;cntu++)
                {
                    stack_mask_investigate[NOI]->val[cntv][cntu]=mask_investigate_next->val[cntv][cntu];
                    if(mask_investigate_next->val[cntv][cntu])
                    {
                        num_grid_process++;
                    }
                }
            }
            printf("NOI=%d,  # grids=%d\n",NOI,num_grid_process);
        }

    }//while(NOI<100 && flag_any_modification && ~flag_fluctuation)


    //de-allocate the arrays
    for(cntn=0;cntn<=NOI;cntn++)
    {
        GMA_uint8_destroy(stack_mask_investigate[cntn]);
    }
    free(stack_mask_investigate);
    GMA_uint8_destroy(mask_investigate_layer[0]);
    GMA_uint8_destroy(mask_investigate_layer[1]);
    uv_neighbor->nrows=ruv_neighbor->nrows; //revert the original size of the array that have altered during the while loop
    duv_neighbor->nrows=ruv_neighbor->nrows;
    w->nrows=ruv_neighbor->nrows;
    GMA_int32_destroy(uv_neighbor);
    GMA_double_destroy(uvoi);
    GMA_double_destroy(duv_neighbor);
    GMA_double_destroy(w);
    GMA_double_destroy(duv_interp);

    GMA_float_destroy(dxy_ps_buffer[0]);
    GMA_float_destroy(dxy_ps_buffer[1]);
    free(dxy_ps_buffer);

}

void quadfit2(GMA_int32 *xy, GMA_double *z, GMA_double *w, GMA_double *xyi, GMA_double *out)
{
    int32_t num_obs=xy->nrows;
    GMA_double *A=GMA_double_create(num_obs,6);
    GMA_double *N=GMA_double_create(6,6);
    GMA_double *IN=GMA_double_create(6,6);
    double terms_quadratic[6];
    double coeff_quadratic[6];

    GMA_double *atwb=GMA_double_create(6,1);
    
    double x,y;

    int32_t cnt_row,cnt_col,cnt_obs,cnt_outcol;

    
    //calculate A
    for(cnt_row=0;cnt_row<num_obs;cnt_row++)
    {
        x=(double)xy->val[cnt_row][0];
        y=(double)xy->val[cnt_row][1];
        A->val[cnt_row][0]=x*x;
        A->val[cnt_row][1]=x*y;
        A->val[cnt_row][2]=y*y;
        A->val[cnt_row][3]=x;
        A->val[cnt_row][4]=y;
        A->val[cnt_row][5]=1;
    }

    //calculate N
    for(cnt_row=0;cnt_row<6;cnt_row++)  //TODO: consider loop unrolling
    {
        for(cnt_col=0;cnt_col<6;cnt_col++)
        {
            N->val[cnt_row][cnt_col]=0;
            for(cnt_obs=0;cnt_obs<num_obs;cnt_obs++)
            {
                N->val[cnt_row][cnt_col]+=A->val[cnt_obs][cnt_row]*w->val[cnt_obs][0]*A->val[cnt_obs][cnt_col];
            }
        }
    }
    
    //calculate IN
    //printf("Calculating inverse of normal matrix\n");
    GMA_double_inv(N,IN);

    //printf("Calculating coefficient\n");
    for(cnt_outcol=0;cnt_outcol<out->ncols;cnt_outcol++)
    {
        //calculate atwb
        for(cnt_row=0;cnt_row<6;cnt_row++)
        {
            atwb->val[cnt_row][0]=0;
            for(cnt_obs=0;cnt_obs<num_obs;cnt_obs++)
            {
                atwb->val[cnt_row][0]+=A->val[cnt_obs][cnt_row]*w->val[cnt_obs][0]*z->val[cnt_obs][cnt_outcol];
            }
        }
        
        //calculate coeff_quadratic
        //NI: 6 by 6
        //atwb: 6 by 1
        for(cnt_row=0;cnt_row<6;cnt_row++)
        {
            coeff_quadratic[cnt_row]=0;
            for(cnt_col=0;cnt_col<6;cnt_col++)
            {
                //coeff_quadratic->val[cnt_row][0]+=IN->val[cnt_row][cnt_col]*atwb->val[cnt_col][0];
                coeff_quadratic[cnt_row]+=IN->val[cnt_row][cnt_col]*atwb->val[cnt_col][0];
            }
        }

        //calculate the interpolated values using the calculated coefficients

        for(cnt_row=0;cnt_row<out->nrows;cnt_row++)
        {
            terms_quadratic[0]=xyi->val[cnt_row][0]*xyi->val[cnt_row][0];
            terms_quadratic[1]=xyi->val[cnt_row][0]*xyi->val[cnt_row][1];
            terms_quadratic[2]=xyi->val[cnt_row][1]*xyi->val[cnt_row][1];
            terms_quadratic[3]=xyi->val[cnt_row][0];
            terms_quadratic[4]=xyi->val[cnt_row][1];
            terms_quadratic[5]=1;

            out->val[cnt_row][cnt_outcol]=0;
            for(cnt_col=0;cnt_col<6;cnt_col++)
            {
                out->val[cnt_row][cnt_outcol]+=terms_quadratic[cnt_col]*coeff_quadratic[cnt_col];
            }
        }
    }

    GMA_double_destroy(A);
    GMA_double_destroy(N);
    GMA_double_destroy(IN);
    GMA_double_destroy(atwb);
}


void GMA_double_mul(GMA_double *a, GMA_double *b, GMA_double *out)
{
    int32_t cnt1,cnt2,cnt3;
    for(cnt1=0;cnt1<a->nrows;cnt1++)  //TODO: consider loop unrolling
    //TODO: consider applying OMP multithreading 
    {
        for(cnt2=0;cnt2<b->ncols;cnt2++)
        {
            out->val[cnt1][cnt2]=0;
            for(cnt3=0;cnt3<a->ncols;cnt3++)
            {
                out->val[cnt1][cnt2]+=a->val[cnt1][cnt3]*b->val[cnt3][cnt2];
            }
        }
    }
}


void GMA_double_inv(GMA_double *a, GMA_double *I)
{
    int32_t cnt1,cnt2,cnt3;
    double pivot,coeff;
    
    GMA_double *b=GMA_double_create(a->nrows,a->ncols);
    //duplicate the matrix a
    for(cnt1=0;cnt1<a->nrows;cnt1++)
    {
        for(cnt2=0;cnt2<a->ncols;cnt2++)
        {
            b->val[cnt1][cnt2]=a->val[cnt1][cnt2];
        }
    }

    //initialize I
    for(cnt1=0;cnt1<b->nrows;cnt1++)
    {
        for(cnt2=0;cnt2<b->nrows;cnt2++)
        {
            I->val[cnt1][cnt2]=(cnt1==cnt2)?1:0;
        }
    }
    
    //Gaussian elimation - forward
    for(cnt1=0;cnt1<b->nrows-1;cnt1++)
    {
        pivot=b->val[cnt1][cnt1];
        for(cnt2=cnt1+1;cnt2<b->ncols;cnt2++)
        {
            coeff=b->val[cnt2][cnt1]/pivot;
            for(cnt3=0;cnt3<b->nrows;cnt3++)
            {
                b->val[cnt2][cnt3]-=b->val[cnt1][cnt3]*coeff;
                I->val[cnt2][cnt3]-=I->val[cnt1][cnt3]*coeff;
            }
        }
    }

    //Backward Elimination
    for(cnt1=b->nrows-1;cnt1>=0;cnt1--)
    {
        pivot=b->val[cnt1][cnt1];
        for(cnt2=cnt1-1;cnt2>=0;cnt2--)
        {
            coeff=b->val[cnt2][cnt1]/pivot;
            for(cnt3=b->nrows-1;cnt3>=0;cnt3--)
            {
                b->val[cnt2][cnt3]-=b->val[cnt1][cnt3]*coeff;
                I->val[cnt2][cnt3]-=I->val[cnt1][cnt3]*coeff;
            }
        }
    }

    //scaling
    for(cnt1=0;cnt1<b->nrows;cnt1++)
    {
        for(cnt2=0;cnt2<b->nrows;cnt2++)
        {
            I->val[cnt1][cnt2]/=b->val[cnt1][cnt1];
        }
        
    }

    GMA_double_destroy(b);

}

void convert_dpf_to_vxy_exy_qual(GMA_int32 *dpf, GMA_float **mvn_dp, GMA_float **vxy_exy_qual)
{
    int32_t cntv,cntu,idx_vec_corr;
    int32_t id_cluster;
    for(cntv=0;cntv<dimy_vmap;cntv++)
    {
        for(cntu=0;cntu<dimy_vmap;cntu++)
        {
            idx_vec_corr=cntv*dimx_vmap+cntu;
            id_cluster=dpf->val[cntv][cntu];
            vxy_exy_qual[0]->val[cntv][cntu]=mvn_dp[idx_vec_corr]->val[id_cluster][0];
            vxy_exy_qual[1]->val[cntv][cntu]=mvn_dp[idx_vec_corr]->val[id_cluster][1];
            vxy_exy_qual[2]->val[cntv][cntu]=mvn_dp[idx_vec_corr]->val[id_cluster][2];
            vxy_exy_qual[3]->val[cntv][cntu]=mvn_dp[idx_vec_corr]->val[id_cluster][3];
            vxy_exy_qual[4]->val[cntv][cntu]=mvn_dp[idx_vec_corr]->val[id_cluster][4];
        }
    }
}

void GMA_float_conv2(GMA_float *in, GMA_float *kernel, GMA_float *out)
{
    int32_t cnt1, cnt2, cnt3,cnt4;
    int32_t dimx_in, dimy_in, dimx_k, dimy_k;
    int32_t ocwx, ocwy;
    float sum_dn,numel_k,dn_in,dn_min;
    const float N_A_N=sqrt(-1);
    
    dimx_in=(int32_t)in->ncols;
    dimy_in=(int32_t)in->nrows;
    dimx_k=(int32_t)kernel->ncols;
    dimy_k=(int32_t)kernel->nrows;

    //printf("Kernel:\n");
    //GMA_float_print(kernel);

    numel_k=(float)(dimy_k*dimx_k);

    ocwx=dimx_k/2;
    ocwy=dimy_k/2;

    //NaN filling
    for(cnt1=ocwy;cnt1<dimy_in-ocwy;cnt1++)
    {
         for(cnt2=ocwx;cnt2<dimx_in-ocwx;cnt2++)
        {
            sum_dn=0;
            for(cnt3=0;cnt3<dimy_k;cnt3++)
            {
                for(cnt4=0;cnt4<dimx_k;cnt4++)
                {
                    dn_in=(int32_t)(in->val[cnt1+cnt3-ocwy][cnt2+cnt4-ocwx]+0.5)?in->val[cnt1+cnt3-ocwy][cnt2+cnt4-ocwx]:N_A_N;
                    sum_dn+=dn_in*kernel->val[cnt3][cnt4];
                }
            }
            out->val[cnt1][cnt2]=sum_dn;
        }
    }

    //find the minimum DN
    dn_min=1e+37;
    for(cnt1=0;cnt1<dimy_in;cnt1++)
    {
         for(cnt2=0;cnt2<dimx_in;cnt2++)
        {
            if(out->val[cnt1][cnt2]<dn_min)
            {
                dn_min=out->val[cnt1][cnt2];
            }
        }
    }
    
    //shift the DN and substitute NaN values into 0
    for(cnt1=ocwy;cnt1<dimy_in-ocwy;cnt1++)
    {
         for(cnt2=ocwx;cnt2<dimx_in;cnt2++)
        {
            if(isnan(out->val[cnt1][cnt2]))
            {
                out->val[cnt1][cnt2]=0;
            }
            else
            {
                out->val[cnt1][cnt2]-=dn_min-1;
            }
        }
    }

}


int32_t GMA_float_find_median(GMA_float *in, int32_t column_of_interest)
{
    int32_t out;
    int32_t idx_order_median=(int32_t)((float)(in->nrows)/2+0.5);

    //copy the 
    //TODO: Keep implementing

}