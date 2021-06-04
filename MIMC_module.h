#ifndef _GMA_
#include "GMA.h"
#endif

#ifndef _STDINT_H_
#include <stdint.h>
#endif
#define _MIMC2_MODULE_

typedef struct param{
    int32_t vec_ocw[4];
    //uint32_t min_num_cp;
    float AW_CRE;
    float AW_SF;
    float spacing_grid;
    float radius_neighbor;
    float radius_neighbor_dpf1;
    float radius_neighbor_ps;
    float meter_per_spacing;
    float mpp;
    int32_t num_cp_max;
    int32_t num_cp_min;
    float ratio_cp;
    float thres_spd_cp;
} param;

//external global varaibles (originally declared in MIMC2_main.c)
extern float dt;   //temporal baseline
extern int32_t num_dp; //number of multiple matching attempts
extern int32_t num_grid,dimx_vmap,dimy_vmap; //# of grids, dimension in mapx and mapy direction
extern param param_mimc2;  //parameters that controlls the software
extern GMA_float **kernel;
//control point determination module
int get_offset_image(GMA_float *i0,GMA_float *i1,GMA_float **kernel, GMA_double *xyuvav, int32_t *offset, GMA_uint8 *flag_cp);
void GMA_double_randperm_row(GMA_double *var);


// feature tracking modules
GMA_int32** get_uv_pivot(GMA_double *xyuvav, float dt, param param_mimc2, int32_t ocw, GMA_float *i1);
char investigate_valid_grid(GMA_float *refchip,GMA_float *sarea);
void find_ncc_peak(GMA_float *refchip, GMA_float *sarea, GMA_int32 *uv_pivot, float *uvncc);
void extract_refchip(GMA_float *i0, int32_t *uv0, int32_t ocw,GMA_float *refchip);
GMA_float* extract_sarea(GMA_float *i1, int32_t *uv0, int32_t ocw, GMA_int32 *uv_pivot);
GMA_float* matching_ncc_dlc_2(GMA_float *i0,GMA_float *i1,GMA_double *xyuvav,int32_t *offset,GMA_int32 **uv_pivot,int32_t ocw,float AW_CRE,float AW_SF);


//routines for postprocessing
GMA_float** mimc2_postprocess(GMA_float **dp, GMA_double *xyuvav, float dt);    //Root function for the postprocessing
GMA_float** calc_mean_var_num_dp_cluster(GMA_float** dp,int32_t num_dpoi);
void cluster_euclidian(GMA_float *dp_stack, float min_dist, GMA_uint8 *id_cluster);
void mark_row(int32_t nrow, uint8_t id_dp, GMA_uint8 *id_cluster, GMA_uint8 *mtrx_dist);   //resursive sub-funciton for "cluster_euclidian"
GMA_uint8* get_mask_nomatching(GMA_float **mvn_dp);
GMA_int32* get_dpf0(GMA_float **mvn_dp, float min_matching_ratio);  //get the most probable dp
GMA_int32* get_ruv_neighbor(GMA_double *xyuvav,float radius_neighbor);

void get_dpf1(GMA_int32 *dpf0, GMA_float *dpf_dx, GMA_float *dpf_dy, GMA_int32 *ruv_neighbor, GMA_float **mvn_dp, GMA_double *xyuvav);  //get the initial dp for pseudosmoothing
void get_dpf1_original(GMA_int32 *dpf0, GMA_float *dpf_dx, GMA_float *dpf_dy, GMA_int32 *ruv_neighbor, GMA_float **mvn_dp, GMA_double *xyuvav);  //get the initial dp for pseudosmoothing
void get_dpf_pseudosmoothing(GMA_int32 *dpf,GMA_float *dpf_dx, GMA_float *dpf_dy,GMA_int32 *ruv_neighbor,GMA_float **mvn_dp, GMA_double *xyuvav);  //pseudosmoothing routine

void quadfit2(GMA_int32 *xy, GMA_double *z, GMA_double *w, GMA_double *xyi, GMA_double *out);
void GMA_double_mul(GMA_double *a, GMA_double *b, GMA_double *out);
void GMA_double_inv(GMA_double *a, GMA_double *inva);

void convert_dpf_to_vxy_exy_qual(GMA_int32 *dpf, GMA_float **mvn_dp, GMA_float **vxy_exy_qual);

//2D convolution function
void GMA_float_conv2(GMA_float *in, GMA_float *kernel, GMA_float *out);

//find median from the column vector
int32_t GMA_float_find_median(GMA_float *in, int32_t column_of_interest);
