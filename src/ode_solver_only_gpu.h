//
// Created by arrouan on 05/12/16.
//

#ifndef RAEVOL_CUDA_ODE_SOLVER_ONLY_GPU_H
#define RAEVOL_CUDA_ODE_SOLVER_ONLY_GPU_H
#include <cuda.h>



#ifdef __CUDACC__
__global__
#endif
void process_delta(int nb_signal, int degradstep, int degradrate, int ***rna_produce_protein_array,
                   int **nb_rna_produce_protein, int *nb_rna_produce,   double **protein_concentration_array,
                   double **rna_basal_concentration_array, int *nb_protein_array, int *nb_rna_array,
                   double ***rna_influence_enhancing_coef_array, double ***rna_influence_operating_coef_array,
                   int **nb_rna_influence_enhancing_coef, int  **nb_rna_influence_operating_coef,
                   double **env_concentration_array, double hill_shape, double hill_shape_n);

#ifdef __CUDACC__
__global__
#endif
void process_delta_float(int nb_signal, int degradstep, int degradrate, int ***rna_produce_protein_array,
                   int **nb_rna_produce_protein, int *nb_rna_produce,   float **protein_concentration_array,
                   float **rna_basal_concentration_array, int *nb_protein_array, int *nb_rna_array,
                   float ***rna_influence_enhancing_coef_array, float ***rna_influence_operating_coef_array,
                   int **nb_rna_influence_enhancing_coef, int  **nb_rna_influence_operating_coef,
                   float **env_concentration_array, float hill_shape, float hill_shape_n);
#ifdef __CUDACC__
__global__
#endif
void process_delta_thrust(int nb_signal, int degradstep, int degradrate,
                          int max_prot, int max_rna,
                          int *nb_rna_produce_protein, int *rna_produce_protein,
                          double *protein_concentration_array,
                          double *rna_basal_concentration_array, int *nb_protein_array,
                          double *rna_influence_enhancing_coef_array, double *rna_influence_operating_coef_array,
                          int *nb_rna_influence_enhancing_coef,
                          double *env_concentration_array, double hill_shape, double hill_shape_n);

double **protein_concentration_array, **rna_basal_concentration_array;
float **f_protein_concentration_array, **f_rna_basal_concentration_array;
int *nb_protein_array, *nb_rna_array;

int ***rna_produce_protein_array;
int **nb_rna_produce_protein, *nb_rna_produce;

double ***rna_influence_enhancing_coef_array,
    ***rna_influence_operating_coef_array;
float ***f_rna_influence_enhancing_coef_array,
    ***f_rna_influence_operating_coef_array;
int **nb_rna_influence_enhancing_coef,
    **nb_rna_influence_operating_coef,
    *nb_rna_influence_enhancing_coef_l1,
    *nb_rna_influence_operating_coef_l1;

double **env_concentration_array;
float **f_env_concentration_array;

double **dev_protein_concentration_array,
    **dev_rna_basal_concentration_array;

float **f_dev_protein_concentration_array,
    **f_dev_rna_basal_concentration_array;
int *dev_nb_protein_array, *dev_nb_rna_array;

int ***dev_rna_produce_protein_array;
int **dev_nb_rna_produce_protein, *dev_nb_rna_produce;

double ***dev_rna_influence_enhancing_coef_array,
    ***dev_rna_influence_operating_coef_array;
float ***f_dev_rna_influence_enhancing_coef_array,
    ***f_dev_rna_influence_operating_coef_array;
int **dev_nb_rna_influence_enhancing_coef,
    **dev_nb_rna_influence_operating_coef,
    *dev_nb_rna_influence_enhancing_coef_l1,
    *dev_nb_rna_influence_operating_coef_l1;

double **dev_env_concentration_array;
float **f_dev_env_concentration_array;

/* dense tab CPU */
double *protein_tab;
double *rna_basal_tab;
int *rna_produce_protein_tab;
double *rna_influence_enhancing_tab;
double *rna_influence_operating_tab;
double *env_concentration_tab;

/* dense tab GPU */

double *dev_protein_tab;
double *dev_rna_basal_tab;
int *dev_rna_produce_protein_tab;
double *dev_rna_influence_enhancing_tab;
double *dev_rna_influence_operating_tab;
double *dev_env_concentration_tab;

size_t pitch_protein_tab;
size_t pitch_rna_basal_tab;
size_t pitch_rna_produce_protein_tab;
size_t pitch_rna_influence_enhancing_tab;
size_t pitch_rna_influence_operating_tab;
size_t pitch_env_concentration_tab;

/* thrust vector */


#endif //RAEVOL_CUDA_ODE_SOLVER_ONLY_GPU_H
