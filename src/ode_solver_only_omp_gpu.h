//
// Created by arrouan on 13/12/16.
//

#ifndef RAEVOL_CUDA_ODE_SOLVER_ONLY_OMP_GPU_H
#define RAEVOL_CUDA_ODE_SOLVER_ONLY_OMP_GPU_H


//#define LOOKUP_TABLE_SIZE 10000000

int max_prot = 0, max_rna = 0;

double *omp_protein_tab;
double *omp_rna_basal_tab;
int *omp_rna_produce_protein_tab;
int *nb_protein;
double *omp_rna_influence_enhancing_tab;
double *omp_rna_influence_operating_tab;
double *omp_env_concentration_tab;
static double lookup_table_pow[LOOKUP_TABLE_SIZE];



#endif //RAEVOL_CUDA_ODE_SOLVER_ONLY_OMP_GPU_H
