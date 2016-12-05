//
// Created by arrouan on 05/12/16.
//

#ifndef RAEVOL_CUDA_ODE_SOLVER_ONLY_GPU_H
#define RAEVOL_CUDA_ODE_SOLVER_ONLY_GPU_H

double degradation_rate = 1;
double hill_shape_n      = 4;
double hill_shape_theta  = 0.5;
double hill_shape        = std::pow( hill_shape_theta, hill_shape_n );

int degradationstep = 0;

double **protein_concentration_array, **rna_basal_concentration_array;
int *nb_protein_array, *nb_rna_array;

int ***rna_produce_protein_array;
int **nb_rna_produce_protein, *nb_rna_produce;

double ***rna_influence_enhancing_coef_array,
    ***rna_influence_operating_coef_array;
int **nb_rna_influence_enhancing_coef,
    **nb_rna_influence_operating_coef;

double **env_concentration_array;

double **dev_protein_concentration_array,
    **dev_rna_basal_concentration_array;
int *dev_nb_protein_array, *dev_nb_rna_array;

int ***dev_rna_produce_protein_array;
int **dev_nb_rna_produce_protein, *dev_nb_rna_produce;

double ***dev_rna_influence_enhancing_coef_array,
    ***dev_rna_influence_operating_coef_array;
int **dev_nb_rna_influence_enhancing_coef,
    **dev_nb_rna_influence_operating_coef,
    *dev_nb_rna_influence_enhancing_coef_l1,
    *dev_nb_rna_influence_operating_coef_l1;

double **dev_env_concentration_array;

#endif //RAEVOL_CUDA_ODE_SOLVER_ONLY_GPU_H
