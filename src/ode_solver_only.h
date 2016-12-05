//
// Created by arrouan on 24/11/16.
//

#ifndef RAEVOL_CUDA_ODE_SOLVER_ONLY_H_H
#define RAEVOL_CUDA_ODE_SOLVER_ONLY_H_H

#include "ode_solver_only_gpu.h"

double degradation_rate = 1;
double hill_shape_n      = 4;
double hill_shape_theta  = 0.5;
double hill_shape        = std::pow( hill_shape_theta, hill_shape_n );

int degradationstep = 0;

void update_env_indiv(int lifestep, int indiv_id);
void solve_one_indiv_one_step(int indiv_id);

void update_env_list_indiv(int lifestep, int start_indiv_id, int end_indiv_id);
void solve_list_indiv_one_step(int start_indiv_id, int end_indiv_id);

std::vector<std::vector<double>*> protein_concentration_list;
std::vector<std::vector<double>*> rna_basal_concentration_list;
std::vector<std::vector<std::vector<int>*>*> rna_produce_protein_list;
std::vector<std::vector<std::vector<double>*>*> rna_influence_enhancing_coef_list;
std::vector<std::vector<std::vector<double>*>*> rna_influence_operating_coef_list;
int nb_signal;
std::vector<std::vector<double>*> env_concentration_list;


#endif //RAEVOL_CUDA_ODE_SOLVER_ONLY_H_H
