//
// Created by arrouan on 24/11/16.
//

#ifndef RAEVOL_CUDA_ODE_SOLVER_ONLY_H_H
#define RAEVOL_CUDA_ODE_SOLVER_ONLY_H_H

//#include "ode_solver_only_omp_gpu.h"

void call_kernel_ode_cuda(int nb_gen,int multiply_population, int max_protein,
                          int nb_signal, int degradationstep,
                          int degradation_rate, double hill_shape_n, double hill_shape);


void call_kernel_ode_cuda_float(int nb_gen,int multiply_population, int max_protein,
                          int nb_signal, int degradationstep,
                          int degradation_rate, double hill_shape_n, double hill_shape);

void call_kernel_ode_cuda_thrust(int multiply_population,
                                 int nb_signal, int degradationstep,
                                 int degradation_rate, double hill_shape_n, double hill_shape);

double degradation_rate = 0.001;
double hill_shape_n      = 4;
double hill_shape_theta  = 0.5;
double hill_shape        = 0.0625; //std::pow( hill_shape_theta, hill_shape_n );

int degradationstep = 0;

void update_env_indiv(int lifestep, int indiv_id);
void solve_one_indiv_one_step(int indiv_id);

void update_env_list_indiv(int lifestep, int start_indiv_id, int end_indiv_id);
void solve_list_indiv_one_step(int start_indiv_id, int end_indiv_id);
void solve_one_indiv_one_step_tl(int indiv_id);
void solve_one_indiv_one_step_tl2(int indiv_id);
void solve_one_indiv_one_step_tl3(int lifestep, int multiply_population, int max_prot);

int transfert_data_to_gpu(int pop_size, int lifestep,
                          std::vector<std::vector<double>*> const &protein_concentration_list,
std::vector<std::vector<double>*> const &rna_basal_concentration_list,
std::vector<std::vector<std::vector<int>*>*> const &rna_produce_protein_list,
std::vector<std::vector<std::vector<double>*>*> const &rna_influence_enhancing_coef_list,
std::vector<std::vector<std::vector<double>*>*> const &rna_influence_operating_coef_list,
int nb_signal,
    std::vector<std::vector<double>*> const &env_concentration_list);

int transfert_data_to_gpu_float(int pop_size, int lifestep,
                          std::vector<std::vector<double>*> const &protein_concentration_list,
std::vector<std::vector<double>*> const &rna_basal_concentration_list,
std::vector<std::vector<std::vector<int>*>*> const &rna_produce_protein_list,
std::vector<std::vector<std::vector<double>*>*> const &rna_influence_enhancing_coef_list,
std::vector<std::vector<std::vector<double>*>*> const &rna_influence_operating_coef_list,
int nb_signal,
    std::vector<std::vector<double>*> const &env_concentration_list);

int transfert_data_to_gpu_dense(int pop_size, int lifestep,
                          std::vector<std::vector<double>*> const &protein_concentration_list,
std::vector<std::vector<double>*> const &rna_basal_concentration_list,
std::vector<std::vector<std::vector<int>*>*> const &rna_produce_protein_list,
std::vector<std::vector<std::vector<double>*>*> const &rna_influence_enhancing_coef_list,
std::vector<std::vector<std::vector<double>*>*> const &rna_influence_operating_coef_list,
int nb_signal,
    std::vector<std::vector<double>*> const &env_concentration_list);

int transfert_data_to_gpu_thrust(
    int nb_gen,
    int multiply_population,
    int nb_signal, int degradationstep,
    int degradation_rate, double hill_shape_n, double hill_shape, int pop_size, int lifestep,
                                 std::vector<std::vector<double>*> const &protein_concentration_list,
std::vector<std::vector<double>*> const &rna_basal_concentration_list,
std::vector<std::vector<std::vector<int>*>*> const &rna_produce_protein_list,
std::vector<std::vector<std::vector<double>*>*> const &rna_influence_enhancing_coef_list,
std::vector<std::vector<std::vector<double>*>*> const &rna_influence_operating_coef_list,
    std::vector<std::vector<double>*> const &env_concentration_list);

void transfer_to_tab(int pop_size, int lifestep);
void compute_openmp_gpu(int lifestep, int degradationstep, int nb_signal,
                        double hill_shape,double hill_shape_n,double degradation_rate );
void compute_openmp_gpu_less_tasks(int lifestep, int degradationstep, int nb_signal,
                        double hill_shape,double hill_shape_n,double degradation_rate );

std::vector<std::vector<double>*> protein_concentration_list;
std::vector<std::vector<double>*> rna_basal_concentration_list;
std::vector<std::vector<std::vector<int>*>*> rna_produce_protein_list;
std::vector<std::vector<std::vector<double>*>*> rna_influence_enhancing_coef_list;
std::vector<std::vector<std::vector<double>*>*> rna_influence_operating_coef_list;
int nb_signal;
std::vector<std::vector<double>*> env_concentration_list;

#endif //RAEVOL_CUDA_ODE_SOLVER_ONLY_H_H
