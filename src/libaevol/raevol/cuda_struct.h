//
// Created by arrouan on 10/02/16.
//


#ifndef AEVOL_CUDA_CUDA_STRUCT_H
#define AEVOL_CUDA_CUDA_STRUCT_H

#include "ExpManager.h"
#include "HybridFuzzy.h"
#include "Individual_R.h"

namespace aevol {
class cuda_struct {
    friend class Individual_R;
    friend class Individual;
 public:
    cuda_struct() {};

    ~cuda_struct() {};

    void init_struct(int max_protein, int max_rna, int max_influence,
                     int nb_signals, int life_time, int nb_eval,
                     double selection_pressure);

    void transfert_to_gpu(ExpManager* exp_m);

    void compute_a_generation(ExpManager* exp_m);

    void get_data_to_cpu();

    void print_dist(ExpManager* exp_m);

    void delete_struct();

 private:
    double* phenotype_inhib;
    double* phenotype_activ;
    double* phenotype;

    double* environment;

    int* eval_step;

    double* protein_concentration;

    double* signals_concentration;

    double* enhance_coef;
    double* operate_coef;
    double* rna_synthesis;
    double* basal_level;

    int* protein_influence; // index in protein concentration

    double* protein_influenced;

    int* protein_triangle_ix0;
    int* protein_triangle_ix1;
    int* protein_triangle_ix2;

    double* protein_triangle_height;

    double* delta;
    double *dist_sum;

    // Max
    int max_protein_;
    int max_rna_;
    int max_influence_;
    int nb_signals_;
    int life_time_;
    int nb_eval_;

    double selection_pressure_;
};
}
#endif //AEVOL_CUDA_CUDA_STRUCT_H
