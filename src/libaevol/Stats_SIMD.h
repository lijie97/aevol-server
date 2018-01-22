//
// Created by arrouan on 19/01/18.
//

#ifndef RAEVOL_CUDA_STATS_SIMD_H
#define RAEVOL_CUDA_STATS_SIMD_H


#include <cstdint>
#include <fstream>

#include "SIMD_Individual.h"

namespace aevol {

class Stats_SIMD {
 public:
    Stats_SIMD(SIMD_Individual* simd_individual, int64_t generation,
               bool best_or_not);

    ~Stats_SIMD() {
      if (is_indiv_) {
        statfile_best_.flush();
        statfile_best_.close();
      } else {
        statfile_mean_.flush();
        statfile_mean_.close();
      }
    }

    void compute_best();
    void compute_average();

    void write_best();
    void write_average();

    bool is_indiv() { return is_indiv_; }


 protected:
    SIMD_Individual* simd_individual_;

    int64_t generation_;

    bool is_indiv_;

    int32_t pop_size_;

    double fitness_ = 0;
    double metabolic_error_ = 0;

    int32_t amount_of_dna_ = 0;
    int32_t nb_coding_rnas_ = 0;
    int32_t nb_non_coding_rnas_ = 0;

    int32_t nb_functional_genes_ = 0;
    int32_t nb_non_functional_genes_ = 0;

    int32_t nb_mut_ = 0;
    int32_t nb_rear_ = 0;
    int32_t nb_switch_ = 0;
    int32_t nb_indels_ = 0;
    int32_t nb_dupl_ = 0;
    int32_t nb_del_ = 0;
    int32_t nb_trans_ = 0;
    int32_t nb_inv_ = 0;

    double dupl_rate_ = 0;
    double del_rate_ = 0;
    double trans_rate_ = 0;
    double inv_rate_ = 0;

    int32_t nb_bases_in_0_CDS_ = 0;
    int32_t nb_bases_in_0_functional_CDS_ = 0;
    int32_t nb_bases_in_0_non_functional_CDS_ = 0;
    int32_t nb_bases_in_0_RNA_ = 0;
    int32_t nb_bases_in_0_coding_RNA_ = 0;
    int32_t nb_bases_in_0_non_coding_RNA_ = 0;

    int32_t nb_bases_non_essential_ = 0;
    int32_t nb_bases_non_essential_including_nf_genes_ = 0;

    bool is_computed_ = false;

    // Stats

    std::ofstream statfile_best_;
    std::ofstream statfile_mean_;
};

}
#endif //RAEVOL_CUDA_STATS_SIMD_H
