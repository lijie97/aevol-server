//
// Created by arrouan on 27/07/17.
//

#ifndef RAEVOL_CUDA_SIMD_INDIVIDUAL_H
#define RAEVOL_CUDA_SIMD_INDIVIDUAL_H

#include "ExpManager.h"

namespace aevol {

constexpr const char* PROM_SEQ_LEAD = "0101011001110010010110";
constexpr const char* PROM_SEQ_LAG  = "1010100110001101101001";

constexpr const char* SHINE_DAL_SEQ_LEAD = "011011000";
constexpr const char* SHINE_DAL_SEQ_LAG  = "100100111";

constexpr const char* PROTEIN_END_LEAD  = "001";
constexpr const char* PROTEIN_END_LAG   = "110";


class promoterStruct {
 public:
    promoterStruct(int32_t t_pos, int8_t t_error, bool lead) {
      pos = t_pos; error = t_error; leading_or_lagging = lead;
    }

    int32_t pos = -1;
    int8_t error = -1;
    bool leading_or_lagging; // TRUE = leading / FALSE = lagging
};

class pRNA {
 public:
    pRNA(int32_t t_begin, int32_t t_end, int8_t t_leading_lagging, double t_e,
         int32_t t_length) {
      begin = t_begin;
      end = t_end;
      leading_lagging = t_leading_lagging;
      e = t_e;
      length = t_length;
    }

    int32_t begin;
    int32_t end;
    int8_t leading_lagging; // 0 = leading, 1 = lagging
    double e;
    std::vector<int32_t> start_prot;
    int32_t length;
};

class pProtein {
 public:
    pProtein(int32_t t_protein_start,
      int32_t t_protein_end,
      int32_t t_protein_length,
      int8_t t_leading_lagging, double t_e) {
      protein_start = t_protein_start;
      protein_end = t_protein_end;
      protein_length = t_protein_length;
      leading_lagging = t_leading_lagging;
      e = t_e;
    }

    int32_t protein_start;
    int32_t protein_end;
    int32_t protein_length;
    int8_t leading_lagging; // 0 = leading, 1 = lagging
    float m;
    float w;
    float h;
    double e;
    bool is_functional;

};

class Internal_SIMD_Struct {
 public:
    Internal_SIMD_Struct() = default;

    std::vector<promoterStruct*> promoters;
    std::set<int> terminator_lag;
    std::set<int> terminator_lead;
    std::vector<pRNA> rnas;
    std::vector<pProtein> proteins;
    float phenotype[300];
    float delta[300];
    double fitness;
    float metaerror;
};

class SIMD_Individual {
 public:
    explicit SIMD_Individual(ExpManager* exp_m);
    ~SIMD_Individual();

    void run_a_step(double w_max, double selection_pressure);

    void start_stop_RNA();
    void compute_RNA();
    void start_protein();
    void compute_protein();
    void translate_protein(double w_max);
    void compute_phenotype();
    void compute_fitness(double selection_pressure);

    void check_result();

 private:
    ExpManager* exp_m_;
    Internal_SIMD_Struct** internal_simd_struct;
    int* dna_size;
    float target[300];
};
}

#endif //RAEVOL_CUDA_SIMD_INDIVIDUAL_H
