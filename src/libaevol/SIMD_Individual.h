//
// Created by arrouan on 27/07/17.
//

#ifndef RAEVOL_CUDA_SIMD_INDIVIDUAL_H
#define RAEVOL_CUDA_SIMD_INDIVIDUAL_H

#include <vector>
#include <map>
#include <set>

//#include "ExpManager.h"
//#include "Stats_SIMD.h"
//#include "Dna_SIMD.h"
#include "Observable.h"
#include "PhenotypicTargetHandler.h"
#include "ae_enums.h"
#include "Stats.h"

namespace aevol {

    class ExpManager;
    class Dna_SIMD;
    class Stats_SIMD;
    class SIMD_Abstract_Metadata;
    class SIMD_Map_Metadata;

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

    promoterStruct(const promoterStruct& clone) {
      pos=clone.pos;error=clone.error;leading_or_lagging=clone.leading_or_lagging;
    }

    promoterStruct(promoterStruct* clone) {
      pos=clone->pos;error=clone->error;leading_or_lagging=clone->leading_or_lagging;
    }

    int32_t pos = -1;
    int8_t error = -1;
    bool leading_or_lagging; // TRUE = leading / FALSE = lagging
};

class pRNA {
 public:
    pRNA() {};
    pRNA(int32_t t_begin, int32_t t_end, int8_t t_leading_lagging, double t_e,
         int32_t t_length) {
      begin = t_begin;
      end = t_end;
      leading_lagging = t_leading_lagging;
      e = t_e;
      length = t_length;
      is_coding_ = false;
      is_init_ = true;
    }

    int32_t begin;
    int32_t end;
    int8_t leading_lagging; // 0 = leading, 1 = lagging
    double e;
    std::vector<int32_t> start_prot;
    int32_t length;
    bool is_coding_;

    bool is_init_ = false;
};

class pProtein {
 public:
    pProtein() {};
    pProtein(int32_t t_protein_start,
      int32_t t_protein_end,
      int32_t t_protein_length,
      int8_t t_leading_lagging, double t_e) {
      protein_start = t_protein_start;
      protein_end = t_protein_end;
      protein_length = t_protein_length;
      leading_lagging = t_leading_lagging;
      e = t_e;
      is_init_ = true;
    }

    int32_t protein_start;
    int32_t protein_end;
    int32_t protein_length;
    int8_t leading_lagging; // 0 = leading, 1 = lagging
    double m;
    double w;
    double h;
    double e;
    bool is_functional;

    bool is_init_ = false;
};

class Internal_SIMD_Struct : public Observable {
 public:
    Internal_SIMD_Struct(ExpManager* exp_m, double w_max);

    Internal_SIMD_Struct(ExpManager* exp_m, Internal_SIMD_Struct* clone, bool copy_dna = true);

    ~Internal_SIMD_Struct();

    double phenotype[300];
    double delta[300];
    double fitness;
    double metaerror;

    Dna_SIMD* dna_;

    int32_t indiv_id;
    int32_t parent_id;

    int32_t usage_count_ = 1;

    ExpManager* exp_m_;

    int global_id = -1;

    double w_max_;

    /** Variables for Tree mgmt **/
    int32_t  nb_genes_activ = 0;
    int32_t  nb_genes_inhib = 0;
    int32_t  nb_func_genes = 0;
    int32_t  nb_non_func_genes = 0;
    int32_t  nb_coding_RNAs = 0;
    int32_t  nb_non_coding_RNAs = 0;
    /** END of Variables for Tree Mgmt **/


    void reset_stats() {
        nb_genes_activ = 0;
        nb_genes_inhib = 0;
        nb_func_genes = 0;
        nb_non_func_genes = 0;
        nb_coding_RNAs = 0;
        nb_non_coding_RNAs = 0;
    }

    void rebuild_index();

    SIMD_Abstract_Metadata* metadata_;
};

class PromoterList {
    PromoterList() = default;

    std::map<int32_t,int> leading_prom_pos;
    std::map<int32_t,int> lagging_prom_pos;
};

class SIMD_Individual : public Observable{
 public:
    SIMD_Individual(ExpManager* exp_m);

    ~SIMD_Individual();

    void run_a_step(double w_max, double selection_pressure,bool optim_prom = false);

    void do_mutation(int indiv_id);

    void start_stop_RNA(int indiv_id);
    void opt_prom_compute_RNA(int indiv_id);
    void compute_RNA(int indiv_id);
    void start_protein(int indiv_id);
    void compute_protein(int indiv_id);
    void translate_protein(int indiv_id, double w_max);
    void compute_phenotype(int indiv_id);
    void compute_fitness(int indiv_id, double selection_pressure);


    void check_result();
    void check_dna();
    void check_struct();
    void check_individual(int indiv_id, int x, int y);
    bool standalone() const { return standalone_; }


    void build_phenotypic_target(PhenotypicTargetHandler* phenotypic_target_handler);


    void set_stats(Stats* stats) { stats_ = stats; }

    Internal_SIMD_Struct** internal_simd_struct;
    Internal_SIMD_Struct** prev_internal_simd_struct;
    Internal_SIMD_Struct* best_indiv;

    int32_t* next_generation_reproducer_;

    int32_t nb_indivs_;
    int32_t nb_clones_;

    static bool standalone_simd;//= true;
    int rna_grain_size = 32;
    int protein_grain_size = 32;

 private:
    ExpManager* exp_m_;
    int* dna_size;
    double* target;
    bool standalone_;

    Stats_SIMD* stats_best = nullptr;
    Stats_SIMD* stats_mean = nullptr;


    Stats* stats_;

    void selection(int indiv_id);

    void check_selection(int indiv_id);



};
}

#endif //RAEVOL_CUDA_SIMD_INDIVIDUAL_H
