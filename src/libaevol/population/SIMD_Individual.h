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

    std::map<int32_t,promoterStruct*> promoters;
    std::map<int32_t,int32_t> leading_prom_pos;
    std::map<int32_t,int32_t> lagging_prom_pos;
    int count_prom = 0;

    std::set<int> terminator_lag;
    std::set<int> terminator_lead;
    std::vector<pRNA*> rnas;
    std::vector<pProtein*> proteins;
    double phenotype[300];
    double delta[300];
    double fitness;
    double metaerror;

    Dna_SIMD* dna_;

    int32_t indiv_id;
    int32_t parent_id;

    int32_t usage_count_ = 1;

    int32_t protein_count_ = 0;
    int32_t rna_count_ = 0;

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

    void remove_promoters_around(int32_t pos_1);
    void remove_promoters_around(int32_t pos_1, int32_t pos_2);
    void remove_all_promoters();

    void look_for_new_promoters_around(int32_t pos_1, int32_t pos_2);
    void look_for_new_promoters_around(int32_t pos);

    void locate_promoters();

    void move_all_promoters_after(int32_t pos, int32_t delta_pos);

    void duplicate_promoters_included_in(int32_t pos_1,
                                         int32_t pos_2,
                                         std::vector<std::list<promoterStruct*>>& duplicated_promoters);
    void extract_promoters_included_in(int32_t pos_1,
                                       int32_t pos_2, std::vector<std::list<promoterStruct*>>& extracted_promoters);
    void insert_promoters(std::vector<std::list<promoterStruct*>>& promoters_to_insert);
    void insert_promoters_at(std::vector<std::list<promoterStruct*>>& promoters_to_insert,
                                                   int32_t pos);

    void invert_promoters_included_in(int32_t pos1,
                                      int32_t pos2);


    static void shift_promoters(
        std::vector<std::list<promoterStruct*>>& promoters_to_shift,
        int32_t delta_pos,
        int32_t seq_length);
    static void invert_promoters(std::vector<std::list<promoterStruct*>>& promoter_lists,
                                 int32_t pos1,
                                 int32_t pos2);

    int8_t is_promoter_leading(int pos);
    int8_t is_promoter_lagging(int pos);
    void lst_promoters(bool lorl,
                                         Position before_after_btw, // with regard to the strand's reading direction
                                         int32_t pos1,
                                         int32_t pos2,
                                         std::list<promoterStruct*>& promoters_list);

 protected:
    void remove_leading_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2);
    void remove_leading_promoters_starting_after(int32_t pos);
    void remove_leading_promoters_starting_before(int32_t pos);

    void remove_lagging_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2);
    void remove_lagging_promoters_starting_after(int32_t pos);
    void remove_lagging_promoters_starting_before(int32_t pos);

    void move_all_leading_promoters_after(int32_t pos, int32_t delta_pos);
    void move_all_lagging_promoters_after(int32_t pos,int32_t delta_pos);

    void look_for_new_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2);
    void look_for_new_leading_promoters_starting_after(int32_t pos);
    void look_for_new_leading_promoters_starting_before(int32_t pos);

    void look_for_new_lagging_promoters_starting_between(int32_t pos_1,int32_t pos_2);
    void look_for_new_lagging_promoters_starting_after(int32_t pos);
    void look_for_new_lagging_promoters_starting_before(int32_t pos);

    void promoters_included_in(int32_t pos_1,
                                                     int32_t pos_2,
                                                     std::vector<std::list<promoterStruct*>>& promoters_list);

    void extract_leading_promoters_starting_between(int32_t pos_1,
                                                                          int32_t pos_2, std::list<promoterStruct*>& extracted_promoters);

    void extract_lagging_promoters_starting_between(int32_t pos_1,
                                                                          int32_t pos_2,
                                                                          std::list<promoterStruct*>& extracted_promoters);

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

    void clear_struct_before_next_step();
    void do_mutation();

    void run_a_step(double w_max, double selection_pressure,bool optim_prom = false);

    void start_stop_RNA();
    void opt_prom_compute_RNA();

    void compute_RNA();

    void start_protein();
    void compute_protein();
    void translate_protein(double w_max);
    void compute_phenotype();
    void compute_fitness(double selection_pressure);


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

    void selection();

    void selection(int indiv_id);

    void check_selection(int indiv_id);



};
}

#endif //RAEVOL_CUDA_SIMD_INDIVIDUAL_H