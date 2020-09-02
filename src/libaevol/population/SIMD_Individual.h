//
// Created by arrouan on 27/07/17.
//

#ifndef RAEVOL_CUDA_SIMD_INDIVIDUAL_H
#define RAEVOL_CUDA_SIMD_INDIVIDUAL_H

#include <map>
#include <raevol/SIMD_PhenotypicTargetHandler_R.h>
#include <set>
#include <vector>

//#include "ExpManager.h"
//#include "Stats_SIMD.h"
//#include "Dna_SIMD.h"
#include "Observable.h"
#include "PhenotypicTargetHandler.h"
#include "ae_enums.h"
#include "Stats.h"
#include "SIMD_DnaFactory.h"
#include "Vector_Fuzzy.h"
#include "raevol/Protein_R.h"

namespace aevol {

    class ExpManager;
    class Dna_SIMD;
    class Stats_SIMD;
    class SIMD_Abstract_Metadata;
    class SIMD_Map_Metadata;
    class SIMD_PhenotypicTargetHandler_R;

constexpr const char* PROM_SEQ_LEAD = "0101011001110010010110";
constexpr const char* PROM_SEQ_LAG  = "1010100110001101101001";

constexpr const char* SHINE_DAL_SEQ_LEAD = "011011000";
constexpr const char* SHINE_DAL_SEQ_LAG  = "100100111";

constexpr const char* PROTEIN_END_LEAD  = "001";
constexpr const char* PROTEIN_END_LAG   = "110";

const int32_t PHENOTYPE_VECTOR_SIZE = 120000;
    constexpr const  double D_PHENOTYPE_VECTOR_SIZE = 120000.0;

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

class pProtein;

class AffinityFactor {
 public:
  AffinityFactor(pProtein* prot, double efactor, double ofactor) {
    protein = prot;
    enhancer_factor = efactor;
    operator_factor = ofactor;
  }

  pProtein* protein;
  double enhancer_factor;
  double operator_factor;

  double concentration();
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
      start_prot_count_ = 0;
    }

    ~pRNA() {
    }

#ifdef __REGUL
    std::list<AffinityFactor> affinity_list;

    int nb_influences_ = 0;

    int32_t enhancer_position(int32_t length) {
      if(leading_lagging == LEADING)
      {
        return (begin - 20)  % ( length ) < 0 ?
               ((begin - 20)  % ( length )) + ( length ) :
               (begin - 20)  % ( length );
      }
      else  // strand_ = LAGGING
      {
        return (begin + 20)  % ( length ) < 0 ?
               ((begin + 20)  % ( length )) + ( length ) :
               (begin + 20)  % ( length );
      }
    }

  int32_t operator_position(int32_t length) {
    if(leading_lagging == LEADING)
    {
      return (begin + PROM_SIZE)  % ( length ) < 0 ?
             (begin + PROM_SIZE)  % ( length ) + (length) :
             (begin + PROM_SIZE)  % ( length );
    }
    else  // strand_ = LAGGING
    {
      return (begin - PROM_SIZE)  % ( length ) < 0 ?
             (begin - PROM_SIZE)  % ( length ) + (length) :
             (begin - PROM_SIZE)  % ( length );
    }
  }

  double affinity_with_protein( int32_t index, pProtein *protein,
                                     Internal_SIMD_Struct* indiv,
                               ExpManager* exp_m);

  double compute_synthesis_rate(Internal_SIMD_Struct* indiv);
#endif

    int32_t begin;
    int32_t end;
    int8_t leading_lagging; // 0 = leading, 1 = lagging
    double e;
    std::list<int32_t > start_prot;
    int32_t start_prot_count_;
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
      int8_t t_leading_lagging, double t_e, pRNA* rna) {
      protein_start = t_protein_start;
      protein_end = t_protein_end;
      protein_length = t_protein_length;
      leading_lagging = t_leading_lagging;
      e = t_e;

      is_init_ = true;

#ifdef __REGUL
      rna_list_.push_back(rna);
      initial_e_ = e;
#endif
    }

#ifdef __REGUL
    pProtein(Protein_R* prot_sig);


  pProtein(pProtein* prot);
#endif

    bool operator<(const pProtein & other);

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

    int8_t codon_list[64] = {};
    int16_t nb_codons_ = 0;

#ifdef __REGUL
  bool is_TF_;

  double initial_e_ = -1;
  double    delta_concentration_ = 0;
  bool      inherited_ = false;
  bool      signal_ = false;
  std::list<pRNA*> rna_list_;
#endif
};

class Internal_SIMD_Struct : public Observable {
 public:
    Internal_SIMD_Struct(ExpManager* exp_m, double w_max, SIMD_DnaFactory* dna_factory);

    Internal_SIMD_Struct(ExpManager* exp_m, Internal_SIMD_Struct* clone, SIMD_DnaFactory* dna_factory);

    ~Internal_SIMD_Struct();

#ifdef PHENOTYPE_VECTOR
    double phenotype[PHENOTYPE_VECTOR_SIZE];
    double delta[PHENOTYPE_VECTOR_SIZE];
#else
        Vector_Fuzzy* phenotype = nullptr;
        Vector_Fuzzy* delta = nullptr;
#endif
    double fitness;
    double metaerror;


  double* fitness_by_env_id_;
  double* metaerror_by_env_id_;

    Dna_SIMD* dna_;

    int32_t indiv_id;
    int32_t parent_id;

    int32_t usage_count_ = 1;

    ExpManager* exp_m_;
    SIMD_DnaFactory* dna_factory_;

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
    void compute_fitness(int indiv_id, double selection_pressure, int env_id = -1);


    void check_result();
    void check_dna();
    void check_struct();
    void check_individual(int indiv_id, int x, int y);
    bool standalone() const { return standalone_; }


    void build_phenotypic_target(PhenotypicTargetHandler* phenotypic_target_handler);

#ifdef __REGUL
  void compute_network(int indiv_id, double selection_pressure);
  void update_network(int indiv_id, double selection_pressure);
  void evaluate_network(int indiv_id, double selection_pressure, int env_id);
  void finalize_network(int indiv_id, double selection_pressure);
  void solve_network(int indiv_id, double selection_pressure);
  void update_phenotype( int indiv_id );


  SIMD_PhenotypicTargetHandler_R* phenotypic_target_handler_;
#endif

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

    SIMD_DnaFactory* dna_factory_;

    double* fitness_sum_tab_;

 private:
    ExpManager* exp_m_;
    int* dna_size;
#ifdef __REGUL
Vector_Fuzzy** targets;
#else
#ifdef PHENOTYPE_VECTOR
    double* target;
#else
    Vector_Fuzzy* target;
#endif
#endif
    bool standalone_;
    bool first_gener_ = true;

    Stats_SIMD* stats_best = nullptr;
    Stats_SIMD* stats_mean = nullptr;



    long apply_mutation[1024];


    Stats* stats_;

    void selection(int indiv_id);

    void check_selection(int indiv_id);



};
}

#endif //RAEVOL_CUDA_SIMD_INDIVIDUAL_H
