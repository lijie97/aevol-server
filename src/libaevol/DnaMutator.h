//
// Created by arrouan on 11/01/18.
//

#ifndef RAEVOL_CUDA_DNAMUTATOR_H
#define RAEVOL_CUDA_DNAMUTATOR_H

#include "Individual.h"
#include "JumpingMT.h"
#include "MutationEvent.h"

namespace aevol {
class DnaMutator {
 public:
    DnaMutator(std::shared_ptr<JumpingMT> mut_prng, int32_t length,
               double duplication_rate,
               double deletion_rate, double translocation_rate,
               double inversion_rate,
               double point_mutation_rate, double small_insertion_rate,
               double small_deletion_rate, int16_t max_indel_size,
               int32_t min_genome_length, int32_t max_genome_length, int indiv_id, int x, int y);

    DnaMutator(Individual* indiv, int x, int y);

    ~DnaMutator() {
      int cpt = 0;
      for (auto repl : mutation_list_) {
        delete repl;
        cpt++;
      }

      mutation_list_.clear();
    }

    void generate_mutations();
    void generate_rearrangements();
    void generate_small_mutations();

    MutationEvent* generate_next_mutation(int32_t length);

    bool mutation_available() { return ((cpt_rear_ + cpt_mut_) > 0); }

    std::list<MutationEvent*> mutation_list_;

    bool hasMutate() {return hasMutate_;}

    bool setMutate(bool mutate) {hasMutate_ = mutate;}

    int x_,y_;

 private:
    std::shared_ptr<JumpingMT> mut_prng_;
    int32_t length_;

    // ------------------------------ Rearrangement rates (without alignements)
    double duplication_rate_;
    double deletion_rate_;
    double translocation_rate_;
    double inversion_rate_;

    // --------------------------------------------------------- Mutation rates
    double point_mutation_rate_;
    double small_insertion_rate_;
    double small_deletion_rate_;
    int16_t max_indel_size_;

    //--------------------------- Mutation counters
    int32_t nb_swi_;
    int32_t nb_ins_;
    int32_t nb_del_;
    int32_t nb_mut_;

    int32_t nb_large_dupl_;
    int32_t nb_large_del_;
    int32_t nb_large_trans_;
    int32_t nb_large_inv_;
    int32_t nb_rear_;

    int32_t cpt_rear_;
    int32_t cpt_mut_;

    int32_t min_genome_length_;
    int32_t max_genome_length_;

    bool hasMutate_ = false;
    unsigned long long int id_;
};
}


#endif //RAEVOL_CUDA_DNAMUTATOR_H
