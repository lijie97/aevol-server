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
               double small_deletion_rate, int16_t max_indel_size);

    DnaMutator(Individual* indiv);

    void generate_mutations();
    void generate_rearrangements();
    void generate_small_mutations();

    std::list<MutationEvent*> mutation_list_;

    bool hasMutate() {return hasMutate_;}

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

    bool hasMutate_ = false;
};
}


#endif //RAEVOL_CUDA_DNAMUTATOR_H
