//
// Created by Jonathan Rouzaud-Cornabas on 24/08/17
//

#ifndef RAEVOL_CUDA_DNA_SIMD_H
#define RAEVOL_CUDA_DNA_SIMD_H


#include <cstdint>

#include "Dna.h"
#include "MutationEvent.h"


namespace aevol {
class Dna_SIMD {
 public:
    Dna_SIMD(Dna* dna, Internal_SIMD_Struct* indiv);
    Dna_SIMD(Dna* dna);
    ~Dna_SIMD();

    void apply_mutations();

    bool do_switch(int32_t pos);
    bool do_small_insertion(int32_t pos,
                                      int16_t nb_insert, char* seq);
    bool do_small_deletion(int32_t pos, int16_t nb_del);
    bool do_duplication(int32_t pos_1, int32_t pos_2, int32_t pos_3);
    bool do_inversion(int32_t pos_1, int32_t pos_2);
    bool do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3,
                               int32_t pos_4, bool invert);
    bool do_deletion(int32_t pos_1, int32_t pos_2);

    void remove(int32_t first, int32_t last);
    void insert(int32_t pos, const char* seq, int32_t seq_length = -1);
    void replace(int32_t pos, char* seq, int32_t seq_length = -1);

    const char* data() const {return data_;}
    int32_t length() const {return length_;}

    static inline int32_t nb_blocks(int32_t length);

    void ABCDE_to_ADCBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                             int32_t pos_E);
    void ABCDE_to_ADBpCpE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                               int32_t pos_E);
    void ABCDE_to_ACpDpBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                               int32_t pos_E);

    std::list<MutationEvent*> mutation_list;
    char* data_;
 private:

    int32_t length_;
    int32_t nb_blocks_;
    Internal_SIMD_Struct* indiv_;

};

int32_t Dna_SIMD::nb_blocks(int32_t length) {
  return length/BLOCK_SIZE + 1;
}
}


#endif //RAEVOL_CUDA_DNA_SIMD_H
