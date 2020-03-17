//
// Created by Jonathan Rouzaud-Cornabas on 24/08/17
//

#ifndef RAEVOL_CUDA_DNA_SIMD_H
#define RAEVOL_CUDA_DNA_SIMD_H


#include <cstdint>
#include <list>

//#include "Dna.h"
#include "MutationEvent.h"
#include "BitSet_SIMD.h"

namespace aevol {

#define BLOCK_SIZE INT32_C(1024)
class Dna;
class Internal_SIMD_Struct;

class Dna_SIMD {
 public:
    Dna_SIMD(Dna* dna, Internal_SIMD_Struct* indiv);
    Dna_SIMD(Dna_SIMD* dna, Internal_SIMD_Struct* indiv, bool copy_dna = true);
    Dna_SIMD(Dna* dna);
    ~Dna_SIMD();

    inline char get_lead(int32_t pos) { return data_[pos - ((unsigned int32_t)(pos - length_) >> 31) * length_];};
    inline char get_lag(int32_t pos) { return data_[pos + ((unsigned int32_t)(pos) >> 31) * length_];};

    void apply_mutations();

    void apply_mutations_standalone();

    bool do_switch(int32_t pos);
#ifdef WITH_BITSET
    bool do_small_insertion(int32_t pos, BitSet_SIMD* seq);
#else
    bool do_small_insertion(int32_t pos,
                                      int16_t nb_insert, char* seq);
#endif
    bool do_small_deletion(int32_t pos, int16_t nb_del);
    bool do_duplication(int32_t pos_1, int32_t pos_2, int32_t pos_3);
    bool do_inversion(int32_t pos_1, int32_t pos_2);
    bool do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3,
                               int32_t pos_4, bool invert);
    bool do_deletion(int32_t pos_1, int32_t pos_2);

#ifndef WITH_BITSET
    void remove(int32_t first, int32_t last);
    void insert(int32_t pos, const char* seq, int32_t seq_length = -1);
    void replace(int32_t pos, char* seq, int32_t seq_length = -1);
#endif

#ifndef WITH_BITSET
    const char* data() const {return data_;}
#endif
    int32_t length() const {
#ifdef WITH_BITSET
      if (bitset_ == nullptr)
        return parent_length_;
      return bitset_->length_;
#else
      return length_;
#endif
    }

    char* to_char() {
#ifdef WITH_BITSET
      return bitset_->to_char();
#else
      return data_;
#endif
    }

    int32_t parent_length() const {return parent_length_;}

    static inline int32_t nb_blocks(int32_t length);

    void ABCDE_to_ADCBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                             int32_t pos_E);
    void ABCDE_to_ADBpCpE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                               int32_t pos_E);
    void ABCDE_to_ACpDpBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                               int32_t pos_E);

    std::list<MutationEvent*> mutation_list;

#ifdef WITH_BITSET
    BitSet_SIMD* bitset_ = nullptr;
#else
    char* data_;
    int32_t length_;
#endif

    // Stats
    int32_t nb_swi_ = 0;
    int32_t nb_indels_= 0;
    int32_t nb_mut_= 0;

    int32_t nb_large_dupl_= 0;
    int32_t nb_large_del_= 0;
    int32_t nb_large_trans_= 0;
    int32_t nb_large_inv_= 0;
    int32_t nb_rear_= 0;

 private:

    int32_t parent_length_;
    int32_t nb_blocks_;
    Internal_SIMD_Struct* indiv_;

};

int32_t Dna_SIMD::nb_blocks(int32_t length) {
  return length/BLOCK_SIZE + 1;
}
}


#endif //RAEVOL_CUDA_DNA_SIMD_H
