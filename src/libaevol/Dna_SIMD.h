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
    Dna_SIMD(Dna* dna);
    ~Dna_SIMD();

    void get(int32_t pos);
    void set(int32_t pos, char *value);
    void flip(int32_t pos);

    void remove(int32_t first, int32_t last);
    void insert(int32_t pos, const char* seq, int32_t seq_length = -1);
    void replace(int32_t pos, char* seq, int32_t seq_length = -1);

    const char* data() const {return data_;}
    int32_t length() const {return length_;}

    static inline int32_t nb_blocks(int32_t length);
 private:
    char* data_;
    int32_t length_;
    int32_t nb_blocks_;

    std::set<MutationEvent> mutation_list;
};

int32_t Dna_SIMD::nb_blocks(int32_t length) {
  return length/BLOCK_SIZE + 1;
}
}


#endif //RAEVOL_CUDA_DNA_SIMD_H
