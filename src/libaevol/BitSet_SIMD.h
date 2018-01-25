//
// Created by arrouan on 24/01/18.
//

#ifndef RAEVOL_CUDA_BITSET_SIMD_H
#define RAEVOL_CUDA_BITSET_SIMD_H

#include <cstdint>

namespace aevol {

#define BITSET_BLOCK_SIZE INT32_C(1024)

class BitSet_SIMD {
 public:
    BitSet_SIMD(BitSet_SIMD* bitset);

    BitSet_SIMD(const char* data, int32_t length);

    inline void set_to_1(int position) {
      data_[(position / sizeof(char))] |= (1 << (position % sizeof(char)));
    }

    inline void set_to_0(int position) {
      data_[(position / sizeof(char))] &= ~(1 << (position % sizeof(char)));
    }

    inline void flip(int position) {
      data_[(position / sizeof(char))] ^= (1 << (position % sizeof(char)));
    }

    inline int get(int position) {
      return data_[(position/sizeof(char))] & (1 << (position%sizeof(char)));
    }

    void append(BitSet_SIMD* bitset);

    void insert_at(BitSet_SIMD* bitset, int position);

    BitSet_SIMD* extract(int pos_1, int pos2);

    void print();

    char *data_;

    int32_t length_;
    int32_t nb_block_;
    int32_t nb_char_block_;
};

}


#endif //RAEVOL_CUDA_BITSET_SIMD_H
