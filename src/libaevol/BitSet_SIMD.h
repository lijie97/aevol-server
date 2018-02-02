//
// Created by arrouan on 24/01/18.
//

#ifndef RAEVOL_CUDA_BITSET_SIMD_H
#define RAEVOL_CUDA_BITSET_SIMD_H

#include <cstdint>
#include <cstdlib>

#include <cstdio>

namespace aevol {

#define BITSET_BLOCK_SIZE INT32_C(1024)
#define CHAR_TO_BITS 8
#define INT32_TO_BITS 32

#define ARRAY_CHECK_BIT(var,pos,type) (var[(pos/type)] & (1 << (pos%type)))
#define CHECK_BIT(var,pos) (var & (1 << (pos)))

constexpr const int32_t BITSET_PROM_SEQ_LEAD = 0b0110100100111001101010;
constexpr const int32_t BITSET_PROM_SEQ_LAG  = 0b1001011011000110010101;


constexpr const int32_t BITSET_SHINE_DAL_SEQ_LEAD = 0b000110110;
constexpr const int32_t BITSET_SHINE_DAL_SEQ_LAG  = 0b111001001;

constexpr const char BITSET_PROTEIN_END_LEAD  = 0b100;
constexpr const char BITSET_PROTEIN_END_LAG   = 0b011;

constexpr int8_t BITSET_CODON_SIZE  = 3;

class BitSet_SIMD {
 public:
    BitSet_SIMD(BitSet_SIMD* bitset);

    BitSet_SIMD(const char* data, int32_t length);

    BitSet_SIMD(int32_t length);

    ~BitSet_SIMD() { free(data_); }

    __inline void set_to_1(int32_t position) {
      __assume_aligned(data_, 64);
      data_[(position / CHAR_TO_BITS)] |= (1 << (position % CHAR_TO_BITS));
    }

    __inline void set_to_0(int32_t position) {
      __assume_aligned(data_, 64);
      data_[(position / CHAR_TO_BITS)] &= ~(1 << (position % CHAR_TO_BITS));
    }

    __inline void flip(int32_t position) {
      __assume_aligned(data_, 64);
      data_[(position / CHAR_TO_BITS)] ^= (1 << (position % CHAR_TO_BITS));
    }

    __inline void flip_all() {
      __assume_aligned(data_, 64);
      for (int32_t position = 0; position < length_; position++)
        data_[(position / CHAR_TO_BITS)] ^= (1 << (position % CHAR_TO_BITS));
    }

    __inline bool get(int32_t position) {
      __assume_aligned(data_, 64);
      return data_[(position/CHAR_TO_BITS)] & (1 << (position%CHAR_TO_BITS));
    }

    __inline static bool get(char* data, int32_t position) {
      return data[(position/CHAR_TO_BITS)] & (1 << (position%CHAR_TO_BITS));
    }

    __inline static bool get(int32_t data, int32_t position) {
      return data & (1 << (position));
    }


    __inline static bool get(char data, int32_t position) {
      return data & (1 << (position));
    }

    void append(BitSet_SIMD* bitset);

    void insert_at(BitSet_SIMD* bitset, int32_t position);
    void insert_at(BitSet_SIMD* bitset, int32_t start, int32_t stop, int32_t position);

    void replace(BitSet_SIMD* bitset, int32_t position);
    void replace(aevol::BitSet_SIMD* bitset, int32_t start, int32_t stop, int32_t position);

    void remove(int32_t pos_1, int32_t pos_2);

    BitSet_SIMD* duplicate(int32_t pos_1, int32_t pos2, bool invert = false);

    __inline int8_t is_promoter(bool LEADING, int pos) {
      __declspec(align(32)) int8_t dist[22];


      for (int8_t motif_id = 0; motif_id < 22; motif_id++) {
          if (LEADING)
            dist[motif_id] = get(BITSET_PROM_SEQ_LEAD,motif_id) == get(pos + motif_id >= length_ ? pos +
                                                                                                   motif_id -
                                                                                                   length_
                                                                                                 : pos +
                                                                                                   motif_id) ? 0 : 1;
          else dist[motif_id] = get(BITSET_PROM_SEQ_LAG,motif_id) == get(pos - motif_id < 0 ? length_ +
                                                                                                               pos -
                                                                                                               motif_id :
                                                                                    pos - motif_id) ? 0 : 1;

      }

      return dist[0] + dist[1]+ dist[2]+ dist[3]+ dist[4]+ dist[5]+ dist[6]
             + dist[7] + dist[8]+ dist[9]+ dist[10]+ dist[11]+ dist[12]+ dist[13]
             + dist[14]+ dist[15]+ dist[16]+ dist[17]+ dist[18]+ dist[19]+ dist[20]
             + dist[21];
    }

    __inline bool is_terminator(bool LEADING, int pos) {
      bool dist[4];

      #pragma ivdep
      for (int8_t motif_id = 0; motif_id < 4; motif_id++) {
        if (LEADING) dist[motif_id] = get(pos + motif_id >= length_ ? pos + motif_id - length_
                                              : pos + motif_id) !=
                                      get(pos - motif_id + 10 >= length_ ? pos - motif_id + 10 - length_
                                              : pos - motif_id + 10) ? true : false;
        else
          dist[motif_id] = get(pos - motif_id < 0 ? pos - motif_id + length_
                                   : pos - motif_id) !=
                           get(pos + motif_id - 10 < 0 ? pos + motif_id - 10 + length_
                                   : pos + motif_id - 10) ? true : false;
      }

      return dist[0] && dist[1] && dist[2] && dist[3];
    }

    __inline bool is_shine_dalgarno_protein_start(bool LEADING, int pos) {
      bool start[9] = {false,false,false,
                       false,false,false,
                       false,false,false};
      for (int k = 0; k < 9; k++) {
        int k_t = k >= 6 ? k + 4 : k;
        if (LEADING)
          (get(pos + k_t >= length_ ? pos + k_t - length_ : pos + k_t) ==
              get(BITSET_SHINE_DAL_SEQ_LEAD,k)) ? (start[k] = true) : (start[k] = false);
        else
          (get(pos - k_t < 0 ? length_ + (pos - k_t) : pos - k_t) ==
           get(BITSET_SHINE_DAL_SEQ_LAG,k)) ? (start[k] = true) : (start[k] = false);
      }
      return start[0]&&start[1]&&start[2]&&
          start[3]&&start[4]&&start[5]&&
          start[6]&&start[7]&&start[8];
    }

    __inline bool is_protein_stop(bool LEADING, int pos) {
      bool is_protein[3] = {false,false,false};
      for (int k = 0; k < 3; k++) {
        if (LEADING)
          (get(pos + k >= length_ ? pos - length_ + k : pos + k) ==
              get(BITSET_PROTEIN_END_LEAD,k)) ? (is_protein[k] = true) :
              (is_protein[k] = false);
        else
          (get(pos - k < 0 ? length_ + (pos - k) : pos - k) ==
              get(BITSET_PROTEIN_END_LAG,k)) ? (is_protein[k] = true) :
              (is_protein[k] = false);
      }
      return is_protein[0]&&is_protein[1]&&is_protein[2];
    }

    __inline int8_t extract_codon(bool LEADING, int pos) {
      int8_t value = 0;

      if (LEADING)
        for (int8_t i = 0; i < 3; i++) {
          if (get(pos + i >= length_ ? pos + i - length_ : pos + i))
            value += 1 << (BITSET_CODON_SIZE - i - 1);
        }
      else
        for (int8_t i = 0; i < 3; i++) {
          if (!get(pos - i < 0 ? length_ + (pos - i) : pos - i))
            value += 1 << (BITSET_CODON_SIZE - i - 1);
        }

      return value;
    }

    void print();

    char *data_;

    int32_t length_;
    int32_t nb_block_;
    int32_t nb_char_block_;


};

}


#endif //RAEVOL_CUDA_BITSET_SIMD_H
