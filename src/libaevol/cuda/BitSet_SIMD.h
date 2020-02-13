//
// Created by arrouan on 24/01/18.
//

#ifndef RAEVOL_CUDA_BITSET_SIMD_H
#define RAEVOL_CUDA_BITSET_SIMD_H

#include <cstdint>
#include <cstdlib>

#include <cstdio>
#include <bitset>
#include <assert.h>
#ifdef _DYNAMIC_BITSET
#include <boost/dynamic_bitset.hpp>
#endif

namespace aevol {

#define BITSET_BLOCK_SIZE INT32_C(1024)
#define CHAR_TO_BITS 8

#ifdef _DYNAMIC_CUSTOM_BITSET
#define INT32_TO_BITS 32

#define ARRAY_CHECK_BIT(var,pos,type) (var[(pos/type)] & (1 << (pos%type)))
#define CHECK_BIT(var,pos) (var & (1 << (pos)))
#endif

constexpr const int32_t BITSET_PROM_SEQ_LEAD = 0b0110100100111001101010;
constexpr const int32_t BITSET_PROM_SEQ_LAG  = 0b1001011011000110010101;


constexpr const int32_t BITSET_SHINE_DAL_SEQ_LEAD = 0b000110110;
constexpr const int32_t BITSET_SHINE_DAL_SEQ_LAG  = 0b111001001;

constexpr const char BITSET_PROTEIN_END_LEAD  = 0b100;
constexpr const char BITSET_PROTEIN_END_LAG   = 0b011;

constexpr int8_t BITSET_CODON_SIZE  = 3;

constexpr int32_t BITSET_STATIC_MAX_SIZE  = 10000000;



class BitSet_SIMD {
 public:
    BitSet_SIMD(BitSet_SIMD* bitset);

    BitSet_SIMD(const char* data, int32_t length);

    BitSet_SIMD(int32_t length);

    ~BitSet_SIMD() {
#ifdef _DYNAMIC_CUSTOM_BITSET
      free(data_);
#endif
       }

#ifdef _DYNAMIC_CUSTOM_BITSET
    __inline void set_to_1(int32_t position) {
      __assume_aligned(data_, 64);
      data_[(position / CHAR_TO_BITS)] |= (1 << (position % CHAR_TO_BITS));
    }

    __inline void set_to_0(int32_t position) {
      __assume_aligned(data_, 64);
      data_[(position / CHAR_TO_BITS)] &= ~(1 << (position % CHAR_TO_BITS));
    }

    __inline bool get(int32_t position) {
      //__assume_aligned(data_, 64);
        if (position >= length_ || position < 0) {
            printf("Position %d Length %d\n",position,length_);

        }
      return data_[(position/CHAR_TO_BITS)] & (1 << (position%CHAR_TO_BITS));
    }


#endif

    __inline void flip(int32_t position) {
#ifdef _DYNAMIC_CUSTOM_BITSET
      __assume_aligned(data_, 64);
      data_[(position / CHAR_TO_BITS)] ^= (1 << (position % CHAR_TO_BITS));
#elif _STATIC_BITSET
      data_[position].flip();
#elif _DYNAMIC_BITSET
      data_[position].flip();
#endif
    }


    __inline void flip_all() {
#ifdef _DYNAMIC_CUSTOM_BITSET
      __assume_aligned(data_, 64);
      for (int32_t position = 0; position < length_; position++)
        data_[(position / CHAR_TO_BITS)] ^= (1 << (position % CHAR_TO_BITS));
#elif _STATIC_BITSET
      data_.flip();
#endif
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
      //__declspec(align(32))
      int8_t dist[22];

#ifdef _STATIC_BITSET
      if (LEADING) {
        dist[0] = get(BITSET_PROM_SEQ_LEAD, 0) == data_[pos + 0 >= length_ ? pos + 0 - length_ : pos + 0] ? 0 : 1;
        dist[1] = get(BITSET_PROM_SEQ_LEAD, 1) == data_[pos + 1 >= length_ ? pos + 1 - length_ : pos + 1] ? 0 : 1;
        dist[2] = get(BITSET_PROM_SEQ_LEAD, 2) == data_[pos + 2 >= length_ ? pos + 2 - length_ : pos + 2] ? 0 : 1;
        dist[3] = get(BITSET_PROM_SEQ_LEAD, 3) == data_[pos + 3 >= length_ ? pos + 3 - length_ : pos + 3] ? 0 : 1;
        dist[4] = get(BITSET_PROM_SEQ_LEAD, 4) == data_[pos + 4 >= length_ ? pos + 4 - length_ : pos + 4] ? 0 : 1;
        dist[5] = get(BITSET_PROM_SEQ_LEAD, 5) == data_[pos + 5 >= length_ ? pos + 5 - length_ : pos + 5] ? 0 : 1;
        dist[6] = get(BITSET_PROM_SEQ_LEAD, 6) == data_[pos + 6 >= length_ ? pos + 6 - length_ : pos + 6] ? 0 : 1;
        dist[7] = get(BITSET_PROM_SEQ_LEAD, 7) == data_[pos + 7 >= length_ ? pos + 7 - length_ : pos + 7] ? 0 : 1;
        dist[8] = get(BITSET_PROM_SEQ_LEAD, 8) == data_[pos + 8 >= length_ ? pos + 8 - length_ : pos + 8] ? 0 : 1;
        dist[9] = get(BITSET_PROM_SEQ_LEAD, 9) == data_[pos + 9 >= length_ ? pos + 9 - length_ : pos + 9] ? 0 : 1;
        dist[10] = get(BITSET_PROM_SEQ_LEAD, 10) == data_[pos + 10 >= length_ ? pos + 10 - length_ : pos + 10] ? 0 : 1;
        dist[11] = get(BITSET_PROM_SEQ_LEAD, 11) == data_[pos + 11 >= length_ ? pos + 11 - length_ : pos + 11] ? 0 : 1;
        dist[12] = get(BITSET_PROM_SEQ_LEAD, 12) == data_[pos + 12 >= length_ ? pos + 12 - length_ : pos + 12] ? 0 : 1;
        dist[13] = get(BITSET_PROM_SEQ_LEAD, 13) == data_[pos + 13 >= length_ ? pos + 13 - length_ : pos + 13] ? 0 : 1;
        dist[14] = get(BITSET_PROM_SEQ_LEAD, 14) == data_[pos + 14 >= length_ ? pos + 14 - length_ : pos + 14] ? 0 : 1;
        dist[15] = get(BITSET_PROM_SEQ_LEAD, 15) == data_[pos + 15 >= length_ ? pos + 15 - length_ : pos + 15] ? 0 : 1;
        dist[16] = get(BITSET_PROM_SEQ_LEAD, 16) == data_[pos + 16 >= length_ ? pos + 16 - length_ : pos + 16] ? 0 : 1;
        dist[17] = get(BITSET_PROM_SEQ_LEAD, 17) == data_[pos + 17 >= length_ ? pos + 17 - length_ : pos + 17] ? 0 : 1;
        dist[18] = get(BITSET_PROM_SEQ_LEAD, 18) == data_[pos + 18 >= length_ ? pos + 18 - length_ : pos + 18] ? 0 : 1;
        dist[19] = get(BITSET_PROM_SEQ_LEAD, 19) == data_[pos + 19 >= length_ ? pos + 19 - length_ : pos + 19] ? 0 : 1;
        dist[20] = get(BITSET_PROM_SEQ_LEAD, 20) == data_[pos + 20 >= length_ ? pos + 20 - length_ : pos + 20] ? 0 : 1;
        dist[21] = get(BITSET_PROM_SEQ_LEAD, 21) == data_[pos + 21 >= length_ ? pos + 21 - length_ : pos + 21] ? 0 : 1;
      } else {
        dist[0] = get(BITSET_PROM_SEQ_LAG, 0) == data_[pos - 0 < 0 ? length_ + pos - 0 : pos - 0] ? 0 : 1;
        dist[1] = get(BITSET_PROM_SEQ_LAG, 1) == data_[pos - 1 < 0 ? length_ + pos - 1 : pos - 1] ? 0 : 1;
        dist[2] = get(BITSET_PROM_SEQ_LAG, 2) == data_[pos - 2 < 0 ? length_ + pos - 2 : pos - 2] ? 0 : 1;
        dist[3] = get(BITSET_PROM_SEQ_LAG, 3) == data_[pos - 3 < 0 ? length_ + pos - 3 : pos - 3] ? 0 : 1;
        dist[4] = get(BITSET_PROM_SEQ_LAG, 4) == data_[pos - 4 < 0 ? length_ + pos - 4 : pos - 4] ? 0 : 1;
        dist[5] = get(BITSET_PROM_SEQ_LAG, 5) == data_[pos - 5 < 0 ? length_ + pos - 5 : pos - 5] ? 0 : 1;
        dist[6] = get(BITSET_PROM_SEQ_LAG, 6) == data_[pos - 6 < 0 ? length_ + pos - 6 : pos - 6] ? 0 : 1;
        dist[7] = get(BITSET_PROM_SEQ_LAG, 7) == data_[pos - 7 < 0 ? length_ + pos - 7 : pos - 7] ? 0 : 1;
        dist[8] = get(BITSET_PROM_SEQ_LAG, 8) == data_[pos - 8 < 0 ? length_ + pos - 8 : pos - 8] ? 0 : 1;
        dist[9] = get(BITSET_PROM_SEQ_LAG, 9) == data_[pos - 9 < 0 ? length_ + pos - 9 : pos - 9] ? 0 : 1;
        dist[10] = get(BITSET_PROM_SEQ_LAG, 10) == data_[pos - 10 < 0 ? length_ + pos - 10 : pos - 10] ? 0 : 1;
        dist[11] = get(BITSET_PROM_SEQ_LAG, 11) == data_[pos - 11 < 0 ? length_ + pos - 11 : pos - 11] ? 0 : 1;
        dist[12] = get(BITSET_PROM_SEQ_LAG, 12) == data_[pos - 12 < 0 ? length_ + pos - 12 : pos - 12] ? 0 : 1;
        dist[13] = get(BITSET_PROM_SEQ_LAG, 13) == data_[pos - 13 < 0 ? length_ + pos - 13 : pos - 13] ? 0 : 1;
        dist[14] = get(BITSET_PROM_SEQ_LAG, 14) == data_[pos - 14 < 0 ? length_ + pos - 14 : pos - 14] ? 0 : 1;
        dist[15] = get(BITSET_PROM_SEQ_LAG, 15) == data_[pos - 15 < 0 ? length_ + pos - 15 : pos - 15] ? 0 : 1;
        dist[16] = get(BITSET_PROM_SEQ_LAG, 16) == data_[pos - 16 < 0 ? length_ + pos - 16 : pos - 16] ? 0 : 1;
        dist[17] = get(BITSET_PROM_SEQ_LAG, 17) == data_[pos - 17 < 0 ? length_ + pos - 17 : pos - 17] ? 0 : 1;
        dist[18] = get(BITSET_PROM_SEQ_LAG, 18) == data_[pos - 18 < 0 ? length_ + pos - 18 : pos - 18] ? 0 : 1;
        dist[19] = get(BITSET_PROM_SEQ_LAG, 19) == data_[pos - 19 < 0 ? length_ + pos - 19 : pos - 19] ? 0 : 1;
        dist[20] = get(BITSET_PROM_SEQ_LAG, 20) == data_[pos - 20 < 0 ? length_ + pos - 20 : pos - 20] ? 0 : 1;
        dist[21] = get(BITSET_PROM_SEQ_LAG, 21) == data_[pos - 21 < 0 ? length_ + pos - 21 : pos - 21] ? 0 : 1;
      }
#endif

#ifdef _DYNAMIC_CUSTOM_BITSET
      for (int8_t motif_id = 0; motif_id < 22; motif_id++) {
          if (LEADING)
            dist[motif_id] = get(BITSET_PROM_SEQ_LEAD,motif_id) == get(pos + motif_id >= length_ ? pos +
                                                                                                   motif_id -
                                                                                                   length_
                                                                                                 : pos +
                                                                                                   motif_id) ? 0 : 1;
          else
            dist[motif_id] = get(BITSET_PROM_SEQ_LAG,motif_id) == get(pos - motif_id < 0 ? length_ +
                                                                                                               pos -
                                                                                                               motif_id :
                                                                                    pos - motif_id) ? 0 : 1;
      }
#endif


      return dist[0] + dist[1]+ dist[2]+ dist[3]+ dist[4]+ dist[5]+ dist[6]
             + dist[7] + dist[8]+ dist[9]+ dist[10]+ dist[11]+ dist[12]+ dist[13]
             + dist[14]+ dist[15]+ dist[16]+ dist[17]+ dist[18]+ dist[19]+ dist[20]
             + dist[21];

    }

    __inline bool is_terminator(bool LEADING, int pos) {
      bool dist[4];


#ifdef _STATIC_BITSET
      if (LEADING) {
        dist[0] = data_[pos + 0 >= length_ ? pos + 0 - length_ : pos + 0] ^ data_[pos - 0 + 10 >= length_ ? pos - 0 + 10 - length_ : pos - 0 + 10];
        dist[1] = data_[pos + 1 >= length_ ? pos + 1 - length_ : pos + 1] ^ data_[pos - 1 + 10 >= length_ ? pos - 1 + 10 - length_ : pos - 1 + 10];
        dist[2] = data_[pos + 2 >= length_ ? pos + 2 - length_ : pos + 2] ^ data_[pos - 2 + 10 >= length_ ? pos - 2 + 10 - length_ : pos - 2 + 10];
        dist[3] = data_[pos + 3 >= length_ ? pos + 3 - length_ : pos + 3] ^ data_[pos - 3 + 10 >= length_ ? pos - 3 + 10 - length_ : pos - 3 + 10];
      } else {
        dist[0] = data_[pos - 0 < 0 ? pos - 0 + length_ : pos - 0] ^ data_[pos + 0 - 10 < 0 ? pos + 0 - 10 + length_ : pos + 0 - 10];
        dist[1] = data_[pos - 1 < 0 ? pos - 1 + length_ : pos - 1] ^ data_[pos + 1 - 10 < 0 ? pos + 1 - 10 + length_ : pos + 1 - 10];
        dist[2] = data_[pos - 2 < 0 ? pos - 2 + length_ : pos - 2] ^ data_[pos + 2 - 10 < 0 ? pos + 2 - 10 + length_ : pos + 2 - 10];
        dist[3] = data_[pos - 3 < 0 ? pos - 3 + length_ : pos - 3] ^ data_[pos + 3 - 10 < 0 ? pos + 3 - 10 + length_ : pos + 3 - 10];
      }
#endif
#ifdef _DYNAMIC_CUSTOM_BITSET
      #pragma ivdep
      for (int8_t motif_id = 0; motif_id < 4; motif_id++) {
        if (LEADING) {
            /*printf("Read %d  and %d -- (length %d)\n",pos + motif_id >= length_ ? pos + motif_id - length_
                                                                     : pos + motif_id,pos - motif_id + 10 >= length_ ? pos - motif_id + 10 - length_
                                                                                                                     : pos - motif_id + 10,length_);
*/
            dist[motif_id] = get(pos + motif_id >= length_ ? pos + motif_id - length_
                                                           : pos + motif_id) !=
                             get(pos - motif_id + 10 >= length_ ? pos - motif_id + 10 - length_
                                                                : pos - motif_id + 10) ? true : false;
        }
        else
          dist[motif_id] = get(pos - motif_id < 0 ? pos - motif_id + length_
                                   : pos - motif_id) !=
                           get(pos + motif_id - 10 < 0 ? pos + motif_id - 10 + length_
                                   : pos + motif_id - 10) ? true : false;

      }
#endif

      return dist[0] && dist[1] && dist[2] && dist[3];
    }

    __inline bool is_shine_dalgarno_protein_start(bool LEADING, int32_t pos) {

      bool start[9] = {false, false, false,
                       false, false, false,
                       false, false, false};

#ifdef _STATIC_BITSET
#pragma omp simd
      for (int32_t k = 0; k < 9; k++) {
        int32_t k_t = k >= 6 ? k + 4 : k;
        if (LEADING)
          (data_[pos + k_t >= length_ ? pos + k_t - length_ : pos + k_t] ==
           get(BITSET_SHINE_DAL_SEQ_LEAD,k)) ? (start[k] = true) : (start[k] = false);
        else
          (data_[pos - k_t < 0 ? length_ + (pos - k_t) : pos - k_t] ==
           get(BITSET_SHINE_DAL_SEQ_LAG,k)) ? (start[k] = true) : (start[k] = false);
      }
      /*
      if (LEADING) {
        (data_[pos + 0 >= length_ ? pos + 0 - length_ : pos + 0] ==
         get(BITSET_SHINE_DAL_SEQ_LEAD, 0)) ? (start[0] = true)
                                            : (start[0] = false);
        (data_[pos + 1 >= length_ ? pos + 1 - length_ : pos + 1] ==
         get(BITSET_SHINE_DAL_SEQ_LEAD, 1)) ? (start[1] = true)
                                            : (start[1] = false);
        (data_[pos + 2 >= length_ ? pos + 2 - length_ : pos + 2] ==
         get(BITSET_SHINE_DAL_SEQ_LEAD, 2)) ? (start[2] = true)
                                            : (start[2] = false);
        (data_[pos + 3 >= length_ ? pos + 3 - length_ : pos + 3] ==
         get(BITSET_SHINE_DAL_SEQ_LEAD, 3)) ? (start[3] = true)
                                            : (start[3] = false);
        (data_[pos + 4 >= length_ ? pos + 4 - length_ : pos + 4] ==
         get(BITSET_SHINE_DAL_SEQ_LEAD, 4)) ? (start[4] = true)
                                            : (start[4] = false);
        (data_[pos + 5 >= length_ ? pos + 5 - length_ : pos + 5] ==
         get(BITSET_SHINE_DAL_SEQ_LEAD, 5)) ? (start[5] = true)
                                            : (start[5] = false);

        (data_[pos + 10 >= length_ ? pos + 10 - length_ : pos + 10] ==
         get(BITSET_SHINE_DAL_SEQ_LEAD, 6)) ? (start[6] = true)
                                            : (start[6] = false);
        (data_[pos + 11 >= length_ ? pos + 11 - length_ : pos + 11] ==
         get(BITSET_SHINE_DAL_SEQ_LEAD, 7)) ? (start[7] = true)
                                            : (start[7] = false);
        (data_[pos + 12 >= length_ ? pos + 12 - length_ : pos + 12] ==
         get(BITSET_SHINE_DAL_SEQ_LEAD, 8)) ? (start[8] = true)
                                            : (start[8] = false);
      } else {
        (data_[pos - 0 < 0 ? length_ + (pos - 0) : pos - 0] ==
         get(BITSET_SHINE_DAL_SEQ_LAG, 0)) ? (start[0] = true)
                                           : (start[0] = false);
        (data_[pos - 1 < 0 ? length_ + (pos - 1) : pos - 1] ==
         get(BITSET_SHINE_DAL_SEQ_LAG, 1)) ? (start[1] = true)
                                           : (start[1] = false);
        (data_[pos - 2 < 0 ? length_ + (pos - 2) : pos - 2] ==
         get(BITSET_SHINE_DAL_SEQ_LAG, 2)) ? (start[2] = true)
                                           : (start[2] = false);
        (data_[pos - 3 < 0 ? length_ + (pos - 3) : pos - 3] ==
         get(BITSET_SHINE_DAL_SEQ_LAG, 3)) ? (start[3] = true)
                                           : (start[3] = false);
        (data_[pos - 4 < 0 ? length_ + (pos - 4) : pos - 4] ==
         get(BITSET_SHINE_DAL_SEQ_LAG, 4)) ? (start[4] = true)
                                           : (start[4] = false);
        (data_[pos - 5 < 0 ? length_ + (pos - 5) : pos - 5] ==
         get(BITSET_SHINE_DAL_SEQ_LAG, 5)) ? (start[5] = true)
                                           : (start[5] = false);

        (data_[pos - 10 < 0 ? length_ + (pos - 10) : pos - 10] ==
         get(BITSET_SHINE_DAL_SEQ_LAG, 6)) ? (start[6] = true)
                                           : (start[6] = false);
        (data_[pos - 11 < 0 ? length_ + (pos - 11) : pos - 11] ==
         get(BITSET_SHINE_DAL_SEQ_LAG, 7)) ? (start[7] = true)
                                           : (start[7] = false);
        (data_[pos - 12 < 0 ? length_ + (pos - 12) : pos - 12] ==
         get(BITSET_SHINE_DAL_SEQ_LAG, 8)) ? (start[8] = true)
                                           : (start[8] = false);
      }*/
#endif

#ifdef _DYNAMIC_CUSTOM_BITSET
      for (int k = 0; k < 9; k++) {
        int k_t = k >= 6 ? k + 4 : k;
        if (LEADING)
          (get(pos + k_t >= length_ ? pos + k_t - length_ : pos + k_t) ==
              get(BITSET_SHINE_DAL_SEQ_LEAD,k)) ? (start[k] = true) : (start[k] = false);
        else
          (get(pos - k_t < 0 ? length_ + (pos - k_t) : pos - k_t) ==
           get(BITSET_SHINE_DAL_SEQ_LAG,k)) ? (start[k] = true) : (start[k] = false);
      }
#endif

      return start[0]&&start[1]&&start[2]&&
          start[3]&&start[4]&&start[5]&&
          start[6]&&start[7]&&start[8];
    }

    __inline bool is_protein_stop(bool LEADING, int pos) {

      bool is_protein[3] = {false,false,false};

#ifdef _STATIC_BITSET
      if (LEADING) {
        (data_[pos + 0 >= length_ ? pos - length_ + 0 : pos + 0] ==
         get(BITSET_PROTEIN_END_LEAD, 0)) ? (is_protein[0] = true) :
        (is_protein[0] = false);
        (data_[pos + 1 >= length_ ? pos - length_ + 1 : pos + 1] ==
         get(BITSET_PROTEIN_END_LEAD, 1)) ? (is_protein[1] = true) :
        (is_protein[1] = false);
        (data_[pos + 2 >= length_ ? pos - length_ + 2 : pos + 2] ==
         get(BITSET_PROTEIN_END_LEAD, 2)) ? (is_protein[2] = true) :
        (is_protein[2] = false);
      } else {
        (data_[pos - 0 < 0 ? length_ + (pos - 0) : pos - 0] ==
         get(BITSET_PROTEIN_END_LAG, 0)) ? (is_protein[0] = true) :
        (is_protein[0] = false);
        (data_[pos - 1 < 0 ? length_ + (pos - 1) : pos - 1] ==
         get(BITSET_PROTEIN_END_LAG, 1)) ? (is_protein[1] = true) :
        (is_protein[1] = false);
        (data_[pos - 2 < 0 ? length_ + (pos - 2) : pos - 2] ==
         get(BITSET_PROTEIN_END_LAG, 2)) ? (is_protein[2] = true) :
        (is_protein[2] = false);
      }
#endif

#ifdef _DYNAMIC_CUSTOM_BITSET
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
#endif
      return is_protein[0]&&is_protein[1]&&is_protein[2];
    }

    __inline int8_t extract_codon(bool LEADING, int pos) {
      int8_t value = 0;

#ifdef _STATIC_BITSET
      if (LEADING) {
        if (data_[pos + 0 >= length_ ? pos + 0 - length_ : pos + 0])
          value += 1 << (BITSET_CODON_SIZE - 0 - 1);
        if (data_[pos + 1 >= length_ ? pos + 1 - length_ : pos + 1])
          value += 1 << (BITSET_CODON_SIZE - 1 - 1);
        if (data_[pos + 2 >= length_ ? pos + 2 - length_ : pos + 2])
          value += 1 << (BITSET_CODON_SIZE - 2 - 1);
      } else {
        if (!data_[pos - 0 < 0 ? length_ + (pos - 0) : pos - 0])
          value += 1 << (BITSET_CODON_SIZE - 0 - 1);
        if (!data_[pos - 1 < 0 ? length_ + (pos - 1) : pos - 1])
          value += 1 << (BITSET_CODON_SIZE - 1 - 1);
        if (!data_[pos - 2 < 0 ? length_ + (pos - 2) : pos - 2])
          value += 1 << (BITSET_CODON_SIZE - 2 - 1);
      }
#endif

#ifdef _DYNAMIC_CUSTOM_BITSET
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
#endif

      return value;
    }

    void print();

    char* to_char() {
#ifdef _DYNAMIC_CUSTOM_BITSET
      char* ret = new char[length_+1];

      for (int32_t i = 0; i < length_; i++) {
        if (get(i))
          ret[i] = '1';
        else
          ret[i] = '0';
      }

      ret[length_] = '\0';
      return ret;
#elif _STATIC_BITSET
      char *ret = (char*)data_.to_string().c_str();
      return ret;
#else
      return nullptr;
#endif
    }

#ifdef _DYNAMIC_CUSTOM_BITSET
    char *data_;
#elif _STATIC_BITSET
    std::bitset<BITSET_STATIC_MAX_SIZE> data_;
#endif

    int32_t length_;
    int32_t nb_block_;
    int32_t nb_char_block_;


};

}


#endif //RAEVOL_CUDA_BITSET_SIMD_H
