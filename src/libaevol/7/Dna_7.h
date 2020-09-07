// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

#ifndef AEVOL_DNA_7_H
#define AEVOL_DNA_7_H


#include <cstdint>
#include <list>

#include "Utils.h"
#include "MutationEvent.h"

namespace aevol {

#define BLOCK_SIZE INT32_C(1024)
class Dna;
class Individual_7;
class DnaFactory;

class Dna_7 {
 public:
  Dna_7(Dna* dna, Individual_7* indiv, DnaFactory* dna_factory);
  Dna_7(Dna_7* dna, Individual_7* indiv, DnaFactory* dna_factory);
  Dna_7(Dna* dna, DnaFactory* dna_factory);
  Dna_7(int length, DnaFactory* dna_factory);
    ~Dna_7();

    void set_indiv(Dna_7* dna, Individual_7* indiv);
    void set_indiv(Dna* dna, DnaFactory* dna_factory);
    void set_indiv(int length, int parent_length, Individual_7* indiv);

    inline char get_lead(int32_t pos) { return data_[Utils::mod(pos,length_ )]; }
    inline char get_lag(int32_t pos) { return data_[Utils::mod(pos,length_)]; }

    inline int nb_block() { return nb_blocks_; }

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
    Individual_7* indiv_;
    DnaFactory* dna_factory_;

};

int32_t Dna_7::nb_blocks(int32_t length) {
  return length/BLOCK_SIZE + 1;
}
}


#endif
