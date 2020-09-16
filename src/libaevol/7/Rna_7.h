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


#ifndef AEVOL_RNA_7_H
#define AEVOL_RNA_7_H


#include <list>
#include <cstdint>

#include "ae_enums.h"
#include "Individual_7.h"
#include "ExpManager.h"
namespace aevol {
class Protein_7;

class AffinityFactor {
 public:
  AffinityFactor(double pconcentration, double efactor, double ofactor) {
    protein_concentration = pconcentration;
    enhancer_factor = efactor;
    operator_factor = ofactor;
  }

  double concentration();

  double protein_concentration;
  double enhancer_factor;
  double operator_factor;
};


class Rna_7 {
 public:
  Rna_7(){};
  Rna_7(int32_t t_begin,
        int32_t t_end,
        int8_t t_leading_lagging,
        double t_e,
        int32_t t_length) {
    begin             = t_begin;
    end               = t_end;
    leading_lagging   = t_leading_lagging;
    e                 = t_e;
    length            = t_length;
    is_coding_        = false;
    is_init_          = true;
    start_prot_count_ = 0;
  }

  ~Rna_7() {}


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
  double compute_synthesis_rate(Individual_7* indiv);

  double affinity_with_protein( int32_t index, Protein_7* protein,
                                     Individual_7* indiv,
                               ExpManager* exp_m);
#endif


  int32_t begin;
  int32_t end;
  int8_t leading_lagging; // 0 = leading, 1 = lagging
  double e;
  std::list<int32_t> start_prot;
  int32_t start_prot_count_;
  int32_t length;
  bool is_coding_;

  bool is_init_ = false;
};

}

#endif //AEVOL_RNA_7_H
