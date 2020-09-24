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

#ifndef AEVOL_PROTEIN_7_H
#define AEVOL_PROTEIN_7_H


#include <cstdint>
#include "Rna_7.h"

namespace aevol {

class Protein_7 {
 public:
  Protein_7(){};
  Protein_7(int32_t t_protein_start,
            int32_t t_protein_end,
            int32_t t_protein_length,
            int8_t t_leading_lagging,
            double t_e, Rna_7* rna) {
    protein_start   = t_protein_start;
    protein_end     = t_protein_end;
    protein_length  = t_protein_length;
    leading_lagging = t_leading_lagging;
    e               = t_e;
    is_init_        = true;

#ifdef __REGUL
    rna_list_.push_back(rna);
      initial_e_ = e;
#endif
  }

#ifdef __REGUL
  Protein_7(Protein_R* prot_sig);


  Protein_7(Protein_7* prot);
#endif

  bool operator<(const Protein_7& other);

  int32_t protein_start;
  int32_t protein_end;
  int32_t protein_length;
  int8_t leading_lagging; // 0 = leading, 1 = lagging
  double m;
  double w;
  double h;
  double e;
  bool is_functional;

  bool is_init_ = false;


  int8_t codon_list[64*3] = {};
  int16_t nb_codons_ = 0;

#ifdef __REGUL
  bool is_TF_;

  double    delta_concentration_;
  bool      inherited_ = false;
  bool      signal_ = false;
#endif
  std::list<Rna_7*> rna_list_;
  double initial_e_ = -1;
};
}

#endif //AEVOL_PROTEIN_7_H
