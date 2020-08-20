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


#include "Rna_7.h"


double Rna_7::affinity_with_protein( int32_t index, pProtein *protein, Internal_SIMD_Struct* indiv, ExpManager* exp_m ) {
  int32_t len = protein->protein_length;

  if (len > 5) {
    double max = 0;
    double temp = 1;

    int32_t quadon_tab[5];

    for (int32_t pos = index; pos < index+5; pos++) {

      int8_t quadon[4];

      if (leading_lagging == LEADING) {
        for (int8_t i = 0; i < QUADON_SIZE; i++) {
          quadon[i] = (indiv->dna_->get_lead(pos + i) == '1')
                      ? 1 << (QUADON_SIZE - i - 1)
                      : 0;
        }
      } else {
        for (int8_t i = 0; i < QUADON_SIZE; i++) {
          quadon[i] = (indiv->dna_->get_lag(pos - i) != '1')
                      ? (QUADON_SIZE - i - 1)
                      : 0;
        }
      }

      quadon_tab[pos - index] = quadon[0] + quadon[1] + quadon[2] + quadon[3];
    }

    for (int32_t i = 0; i < len - 4; i++) {
      temp = 1;

      for (int8_t j = 0; j < 5; j++) {
        temp *= exp_m->exp_s()->get_binding_matrix(quadon_tab[j],
                                                   protein->codon_list[i + j]);
      }

      max = (max < temp) ? temp : max;
    }

    return max;
  } else {
    return 0.0;
  }
}
