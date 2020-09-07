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

namespace aevol {
class Protein_7 {
 public:
  Protein_7(){};
  Protein_7(int32_t t_protein_start,
            int32_t t_protein_end,
            int32_t t_protein_length,
            int8_t t_leading_lagging,
            double t_e) {
    protein_start   = t_protein_start;
    protein_end     = t_protein_end;
    protein_length  = t_protein_length;
    leading_lagging = t_leading_lagging;
    e               = t_e;
    is_init_        = true;
  }

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
};
}

#endif //AEVOL_PROTEIN_7_H
