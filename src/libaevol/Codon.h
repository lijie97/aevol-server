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

#ifndef AEVOL_CODON_H_
#define AEVOL_CODON_H_

// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "Dna.h"
#include "macros.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================

class Codon {
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  Codon();
  explicit Codon(const Codon &model);
  explicit Codon(int8_t value);
  Codon(Dna* genome, Strand strand, int32_t index);
  explicit Codon(gzFile backup_file);

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~Codon();

  // =================================================================
  //                              Accessors
  // =================================================================
  int8_t value() { return value_; }

  // =================================================================
  //                            Public Methods
  // =================================================================
  bool is_start() { return value_ == CODON_START; }
  bool is_stop() { return value_ == CODON_STOP; }
  Codon* copy() { return new Codon(value_); } // TODO(dpa) use copy ctor instead!
  void save(gzFile backup_file);

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  int8_t value_;
};

// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_CODON_H_
