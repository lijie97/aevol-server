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
//*****************************************************************************




// ============================================================================
//                                   Includes
// ============================================================================
#include "habitat.h"


namespace aevol {


//##############################################################################
//                                                                             #
//                                Class Habitat                                #
//                                                                             #
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
Habitat::Habitat(void) {
  compound_amount_ = 0.0;
}

Habitat::Habitat(gzFile backup_file) {
  gzread(backup_file, &compound_amount_, sizeof(compound_amount_));
}

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
void Habitat::save(gzFile backup_file) const {
  gzwrite(backup_file, &compound_amount_, sizeof(compound_amount_));
}

// ============================================================================
//                            Non inline accessors
// ============================================================================
} // namespace aevol
