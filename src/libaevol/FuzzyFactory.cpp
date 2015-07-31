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
#include "FuzzyFactory.h"
#include "ExpSetup.h"

namespace aevol {


//##############################################################################
//                                                                             #
//                           Class IndividualFactory                           #
//                                                                             #
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
/**
 * Create an individual with random sequences
 */
AbstractFuzzy* FuzzyFactory::create_fuzzy(ExpSetup* exp_s)
{
  AbstractFuzzy* fuzzy;

  if (exp_s->get_fuzzy_flavor() == 0) {
    fuzzy = new Fuzzy();
  } else {
    fuzzy = new HybridFuzzy();
  }

  return fuzzy;
}

// ============================================================================
//                            Non inline accessors
// ============================================================================
} // namespace aevol
