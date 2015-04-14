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


#ifndef __AE_PHENOTYPIC_TARGET_HANDLER_H__
#define __AE_PHENOTYPIC_TARGET_HANDLER_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>
#include <list>

#include "phenotypic_target.h"
#include "ae_gaussian.h"
#include "ae_enums.h"
#include "ae_jumping_mt.h"
#include "fuzzy.h"

using std::list;


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================





/**
 * Manages a phenotypic target and its "evolution" over time
 *
 * Handles a phenotypic target, the variation and/or noise that may be applied
 * to it as well as the set of possible phenotypic targets and the rules that
 * define how and when we switch from one to another
 */
class PhenotypicTargetHandler
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  PhenotypicTargetHandler(void) = default; //< Default ctor
  PhenotypicTargetHandler(const PhenotypicTargetHandler&) = delete; //< Copy ctor
  PhenotypicTargetHandler(PhenotypicTargetHandler&&) = delete; //< Move ctor

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~PhenotypicTargetHandler(void) = default; //< Destructor

  // ==========================================================================
  //                                 Getters
  // ==========================================================================

  // ==========================================================================
  //                                 Setters
  // ==========================================================================

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================





 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  // ------------------------------------------------ Current Phenotypic Target
  std::unique_ptr<PhenotypicTarget> phenotypic_target_;

  // ---------------------------------------------------------------- Gaussians
  /// Phenotypic target's constitutive Gaussians in their initial state
  std::list<ae_gaussian> initial_gaussians_;
    /// Phenotypic target's constitutive Gaussians in their current state
  std::list<ae_gaussian> current_gaussians_;

  // ----------------------------------------------------------------- Sampling
  /// Number of points to be generated from the gaussians.
  size_t sampling_;

  // ---------------------------------------------------------------- Variation
  /// Variation method
  ae_env_var var_method_;
  /// PRNG used for variation
  std::shared_ptr<ae_jumping_mt> var_prng_;
  /// Autoregressive mean variation sigma parameter
  double var_sigma_;
  /// Autoregressive mean variation tau parameter
  size_t var_tau_; // TODO why size_t ?

  // -------------------------------------------------------------------- Noise
  /// Current noise (pure noise that is added to the phenotypic target)
  Fuzzy* cur_noise_;
  /// PRNG used for noise
  std::shared_ptr<ae_jumping_mt> noise_prng_;
  ae_env_noise noise_method_;
  /// Alpha value (variance coefficient)
  double noise_alpha_;
  /// Variance of the noise
  double noise_sigma_;
  /// Probability of variation.
  double noise_prob_;
  /// Log2 of the number of points in the noise fuzzy_set
  size_t noise_sampling_log_; // TODO why size_t ?
};


// ============================================================================
//                           Getters' definitions
// ============================================================================

// ============================================================================
//                           Setters' definitions
// ============================================================================

// ============================================================================
//                          Operators' definitions
// ============================================================================

// ============================================================================
//                       Inline functions' definition
// ============================================================================

} // namespace aevol

#endif // __AE_PHENOTYPIC_TARGET_HANDLER_H__
