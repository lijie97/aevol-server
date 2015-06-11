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
#include "Gaussian.h"
#include "ae_enums.h"
#include "ae_jumping_mt.h"
#include "Fuzzy.h"

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
  PhenotypicTargetHandler(void); //< Default ctor
  PhenotypicTargetHandler(const PhenotypicTargetHandler&); //< Copy ctor
  PhenotypicTargetHandler(PhenotypicTargetHandler&&) = delete; //< Move ctor
  PhenotypicTargetHandler(gzFile backup_file);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~PhenotypicTargetHandler(void); //< Destructor

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  const PhenotypicTarget& phenotypic_target() const {
    return *phenotypic_target_;
  };
  double get_geometric_area() const {
    return phenotypic_target_->get_geometric_area();
  };
  double area_by_feature(int8_t feature) const {
    return phenotypic_target_->area_by_feature(feature);
  }

  // ==========================================================================
  //                                 Setters
  // ==========================================================================
  void set_gaussians(const list<Gaussian>& gaussians) {
    current_gaussians_ = initial_gaussians_ = gaussians;
  }
  void set_sampling(int16_t val){
    sampling_ = val;
  }
  void set_segmentation(int8_t nb_segments,
                        double* boundaries,
                        ae_env_axis_feature* features,
                        bool separate_segments = false) {
    phenotypic_target_->set_segmentation(nb_segments,
        boundaries,
        features,
        separate_segments);
  };
  void set_var_method(ae_env_var var_method) {
    var_method_ = var_method;
  }
  void set_var_prng(std::shared_ptr<ae_jumping_mt> prng) {
    var_prng_ = prng;
  }
  void set_var_sigma(double sigma) {
    var_sigma_ = sigma;
  }
  void set_var_tau(int32_t tau) {
    var_tau_ = tau;
  }
  void set_var_sigma_tau(double sigma, int32_t tau) {
    var_sigma_  = sigma;
    var_tau_    = tau;
  }
  void set_noise_method(ae_env_noise noise_method) {
    noise_method_ = noise_method;
  }
  void set_noise_prng(std::shared_ptr<ae_jumping_mt> prng) {
    noise_prng_ = prng;
  }
  void set_noise_sigma(double sigma) {
    noise_sigma_ = sigma;
  }
  void set_noise_alpha(double alpha) {
    noise_alpha_ = alpha;
  }
  void set_noise_prob(double prob) {
    noise_prob_ = prob;
  }
  void set_noise_sampling_log(int8_t sampling_log) {
    noise_sampling_log_ = sampling_log;
  }

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void build_phenotypic_target();

  void save(gzFile backup_file) const;
  void load(gzFile backup_file);




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
  std::list<Gaussian> initial_gaussians_;
    /// Phenotypic target's constitutive Gaussians in their current state
  std::list<Gaussian> current_gaussians_;

  // ----------------------------------------------------------------- Sampling
  /// Number of points to be generated from the gaussians.
  int16_t sampling_;

  // ---------------------------------------------------------------- Variation
  /// Variation method
  ae_env_var var_method_;
  /// PRNG used for variation
  std::shared_ptr<ae_jumping_mt> var_prng_;
  /// Autoregressive mean variation sigma parameter
  double var_sigma_;
  /// Autoregressive mean variation tau parameter
  int16_t var_tau_;

  // -------------------------------------------------------------------- Noise
  /// Current noise (pure noise that is added to the phenotypic target)
  Fuzzy* cur_noise_ = NULL;
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
  int8_t noise_sampling_log_;
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
