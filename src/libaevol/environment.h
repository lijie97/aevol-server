// Aevol - An in silico experimental evolution platform
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

#ifndef AEVOL_ENVIRONMENT_H
#define AEVOL_ENVIRONMENT_H

#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <list>
#include <memory>

#include <zlib.h>

#include "ae_env_segment.h"
#include "ae_jumping_mt.h"
#include "fuzzy.h"
#include "ae_gaussian.h"

namespace aevol {

class Environment : public Fuzzy
{
 public :
  // =========================================================================
  //                          Constructors & Destructor
  // =========================================================================
  Environment();
  Environment(const Environment& model);
  Environment(gzFile backup_file);
  virtual ~Environment();

  // =========================================================================
  //                                  Getters
  // =========================================================================
  const std::list<ae_gaussian>& get_gaussians() const {
    return gaussians;
  }
  bool gaussians_provided() const {
    return not gaussians.empty();
  }
  size_t get_nb_segments() const {
    return _nb_segments;
  }
  ae_env_segment** get_segments() const {
    return _segments;
  }
  inline double get_segment_boundaries(size_t i) const;
  inline ae_env_axis_feature get_axis_feature(size_t i) const;
  double get_area_by_feature(size_t feature) const {
    return _area_by_feature[feature];
  }
  ae_env_noise get_noise_method() const {
    return _noise_method;
  }
  bool is_noise_allowed() const {
    return (_noise_method != NO_NOISE);
  }

  // =========================================================================
  //                                  Setters
  // =========================================================================
  void set_gaussians(const std::list<ae_gaussian>& gaussians_model) {
    gaussians = gaussians_model;
  }
  void set_sampling(size_t val){
    _sampling = val;
  }
  void set_segmentation(size_t nb_segments,
                        double* boundaries,
                        ae_env_axis_feature*
                        features,
                        bool separate_segments = false);
  void set_var_method(ae_env_var var_method) {
    _var_method = var_method;
  }
  void set_var_prng(std::shared_ptr<ae_jumping_mt> prng) {
    _var_prng = prng;
  }
  void set_var_sigma(double sigma) {
    _var_sigma = sigma;
  }
  void set_var_tau(size_t tau) {
    _var_tau = tau;
  }
  void set_var_sigma_tau(double sigma, size_t tau) {
    _var_sigma  = sigma;
    _var_tau    = tau;
  }
  void set_noise_method(ae_env_noise noise_method) {
    _noise_method = noise_method;
  }
  void set_noise_prng(std::shared_ptr<ae_jumping_mt> prng) {
    _noise_prng = prng;
  }
  void set_noise_sigma(double sigma) {
    _noise_sigma = sigma;
  }
  void set_noise_alpha(double alpha) {
    _noise_alpha = alpha;
  }
  void set_noise_prob(double prob) {
    _noise_prob = prob;
  }
  void set_noise_sampling_log(size_t sampling_log) {
    _noise_sampling_log = sampling_log;
  }

  // =========================================================================
  //                             Public Methods
  // =========================================================================
  void save(gzFile backup_file) const;
  void load(gzFile backup_file);

  void add_gaussian(double a, double b, double c);
  void add_initial_gaussian(double a, double b, double c);
  void build();
  void clear_gaussians() {
    gaussians.clear();
  }
  void clear_initial_gaussians() {
    initial_gaussians.clear();
  }

  void apply_variation();
  void apply_noise();

 protected :
  // =========================================================================
  //                            Protected Methods
  // =========================================================================
  void _apply_autoregressive_mean_variation();
  void _apply_autoregressive_height_variation();
  void _apply_local_gaussian_variation();
  void _compute_area();


  // =========================================================================
  //                               Attributes
  // =========================================================================
  std::list<ae_gaussian> initial_gaussians; // List containing all the gaussians of the environment in their initial state
  std::list<ae_gaussian> gaussians;         // List containing all the gaussians of the environment
  size_t _sampling;                             // Number of points to be generated from the gaussians.

  size_t _nb_segments;
  ae_env_segment** _segments; // Ordered table of segments.
  // Each ae_env_segment knows its boundaries and corresponding feature.
  // When the environment is not segmented, this table contains a single
  // segment with feature METABOLIC and boundaries MIN_X and MAX_X
  double* _area_by_feature; // Geometric area of each feature

  // Variation management (compatible only with gaussians)
  ae_env_var _var_method;   // Variation method
  std::shared_ptr<ae_jumping_mt> _var_prng; // PRNG used for variation
  double _var_sigma;        // Autoregressive mean variation sigma parameter
  size_t _var_tau;          // Autoregressive mean variation tau parameter

  // Noise management
  Fuzzy* _cur_noise;           // Current noise (pure noise that is added to the environment fuzzy set)
  std::shared_ptr<ae_jumping_mt> _noise_prng;  // PRNG used for noise
  ae_env_noise _noise_method;  // Probability of variation.
  double _noise_alpha;         // Alpha value (variance coefficient)
  double _noise_sigma;         // Variance of the noise
  double _noise_prob;          // Probability of variation.
  size_t _noise_sampling_log;  // Log2 of the number of points in the noise fuzzy_set
};

inline double Environment::get_segment_boundaries(size_t i) const {
  assert(i <= _nb_segments);
  if (i == _nb_segments) return _segments[i-1]->stop;
  else return _segments[i]->start;
}

inline ae_env_axis_feature Environment::get_axis_feature(size_t i) const {
  assert(i < _nb_segments);
  return _segments[i]->feature;
}

inline void Environment::set_segmentation(size_t nb_segments,
                                   double* boundaries,
                                   ae_env_axis_feature* features,
                                   bool separate_segments /*= false*/) {
  // Delete the data to be replaced
  for (size_t i = 0 ; i < _nb_segments ; i++)
    delete _segments[i];
  delete _segments;

  // Now replace with the new data
  _nb_segments  = nb_segments;
  _segments     = new ae_env_segment* [_nb_segments];

  for (size_t i = 0 ; i < _nb_segments; i++)
    _segments[i] = new ae_env_segment(boundaries[i], boundaries[i+1], features[i]);

  // TODO : Manage separate_segments
}

inline void Environment::apply_variation() {
  switch (_var_method)
  {
    case NO_VAR :
      return;
    case AUTOREGRESSIVE_MEAN_VAR :
      _apply_autoregressive_mean_variation();
      break;
    case AUTOREGRESSIVE_HEIGHT_VAR :
      _apply_autoregressive_height_variation();
      break;
    case LOCAL_GAUSSIANS_VAR :
      _apply_local_gaussian_variation();
      break;
    default :
      printf("ERROR : in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
  }

  // Environment has changed, recompute its area
  _compute_area();
}
} // namespace aevol
#endif // AEVOL_ENVIRONMENT_H
