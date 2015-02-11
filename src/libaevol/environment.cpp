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

#include <math.h>

#include <list>

#include "environment.h"
#include "point.h"
#include "ae_gaussian.h"

namespace aevol {

Environment::Environment() :
#ifdef __NO_X
    ae_fuzzy_set()
#elif defined __X11
    ae_fuzzy_set_X11()
#else
#error You must specify a graphic option
#endif
{
  // Environment "shape"
  _sampling   = 0;
  _total_area = 0.0;

  // Environment segmentation
  _nb_segments     = 1;
  _segments        = new ae_env_segment* [1];
  _segments[0]     = new ae_env_segment(X_MIN, X_MAX, METABOLISM);
  _area_by_feature = new double [NB_FEATURES];

  // Variation management
  _var_prng   = NULL;
  _var_method = NO_VAR;
  _var_sigma  = 0.0;
  _var_tau    = 0;

  // Noise management
  _cur_noise          = NULL;
  _noise_method       = NO_NOISE;
  _noise_prng         = NULL;
  _noise_prob         = 0.0;
  _noise_alpha        = 0.0;
  _noise_sigma        = 0.0;
  _noise_sampling_log = 8;
}

Environment::Environment(const Environment &model) :
#ifdef __NO_X
    ae_fuzzy_set(model)
#elif defined __X11
    ae_fuzzy_set_X11(model)
#else
#error You must specify a graphic option
#endif
{
  // Environment "shape"
  _sampling           = model._sampling;

  initial_gaussians = model.initial_gaussians;
  gaussians = model.gaussians;

  _total_area = model._total_area;

  // Environment segmentation
  _nb_segments      = model._nb_segments;

  if (_nb_segments == 0)
    _segments = NULL;
  else {
    _segments = new ae_env_segment* [_nb_segments];
    for (size_t i = 0; i < _nb_segments; i++)
      _segments[i] = new ae_env_segment(*(model._segments[i]));
  }

  _area_by_feature  = new double [NB_FEATURES];
  for (size_t i = 0; i < NB_FEATURES; i++)
    _area_by_feature[i] = model._area_by_feature[i];

  // Variation management
  _var_method = model._var_method;
  if (model._var_prng == NULL)
    _var_prng = NULL;
  else
    _var_prng = new ae_jumping_mt(*(model._var_prng));

  _var_sigma  = model._var_sigma;
  _var_tau    = model._var_tau;

  // Noise management
  _noise_method       = model._noise_method;
  if (model._cur_noise == NULL)
    _cur_noise = NULL;
  else
    _cur_noise = new ae_fuzzy_set(*(model._cur_noise));
  if (model._noise_prng == NULL)
    _noise_prng = NULL;
  else
    _noise_prng = new ae_jumping_mt(*(model._noise_prng));

  _noise_prob         = model._noise_prob;
  _noise_alpha        = model._noise_alpha;
  _noise_sigma        = model._noise_sigma;
  _noise_sampling_log = model._noise_sampling_log;
};

Environment::~Environment() {
  delete _var_prng;
  _var_prng = NULL;
  delete _noise_prng;
  _noise_prng = NULL;

  if (_segments != NULL) {
    for (size_t i = 0 ; i < _nb_segments; i++)
      delete _segments[i];
    delete [] _segments;
  }

  delete [] _area_by_feature;

  delete _cur_noise;
}

void Environment::save(gzFile backup_file) const {
  // ---------------------
  //  Write gaussians
  // ---------------------
  size_t nb_gaussians = gaussians.size();
  gzwrite(backup_file, &nb_gaussians, sizeof(nb_gaussians));

  if (not gaussians.empty())
    for (const ae_gaussian& g: gaussians)
      g.save(backup_file);

  // ---------------------
  //  Write sampling
  // ---------------------
  gzwrite(backup_file, &_sampling, sizeof(_sampling));

  // -------------------------------
  //  Write x-axis segmentation
  // -------------------------------
  gzwrite(backup_file, &_nb_segments, sizeof(_nb_segments));

  for (size_t i = 0 ; i < _nb_segments; i++) {
    _segments[i]->save(backup_file);
  }


  // -----------------------------------
  //  Write environmental variation data
  // -----------------------------------
  int8_t tmp_var_method = _var_method;
  gzwrite(backup_file, &tmp_var_method,  sizeof(tmp_var_method));

  if (_var_method != NO_VAR) {
    _var_prng->save(backup_file);
    gzwrite(backup_file, &_var_sigma, sizeof(_var_sigma));
    gzwrite(backup_file, &_var_tau,   sizeof(_var_tau));
  }

  // ---------------------
  //  Write noise data
  // ---------------------
  int8_t tmp_noise_method = _noise_method;
  gzwrite(backup_file, &tmp_noise_method, sizeof(tmp_noise_method));

  if (_noise_method != NO_NOISE) {
    int8_t tmp_save_cur_noise = (_cur_noise != NULL);
    gzwrite(backup_file, &tmp_save_cur_noise,  sizeof(tmp_save_cur_noise));
    if (tmp_save_cur_noise) _cur_noise->save(backup_file);

    _noise_prng->save(backup_file);
    gzwrite(backup_file, &_noise_alpha,  sizeof(_noise_alpha));
    gzwrite(backup_file, &_noise_sigma,  sizeof(_noise_sigma));
    gzwrite(backup_file, &_noise_prob,   sizeof(_noise_prob));
    gzwrite(backup_file, &_noise_sampling_log, sizeof(_noise_sampling_log));
  }

  // ---------------------------------------------------------------
  //  If needed, keep a copy of the initial state of the gaussians
  // ---------------------------------------------------------------
  if (_var_method != NO_VAR || is_noise_allowed()) {
    size_t nb_gaussians = initial_gaussians.size();
    gzwrite(backup_file, &nb_gaussians, sizeof(nb_gaussians));

    if (not initial_gaussians.empty())
      for (const ae_gaussian& g: initial_gaussians)
        g.save(backup_file);
  }
}

void Environment::load(gzFile backup_file) {
  // ---------------------
  //  Retreive gaussians
  // ---------------------
  size_t nb_gaussians;
  gzread(backup_file, &nb_gaussians, sizeof(nb_gaussians));

  for (size_t i = 0; i < static_cast<size_t>(nb_gaussians); ++i)
    gaussians.push_back(ae_gaussian(backup_file));

  // ------------------------------
  //  Retrieve sampling
  // ------------------------------
  gzread(backup_file, &_sampling, sizeof(_sampling));

  // -------------------------------
  //  Retrieve x-axis segmentation
  // -------------------------------
  // Delete old data
  for (size_t i = 0 ; i < _nb_segments ; i++)
    delete _segments[i];

  delete [] _segments;

  // Replace by data from the backup
  gzread(backup_file, &_nb_segments, sizeof(_nb_segments));

  _segments = new ae_env_segment* [_nb_segments];

  for (size_t i = 0 ; i < _nb_segments; i++)
    _segments[i] = new ae_env_segment(backup_file);

  // ----------------------------------------
  //  Retrieve environmental variation data
  // ----------------------------------------
  int8_t tmp_var_method;
  gzread(backup_file, &tmp_var_method, sizeof(tmp_var_method));
  _var_method = (ae_env_var) tmp_var_method;

  if (_var_method != NO_VAR) {
    _var_prng = new ae_jumping_mt(backup_file);
    gzread(backup_file, &_var_sigma, sizeof(_var_sigma));
    gzread(backup_file, &_var_tau,   sizeof(_var_tau));
  }

  // ------------------------------------
  //  Retrieve environmental noise data
  // ------------------------------------);
  int8_t tmp_noise_method;
  gzread(backup_file, &tmp_noise_method, sizeof(tmp_noise_method));
  _noise_method = (ae_env_noise) tmp_noise_method;

  if (_noise_method != NO_NOISE) {
    int8_t tmp_cur_noise_saved;
    gzread(backup_file, &tmp_cur_noise_saved,  sizeof(tmp_cur_noise_saved));
    if (tmp_cur_noise_saved) {
      _cur_noise  = new ae_fuzzy_set(backup_file);
    }

    _noise_prng = new ae_jumping_mt(backup_file);
    gzread(backup_file, &_noise_alpha, sizeof(_noise_alpha));
    gzread(backup_file, &_noise_sigma, sizeof(_noise_sigma));
    gzread(backup_file, &_noise_prob,  sizeof(_noise_prob));
    gzread(backup_file, &_noise_sampling_log, sizeof(_noise_sampling_log));
  }

  // --------------------------------------------------------------------
  //  If needed, retreive the copy of the initial state of the gaussians
  // --------------------------------------------------------------------
  if (_var_method != NO_VAR || is_noise_allowed()) {
    size_t nb_gaussians;
    gzread(backup_file, &nb_gaussians, sizeof(nb_gaussians));
    for (size_t i = 0 ; i < static_cast<size_t>(nb_gaussians) ; ++i)
      initial_gaussians.push_back(ae_gaussian(backup_file));
  }

  // ------------------------------
  //  Build
  // ------------------------------
  build();
}

void Environment::add_gaussian(double a, double b, double c) {
  gaussians.push_back(ae_gaussian(a, b, c));
}

void Environment::add_initial_gaussian(double a, double b, double c) {
  initial_gaussians.push_back(ae_gaussian(a, b, c));
}

void Environment::build() {
  // NB : Extreme points (at abscissa X_MIN and X_MAX) will be generated, we need to erase the list first
  points.clear();

  // 1) Generate sample points from gaussians
  if (not gaussians.empty())
    for (size_t i = 0 ; i <= _sampling ; i++) {
      Point new_point = Point(X_MIN + (double)i * (X_MAX - X_MIN) / (double)_sampling, 0.0);
      for (const ae_gaussian& g: gaussians)
        new_point.y += g.compute_y(new_point.x);
      points.push_back(new_point);
    }

  // 2) Add custom points
  // custom points were unused: removed

  // 3) Simplify (get rid of useless points)
  clip(min, Y_MIN);
  clip(max, Y_MAX);
  simplify();

  // 4) Compute areas (total and by feature)
  _compute_area();

  //  5) If needed, create a copy of the initial state of the gaussians
  if (initial_gaussians.empty() and (_var_method != NO_VAR or is_noise_allowed()))
    for (ae_gaussian& g: gaussians)
      initial_gaussians.push_back(ae_gaussian(g));
}

/*!
  Apply a gaussian noise to the set of points with a variance sigma and a coefficient alpha :

  alpha = [ -1.0 ------------------ 0.0 -------------------- 1.0 ]
  white noise          uniform fractal             unique draw
  noise
*/
void Environment::apply_noise() {
  if (_noise_method != NO_NOISE && _noise_prng->random() < _noise_prob && _noise_sampling_log > 0) {
    // =====================================================================================
    // Compute a fractal noise in a new fuzzy set and apply it to the (unnoised) environment
    // =====================================================================================

    // Clear previous noise (get an unnoised state of th ecurrent environment)
    build();

    // Initialize the cur_noise (current noise) fuzzy set to a set of points with y=0
    // The number of points is determined by _noise_sampling_log (2^_noise_sampling_log)
    if (_cur_noise == NULL) {
      _cur_noise = new ae_fuzzy_set();

      // Add points to reflect the sampling
      size_t nb_points = 1 << (_noise_sampling_log - 1);
      double interval = (X_MAX - X_MIN) / (nb_points - 1);
      double half_interval = interval / 2;
      double cur_x = interval;
      while (cur_x - X_MAX < -half_interval) { // while (cur_x < X_MAX) but precision-problems-proof
        printf("infinite loop\n");
        _cur_noise->create_interpolated_point(cur_x);
        cur_x += interval;
      }
    }
    else { // _cur_noise has already been created -> reinitialize all its points to 0
      _cur_noise->reset();
    }

    // Compute a fractal noise:
    // Add a random noise to the whole fuzzy set, then cut it in 2 and add
    // another noise to each half and so on (apply noise to the 4 quarters...)
    // until each zone contains only one single point.
    size_t fractal_step = 0;
    while (fractal_step < _noise_sampling_log) {
      size_t num_zone;
      size_t nb_zones = 1 << fractal_step;
      size_t nb_points_in_each_zone = _cur_noise->get_points().size() / nb_zones;

      // Compute current noise intensity
      // We first test the trivial (most common) cases, then the general (positive or negative) cases
      double noise_intensity;
      if (_noise_alpha == 0.0)
        noise_intensity = _noise_sigma;
      else if (_noise_alpha > 0.0)
        noise_intensity = _noise_sigma * pow(1.0 - _noise_alpha, fractal_step);
      else // (_noise_alpha < 0.0)
        noise_intensity = _noise_sigma * pow(1.0 + _noise_alpha, _noise_sampling_log - fractal_step);

      // For each zone in the current fractal step, compute a random noise to
      // be applied to all the points in the zone
      double* noise_component = new double[nb_zones];
      for (num_zone = 0 ; num_zone < nb_zones ; num_zone++)
        // Compute noise component for the current region
        noise_component[num_zone] = _noise_prng->gaussian_random() * noise_intensity;

      // For each point in the noise fuzzy set, apply the noise computed for the corresponding zone
      // Rewritten along with fuzzy sets
      // TODO: test (vld, 2014-12-17)
      {
        std::list<Point> points = _cur_noise->get_points();
        size_t i = 0;
        for (std::list<Point>::iterator p = points.begin() ; p != points.end() ; ++p, ++i)
          p->y += noise_component[static_cast<size_t>(floor(i / nb_points_in_each_zone))];
      }
      delete noise_component;

      fractal_step++;
    }

    // Apply the fractal noise to the environment
    this->add(*_cur_noise);

    // Bind Y values in [Y_MIN, Y_MAX]
    clip(min, Y_MIN);
    clip(max, Y_MAX);

    // Environment has changed, recompute its area
    _compute_area();
  }
}

void Environment::_apply_autoregressive_mean_variation() {
  // For each gaussian :
  // current_mean = ref_mean + delta_m, where
  // delta_m follows an autoregressive stochastic process
  // with the parameters _var_sigma and _var_tau

  auto ref = initial_gaussians.begin();
  for (ae_gaussian& g: gaussians) {
    // Find the current delta_mean = current_mean - ref_mean
    double delta_mean = g.get_mean() - ref->get_mean();

    // Compute the next value:
    // Dm(t+1) = Dm(t) × (1 - 1/tau) + ssd/tau × sqrt(2 tau - 1) × normal_random()
    delta_mean =  delta_mean * (1.0 - 1.0/_var_tau) + (_var_sigma/_var_tau) * sqrt(2*_var_tau- 1.0) * _var_prng->gaussian_random();

    // Deduce the new value of the mean: ref_mean + delta_m
    g.set_mean(ref->get_mean() + delta_mean);
    ++ref;
  }

  build();
}

void Environment::_apply_autoregressive_height_variation() {
  // For each gaussian :
  // current_height = ref_height + delta_h, where
  // delta_m follows an autoregressive stochastic process
  // with the parameters _var_sigma and _var_tau

  auto ref = initial_gaussians.begin();
  for (ae_gaussian& g: gaussians) {
    // Find the current delta_height = current_height - ref_height
    double delta_height = g.get_height() - ref->get_height();

    // Compute the next value :
    // Dh(t+1) = Dh(t) * (1 - 1/tau) + ssd/tau * sqrt(2 tau - 1) * normal_random()
    delta_height =  delta_height * (1.0 - 1.0/_var_tau) + (_var_sigma/_var_tau) * sqrt(2*_var_tau- 1.0) * _var_prng->gaussian_random();

    // Deduce the new value of the height : ref_height + delta_h
    g.set_height(ref->get_height() + delta_height);
    ++ref;
  }

  build();
}

void Environment::_apply_local_gaussian_variation() {
  printf("ERROR, _apply_local_gaussian_variation has not yet been implemented. in file %s:%d\n", __FILE__, __LINE__);
  exit(EXIT_FAILURE);
}

void Environment::_compute_area() {
  _total_area = 0.0;

  for (size_t i = 0 ; i < NB_FEATURES ; i++)
    _area_by_feature[i] = 0.0;

  // TODO : We should take into account that we compute the areas in a specific order (from the leftmost segment, rightwards)
  //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have
  //      already been through them!)
  for (size_t i = 0 ; i < _nb_segments ; i++) {
    _area_by_feature[_segments[i]->feature] += get_geometric_area(_segments[i]->start, _segments[i]->stop);
    _total_area += _area_by_feature[_segments[i]->feature];
  }
}
} // namespace aevol
