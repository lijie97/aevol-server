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




// =================================================================
//                              Libraries
// =================================================================
#include <math.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_environment.h>
#include <point.h>
#include <ae_gaussian.h>
#include <ae_list.h>
#include <list>

namespace aevol {



//##############################################################################
//                                                                             #
//                            Class ae_environment                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_environment::ae_environment( void ) :
#ifdef __NO_X
  ae_fuzzy_set()
#elif defined __X11
  ae_fuzzy_set_X11()
#else
#error You must specify a graphic option
#endif
{
  // Environment "shape"
  _initial_gaussians  = NULL;
  _gaussians          = NULL;
  _sampling           = 0;
  _custom_points      = NULL;

  _total_area = 0.0;

  // Environment segmentation
  _nb_segments      = 1;
  _segments         = new ae_env_segment* [1];
  _segments[0]      = new ae_env_segment( X_MIN, X_MAX, METABOLISM );
  _area_by_feature  = new double [NB_FEATURES];

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



ae_environment::ae_environment( const ae_environment &model ) :
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

  if (model._initial_gaussians == NULL) { _initial_gaussians  = NULL; }
  else
    {
      _initial_gaussians = new ae_list<ae_gaussian *>;
      ae_list_node<ae_gaussian *> * node = model._initial_gaussians->get_first();
      ae_gaussian * gaussian = NULL;
      while (node != NULL)
        {
          gaussian = node->get_obj();
          _initial_gaussians->add( new ae_gaussian(*gaussian));
          node = node->get_next();
        }
    }

  std_initial_gaussians = model.std_initial_gaussians;

  if (model._gaussians == NULL) { _gaussians  = NULL; }
  else
    {
      _gaussians = new ae_list<ae_gaussian *>;
      ae_list_node<ae_gaussian *> * node = model._gaussians->get_first();
      ae_gaussian * gaussian = NULL;
      while (node != NULL)
        {
          gaussian = node->get_obj();
          _gaussians->add( new ae_gaussian(*gaussian));
          node = node->get_next();
        }
    }

  std_gaussians = model.std_gaussians;

  if (model._custom_points == NULL) { _custom_points  = NULL; }
  else
    {
      _custom_points = new ae_list<Point *>;
      ae_list_node<Point *> * node = model._custom_points->get_first();
      Point * m = NULL;
      while (node != NULL)
        {
          m = node->get_obj();
          _custom_points->add( new Point(*m));
          node = node->get_next();
        }
    }

  std_custom_points = model.std_custom_points;

  _total_area = model._total_area;

  // Environment segmentation
  _nb_segments      = model._nb_segments;

  if (_nb_segments == 0) { _segments = NULL; }
  else
    {
      _segments = new ae_env_segment* [_nb_segments];
      for(int32_t i = 0; i < _nb_segments; i++)
        {
          _segments[i] = new ae_env_segment( *(model._segments[i]) );
        }
    }

  _area_by_feature  = new double [NB_FEATURES];
  for (int32_t i = 0; i < NB_FEATURES; i++)
    {
      _area_by_feature[i] = model._area_by_feature[i];
    }


  // Variation management
  _var_method = model._var_method;
  if (model._var_prng == NULL) { _var_prng = NULL; }
  else  { _var_prng   = new ae_jumping_mt( *(model._var_prng) ); }
  _var_sigma  = model._var_sigma;
  _var_tau    = model._var_tau;

  // Noise management
  _noise_method       = model._noise_method;
  if (model._cur_noise == NULL) { _cur_noise = NULL; }
  else { _cur_noise = new ae_fuzzy_set( *(model._cur_noise) ); }
  if (model._noise_prng == NULL) { _noise_prng = NULL; }
  else  { _noise_prng   = new ae_jumping_mt( *(model._noise_prng) ); }
  _noise_prob         = model._noise_prob;
  _noise_alpha        = model._noise_alpha;
  _noise_sigma        = model._noise_sigma;
  _noise_sampling_log = model._noise_sampling_log;


};



// =================================================================
//                             Destructors
// =================================================================
ae_environment::~ae_environment( void )
{
  if ( _gaussians != NULL )
  {
    _gaussians->erase( true );
    delete _gaussians;
  }

  if ( _initial_gaussians != NULL )
  {
    _initial_gaussians->erase( true );
    delete _initial_gaussians;
  }

  if ( _custom_points != NULL )
  {
    _custom_points->erase( true );
    delete _custom_points;
  }

  if (_var_prng != NULL)   delete _var_prng;
  if (_noise_prng != NULL) delete _noise_prng;

  if ( _segments != NULL )
  {
    for ( int16_t i = 0 ; i < _nb_segments; i++ )
    {
      delete _segments[i];
    }

    delete [] _segments;
  }

  delete [] _area_by_feature;

  delete _cur_noise;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_environment::save( gzFile backup_file ) const
{
  // ---------------------
  //  Write gaussians
  // ---------------------
  int16_t nb_gaussians = ( _gaussians == NULL ) ? 0 : _gaussians->get_nb_elts();
  gzwrite( backup_file, &nb_gaussians, sizeof(nb_gaussians) );

  if ( _gaussians != NULL )
  {
    ae_list_node<ae_gaussian*>* gaussian_node = _gaussians->get_first();
    ae_gaussian*  gaussian;
    for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
    {
      gaussian = gaussian_node->get_obj();

      gaussian->save( backup_file );

      gaussian_node = gaussian_node->get_next();
    }
  }

  if (not std_gaussians.empty())
    for (ae_gaussian* g: std_gaussians)
      g->save(backup_file);

  // ---------------------
  //  Write sampling
  // ---------------------
  gzwrite( backup_file, &_sampling, sizeof(_sampling) );

  // ---------------------
  //  Write custom points
  // ---------------------
  int16_t nb_custom_points = (_custom_points == NULL) ? 0 : _custom_points->get_nb_elts();
  gzwrite( backup_file, &nb_custom_points, sizeof(nb_custom_points) );

  if ( _custom_points != NULL )
  {
    ae_list_node<Point*>* custom_point_node = _custom_points->get_first();
    Point*  custom_point;
    for ( int16_t i = 0 ; i < nb_custom_points ; i++ )
    {
      custom_point = custom_point_node->get_obj();

      writepoint(*custom_point, backup_file);

      custom_point_node = custom_point_node->get_next();
    }
  }

  if (not std_custom_points.empty())
    for (Point* p: std_custom_points)
      writepoint(*p, backup_file);

  // -------------------------------
  //  Write x-axis segmentation
  // -------------------------------
  gzwrite( backup_file, &_nb_segments, sizeof(_nb_segments) );

  for ( int16_t i = 0 ; i < _nb_segments; i++ )
  {
    _segments[i]->save( backup_file );
  }


  // -----------------------------------
  //  Write environmental variation data
  // -----------------------------------
  int8_t tmp_var_method = _var_method;
  gzwrite( backup_file, &tmp_var_method,  sizeof(tmp_var_method) );

  if ( _var_method != NO_VAR )
  {
    _var_prng->save( backup_file );
    gzwrite( backup_file, &_var_sigma, sizeof(_var_sigma) );
    gzwrite( backup_file, &_var_tau,   sizeof(_var_tau)   );
  }

  // ---------------------
  //  Write noise data
  // ---------------------
  int8_t tmp_noise_method = _noise_method;
  gzwrite( backup_file, &tmp_noise_method, sizeof(tmp_noise_method) );

  if ( _noise_method != NO_NOISE )
  {
    int8_t tmp_save_cur_noise = ( _cur_noise != NULL );
    gzwrite( backup_file, &tmp_save_cur_noise,  sizeof(tmp_save_cur_noise) );
    if ( tmp_save_cur_noise ) _cur_noise->save( backup_file );

    _noise_prng->save( backup_file );
    gzwrite( backup_file, &_noise_alpha,  sizeof(_noise_alpha) );
    gzwrite( backup_file, &_noise_sigma,  sizeof(_noise_sigma) );
    gzwrite( backup_file, &_noise_prob,   sizeof(_noise_prob) );
    gzwrite( backup_file, &_noise_sampling_log, sizeof(_noise_sampling_log) );
  }

  // ---------------------------------------------------------------
  //  If needed, keep a copy of the initial state of the gaussians
  // ---------------------------------------------------------------
  if ( _var_method != NO_VAR || is_noise_allowed() )
  {
    int16_t nb_gaussians = ( _initial_gaussians == NULL ) ? 0 : _initial_gaussians->get_nb_elts();
    gzwrite( backup_file, &nb_gaussians, sizeof(nb_gaussians) );

    if ( _initial_gaussians != NULL )
    {
      ae_list_node<ae_gaussian*>* gaussian_node = _initial_gaussians->get_first();
      ae_gaussian*  gaussian;
      for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
      {
        gaussian = gaussian_node->get_obj();

        gaussian->save( backup_file );

        gaussian_node = gaussian_node->get_next();
      }
    }

    if (not std_initial_gaussians.empty())
      for (ae_gaussian* g: std_initial_gaussians)
        g->save(backup_file);
  }
}

void ae_environment::load( gzFile backup_file )
{
  // ---------------------
  //  Retreive gaussians
  // ---------------------
  int16_t nb_gaussians;
  gzread( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
  if ( nb_gaussians > 0 ) _gaussians = new ae_list<ae_gaussian*>();
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    _gaussians->add( new ae_gaussian( backup_file ) );
  }

  for (size_t i = 0; i < static_cast<size_t>(nb_gaussians); ++i)
    std_gaussians.push_back(new ae_gaussian(backup_file));

  // ------------------------------
  //  Retrieve sampling
  // ------------------------------
  gzread( backup_file, &_sampling, sizeof(_sampling) );


  // -------------------------
  //  Retreive custom points
  // -------------------------
  int16_t nb_custom_points;
  gzread( backup_file, &nb_custom_points, sizeof(nb_custom_points) );
  if ( nb_custom_points > 0 ) _custom_points = new ae_list<Point*>();
  for ( int16_t i = 0 ; i < nb_custom_points ; i++ )
  {
    _custom_points->add(new Point(readpoint(backup_file)));
  }

  for (size_t i = 0; i < static_cast<size_t>(nb_custom_points); ++i)
    std_custom_points.push_back(new Point(readpoint(backup_file)));

  // -------------------------------
  //  Retrieve x-axis segmentation
  // -------------------------------
  // Delete old data
  for ( int16_t i = 0 ; i < _nb_segments ; i++ )
  {
    delete _segments[i];
  }
  delete [] _segments;


  // Replace by data from the backup
  gzread( backup_file, &_nb_segments, sizeof(_nb_segments) );

  _segments = new ae_env_segment* [_nb_segments];

  for ( int16_t i = 0 ; i < _nb_segments; i++ )
  {
    _segments[i] = new ae_env_segment( backup_file );
    //~ _segments[i] = new ae_env_segment(  ae_common::init_params->get_env_axis_segment_boundaries()[i],
                                        //~ ae_common::init_params->get_env_axis_segment_boundaries()[i+1],
                                        //~ ae_common::init_params->get_env_axis_features()[i] );
  }


  // ----------------------------------------
  //  Retrieve environmental variation data
  // ----------------------------------------
  int8_t tmp_var_method;
  gzread( backup_file, &tmp_var_method, sizeof(tmp_var_method) );
  _var_method = (ae_env_var) tmp_var_method;

  if ( _var_method != NO_VAR )
  {
    _var_prng = new ae_jumping_mt( backup_file );
    gzread( backup_file, &_var_sigma, sizeof(_var_sigma) );
    gzread( backup_file, &_var_tau,   sizeof(_var_tau) );
  }


  // ------------------------------------
  //  Retrieve environmental noise data
  // ------------------------------------);
  int8_t tmp_noise_method;
  gzread( backup_file, &tmp_noise_method, sizeof(tmp_noise_method) );
  _noise_method = (ae_env_noise) tmp_noise_method;

  if ( _noise_method != NO_NOISE )
  {
    int8_t tmp_cur_noise_saved;
    gzread( backup_file, &tmp_cur_noise_saved,  sizeof(tmp_cur_noise_saved) );
    if ( tmp_cur_noise_saved )
    {
      _cur_noise  = new ae_fuzzy_set( backup_file );
    }

    _noise_prng = new ae_jumping_mt( backup_file );
    gzread( backup_file, &_noise_alpha, sizeof(_noise_alpha) );
    gzread( backup_file, &_noise_sigma, sizeof(_noise_sigma) );
    gzread( backup_file, &_noise_prob,  sizeof(_noise_prob)  );
    gzread( backup_file, &_noise_sampling_log, sizeof(_noise_sampling_log) );
  }


  // --------------------------------------------------------------------
  //  If needed, retreive the copy of the initial state of the gaussians
  // --------------------------------------------------------------------
  if ( _var_method != NO_VAR || is_noise_allowed() )
  {
    int16_t nb_gaussians;
    gzread( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
    if ( nb_gaussians > 0 ) _initial_gaussians = new ae_list<ae_gaussian*>();
    for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
    {
      _initial_gaussians->add( new ae_gaussian( backup_file ) );
    }
    for (size_t i = 0 ; i < static_cast<size_t>(nb_gaussians) ; ++i)
      std_initial_gaussians.push_back(new ae_gaussian(backup_file));
  }

  // ------------------------------
  //  Build
  // ------------------------------
  build();
}

void ae_environment::add_gaussian( double a, double b, double c )
{
  _gaussians->add( new ae_gaussian( a, b, c ) );
}

void ae_environment::add_gaussian2(double a, double b, double c) {
  std_gaussians.push_back(new ae_gaussian(a, b, c));
}

void ae_environment::add_initial_gaussian( double a, double b, double c )
{
  _initial_gaussians->add( new ae_gaussian( a, b, c ) );
}

void ae_environment::add_initial_gaussian2(double a, double b, double c) {
  std_initial_gaussians.push_back(new ae_gaussian(a, b, c));
}

void ae_environment::add_custom_point( double x, double y )
{
  _custom_points->add(new Point(x, y));
}

void ae_environment::add_custom_point2(double x, double y) {
  std_custom_points.push_back(new Point(x, y));
}

void ae_environment::build( void )
{
  // NB : Extreme points (at abscissa X_MIN and X_MAX) will be generated, we need to erase the list first
  points.clear();

  // 1) Generate sample points from gaussians
  if ( _gaussians != NULL) {
    ae_list_node<ae_gaussian*>* node = NULL;

    for ( int16_t i = 0 ; i <= _sampling ; i++ ) {
      Point new_point = Point( X_MIN + (double)i * (X_MAX - X_MIN) / (double)_sampling, 0.0 );
      node = _gaussians->get_first();

      while ( node ) {
        new_point.y += node->get_obj()->compute_y( new_point.x );
        node = node->get_next();
      }
      points.push_back(new_point);
    }
  }

  if (not std_gaussians.empty()) {
    for ( int16_t i = 0 ; i <= _sampling ; i++ ) {
      Point new_point = Point(X_MIN + (double)i * (X_MAX - X_MIN) / (double)_sampling, 0.0);
      for (ae_gaussian* g: std_gaussians)
        new_point.y += g->compute_y(new_point.x);
      points.push_back(new_point);
    }
  }

  // 2) Add custom points
  if ( _custom_points != NULL) {
    ae_list_node<Point*>* pt_node = _custom_points->get_first();
    Point *custom_point = pt_node->get_obj();
    Point new_point;

    if ( custom_point->x > X_MIN) {
      // Add the point (X_MIN, Y_MIN) in front of the list of points
      new_point = Point( X_MIN, Y_MIN );
      points.push_front( new_point );
    }

    while ( pt_node != NULL ) {
      custom_point = pt_node->get_obj();
      new_point = Point( *custom_point );
      points.push_back(new_point);
      pt_node = pt_node->get_next();
    }

    if ( custom_point->x < X_MAX ) {
      // Add the point (X_MAX, Y_MIN) at the end of the list of points
      new_point = Point( X_MAX, Y_MIN );
      points.push_back(new_point);
    }
  }

  if (not std_custom_points.empty()) {
    auto p = std_custom_points.begin();

    if ((*p)->x > X_MIN)
      // Add the point (X_MIN, Y_MIN) in front of the list of points
      points.push_front(Point(X_MIN, Y_MIN));

    while (p != std_custom_points.end()) {
      points.push_back(Point(*p));
      ++p;
    }

    if ((*prev(p))->x < X_MAX ) {
      // Add the point (X_MAX, Y_MIN) at the end of the list of points
      points.push_back(Point(X_MAX, Y_MIN));
    }
  }

  // 3) Simplify (get rid of useless points)
  add_lower_bound( Y_MIN );
  add_upper_bound( Y_MAX );
  simplify();

  // 4) Compute areas (total and by feature)
  _compute_area();

  //  5) If needed, create a copy of the initial state of the gaussians
  if ( _initial_gaussians == NULL && (_var_method != NO_VAR || is_noise_allowed()) ) {
    _initial_gaussians = new ae_list<ae_gaussian*>();

    int32_t nb_gaussians = _gaussians->get_nb_elts();
    ae_list_node<ae_gaussian*>* gaussian_node = _gaussians->get_first();
    ae_gaussian*  gaussian      = NULL;
    for ( int16_t i = 0 ; i < nb_gaussians ; i++ ) {
      gaussian = gaussian_node->get_obj();
      _initial_gaussians->add( new ae_gaussian( *gaussian ) );
      gaussian_node = gaussian_node->get_next();
    }
  }

  if (std_initial_gaussians.empty() and (_var_method != NO_VAR or is_noise_allowed()))
    for (ae_gaussian* g: std_gaussians)
      std_initial_gaussians.push_back(new ae_gaussian(*g));
}

/*!
  Apply a gaussian noise to the set of points with a variance sigma and a coefficient alpha :

    alpha = [ -1.0 ------------------ 0.0 -------------------- 1.0 ]
           white noise          uniform fractal             unique draw
                                     noise
 */
void ae_environment::apply_noise( void )
{
  if ( _noise_method != NO_NOISE && _noise_prng->random() < _noise_prob && _noise_sampling_log > 0 )
  {
    // =====================================================================================
    // Compute a fractal noise in a new fuzzy set and apply it to the (unnoised) environment
    // =====================================================================================


    // Clear previous noise (get an unnoised state of th ecurrent environment)
    build();

    // Initialize the cur_noise (current noise) fuzzy set to a set of points with y=0
    // The number of points is determined by _noise_sampling_log (2^_noise_sampling_log)
    if ( _cur_noise == NULL )
    {
      _cur_noise = new ae_fuzzy_set();

      // Add points to reflect the sampling
      int32_t nb_points = 1 << (_noise_sampling_log - 1);
      double interval = (X_MAX - X_MIN) / (nb_points - 1);
      double half_interval = interval / 2;
      double cur_x = interval;
      while ( cur_x - X_MAX < -half_interval ) // while (cur_x < X_MAX) but precision-problems-proof
      {
				printf("infinite loop\n");
        _cur_noise->create_interpolated_point( cur_x );
        cur_x += interval;
      }
    }
    else // _cur_noise has already been created -> reinitialize all its points to 0
    {
      _cur_noise->reset();
    }


    // Compute a fractal noise:
    // Add a random noise to the whole fuzzy set, then cut it in 2 and add
    // another noise to each half and so on (apply noise to the 4 quarters...)
    // until each zone contains only one single point.
    int8_t fractal_step = 0;
    while ( fractal_step < _noise_sampling_log )
    {
      int32_t num_zone;
      int32_t nb_zones = 1 << fractal_step;
      int32_t nb_points_in_each_zone = _cur_noise->get_points().size() / nb_zones;

      // Compute current noise intensity
      // We first test the trivial (most common) cases, then the general (positive or negative) cases
      double noise_intensity;
      if (_noise_alpha == 0.0)
      {
        noise_intensity = _noise_sigma;
      }
      //~ else if (_noise_alpha == 1.0)
      //~ {
        //~ noise_intensity = _noise_sigma * ( 1 << (fractal_step - 1) );
      //~ }
      //~ else if (_noise_alpha == -1.0)
      //~ {
        //~ noise_intensity = _noise_sigma * ( 1 << (_noise_sampling_log - fractal_step - 1) );
      //~ }
      else if ( _noise_alpha > 0.0 )
      {
        noise_intensity = _noise_sigma * pow( 1.0 - _noise_alpha, fractal_step );
      }
      else // ( _noise_alpha < 0.0 )
      {
        noise_intensity = _noise_sigma * pow( 1.0 + _noise_alpha, _noise_sampling_log - fractal_step );
      }
      //~ printf( "_noise_sigma %e\tfractal_step : %"PRId8"\tnoise_intensity : %e\n", _noise_sigma, fractal_step, noise_intensity );

      // For each zone in the current fractal step, compute a random noise to
      // be applied to all the points in the zone
      double* noise_component = new double[nb_zones];
      for ( num_zone = 0 ; num_zone < nb_zones ; num_zone++ )
      {
        // Compute noise component for the current region
        noise_component[num_zone] = _noise_prng->gaussian_random() * noise_intensity;
      }

      // For each point in the noise fuzzy set, apply the noise computed for the corresponding zone
      // Rewritten along with fuzzy sets
      // TODO: test (vld, 2014-12-17)
      {
        std::list<Point> points = _cur_noise->get_points();
        size_t i = 0;
        for (std::list<Point>::iterator p = points.begin() ; p != points.end() ; ++p, ++i) {
          p->y += noise_component[static_cast<size_t>(floor(i / nb_points_in_each_zone))];
        }
      }
      delete noise_component;

      fractal_step++;
    }


    // <DEBUG>
    //~ ae_list_node* point_node  = _cur_noise->get_points()->get_first();
    //~ point*  point       = NULL;
    //~ while ( point_node != NULL )
    //~ {
      //~ point = point_node->get_obj();

      //~ printf( "  x: %f\ty: %e\n", point->x, point->y );

      //~ point_node = point_node->get_next();
    //~ }
    // </DEBUG>


    // Apply the fractal noise to the environment
    this->add( *_cur_noise );

    // Bind Y values in [Y_MIN, Y_MAX]
    add_lower_bound( Y_MIN );
    add_upper_bound( Y_MAX );

    // Environment has changed, recompute its area
    _compute_area();
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
void ae_environment::_apply_autoregressive_mean_variation( void )
{
  // For each gaussian :
  // current_mean = ref_mean + delta_m, where
  // delta_m follows an autoregressive stochastic process
  // with the parameters _var_sigma and _var_tau

  int16_t nb_gaussians = _gaussians->get_nb_elts();

  ae_list_node<ae_gaussian*>* gaussian_node = _gaussians->get_first();
  ae_gaussian* gaussian;
  ae_list_node<ae_gaussian*>* ref_gaussian_node = _initial_gaussians->get_first();
  ae_gaussian* ref_gaussian;
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    gaussian      = gaussian_node->get_obj();
    ref_gaussian  = ref_gaussian_node->get_obj();

    // Find the current delta_mean = current_mean - ref_mean
    double delta_mean = gaussian->get_mean() - ref_gaussian->get_mean();
    //double delta_height = gaussian->get_height() - ref_gaussian->get_height();

    // Compute the next value :
    // Dm(t+1) = Dm(t)*(1-1/tau) + ssd/tau*sqrt(2*tau-1)*normal_random()
    delta_mean =  delta_mean * (1.0 - 1.0/_var_tau) + (_var_sigma/_var_tau) * sqrt(2*_var_tau- 1.0) * _var_prng->gaussian_random();
    //delta_height =  delta_height * (1.0 - 1.0/_var_tau) + (_var_sigma/_var_tau) * sqrt(2*_var_tau- 1.0) * _var_prng->gaussian_random();

    // Deduce the new value of the mean : ref_mean + delta_m
    gaussian->set_mean( ref_gaussian->get_mean() + delta_mean );
    //gaussian->set_height( ref_gaussian->get_height() + delta_height );

    gaussian_node = gaussian_node->get_next();
    ref_gaussian_node = ref_gaussian_node->get_next();
  }

  auto ref = std_initial_gaussians.begin();
  for (ae_gaussian* g: std_gaussians) {
    // Find the current delta_mean = current_mean - ref_mean
    double delta_mean = g->get_mean() - (*ref)->get_mean();

    // Compute the next value:
    // Dm(t+1) = Dm(t) × (1 - 1/tau) + ssd/tau × sqrt(2 tau - 1) × normal_random()
    delta_mean =  delta_mean * (1.0 - 1.0/_var_tau) + (_var_sigma/_var_tau) * sqrt(2*_var_tau- 1.0) * _var_prng->gaussian_random();

    // Deduce the new value of the mean: ref_mean + delta_m
    g->set_mean((*ref)->get_mean() + delta_mean);
    ++ref;
  }


  build();
}


void ae_environment::_apply_autoregressive_height_variation( void )
{
  // For each gaussian :
  // current_height = ref_height + delta_h, where
  // delta_m follows an autoregressive stochastic process
  // with the parameters _var_sigma and _var_tau

  int16_t nb_gaussians = _gaussians->get_nb_elts();

  ae_list_node<ae_gaussian*>* gaussian_node = _gaussians->get_first();
  ae_gaussian* gaussian;
  ae_list_node<ae_gaussian*>* ref_gaussian_node = _initial_gaussians->get_first();
  ae_gaussian* ref_gaussian;
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    gaussian      = gaussian_node->get_obj();
    ref_gaussian  = ref_gaussian_node->get_obj();

    // Find the current delta_height = current_height - ref_height
    double delta_height = gaussian->get_height() - ref_gaussian->get_height();

    // Compute the next value :
    // Dh(t+1) = Dh(t)*(1-1/tau) + ssd/tau*sqrt(2*tau-1)*normal_random()
    delta_height =  delta_height * (1.0 - 1.0/_var_tau) + (_var_sigma/_var_tau) * sqrt(2*_var_tau- 1.0) * _var_prng->gaussian_random();

    // Deduce the new value of the height : ref_height + delta_h
    gaussian->set_height( ref_gaussian->get_height() + delta_height );

    gaussian_node = gaussian_node->get_next();
    ref_gaussian_node = ref_gaussian_node->get_next();
  }

  auto ref = std_initial_gaussians.begin();
  for (ae_gaussian* g: std_gaussians) {
    // Find the current delta_height = current_height - ref_height
    double delta_height = g->get_height() - (*ref)->get_height();

    // Compute the next value :
    // Dh(t+1) = Dh(t) * (1 - 1/tau) + ssd/tau * sqrt(2 tau - 1) * normal_random()
    delta_height =  delta_height * (1.0 - 1.0/_var_tau) + (_var_sigma/_var_tau) * sqrt(2*_var_tau- 1.0) * _var_prng->gaussian_random();

    // Deduce the new value of the height : ref_height + delta_h
    g->set_height((*ref)->get_height() + delta_height );
    ++ref;
  }

  build();
}

void ae_environment::_apply_local_gaussian_variation( void )
{
  printf( "ERROR, _apply_local_gaussian_variation has not yet been implemented. in file %s:%d\n", __FILE__, __LINE__ );
  exit( EXIT_FAILURE );
}

void ae_environment::_compute_area( void )
{
  _total_area = 0.0;

  for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
  {
    _area_by_feature[i] = 0.0;
  }

  // TODO : We should take into account that we compute the areas in a specific order (from the leftmost segment, rightwards)
  //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have
  //      already been through them!)
  for ( int16_t i = 0 ; i < _nb_segments ; i++ )
  {
    _area_by_feature[_segments[i]->feature] += get_geometric_area( _segments[i]->start, _segments[i]->stop );
    _total_area += _area_by_feature[_segments[i]->feature];
  }
}

} // namespace aevol
