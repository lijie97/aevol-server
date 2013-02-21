//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, Virginie Lefort
//                    David Parsons
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************


/** \class
 *  \brief
 */
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <math.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_environment.h>
#include <ae_point_2d.h>
#include <ae_gaussian.h>
#include <ae_list.h>




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
  _noise_prng         = NULL;
  _noise_prob         = 0.0;
  _noise_alpha        = 0.0;
  _noise_sigma        = 0.0;
  _noise_sampling_log = 8;
}

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
  
  if ( _custom_points != NULL )
  {
    _custom_points->erase( true );
    delete _custom_points;
  }
  
  delete _var_prng;
  delete _noise_prng;
  
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
/*void load_params( void )
{
  // --------------------------------
  // 1) Copy gaussians from ae_common
  // --------------------------------
  _gaussians = new ae_list();
  int16_t nb_gaussians = ae_common::init_params->get_env_gaussians()->get_nb_elts();
  ae_list_node* ref_gaussian_node = ae_common::init_params->get_env_gaussians()->get_first();
  ae_gaussian*  ref_gaussian;
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    ref_gaussian = ( ae_gaussian* ) ref_gaussian_node->get_obj();

    _gaussians->add( new ae_gaussian( *ref_gaussian ) );

    ref_gaussian_node = ref_gaussian_node->get_next();
  }
  
  
  // ------------------------------------
  // 2) Copy custom points from ae_common
  // ------------------------------------
  _custom_points  = new ae_list();
  int16_t nb_points = ae_common::init_params->get_env_custom_points()->get_nb_elts();
  ae_list_node* ref_point_node = ae_common::init_params->get_env_custom_points()->get_first();
  ae_point_2d*  ref_point;
  for ( int16_t i = 0 ; i < nb_points ; i++ )
  {
    ref_point = ( ae_point_2d* ) ref_point_node->get_obj();

    _custom_points->add( new ae_point_2d( *ref_point ) );

    ref_point_node = ref_point_node->get_next();
  }
  
  
  // -----------------------------
  // 3) Manage x-axis segmentation
  // -----------------------------
  // Check parameters consistency
  if ( ( ae_common::init_params->get_env_axis_features() == NULL && ae_common::init_params->get_env_axis_segment_boundaries() != NULL ) ||
       ( ae_common::init_params->get_env_axis_features() != NULL && ae_common::init_params->get_env_axis_segment_boundaries() == NULL ) )
  {
    printf( "ERROR in param file : you must specifies both the boundaries and the features of the x-axis segments\n" );
    exit( EXIT_FAILURE );
  }
  
  // Copy segmentation parameters from ae_common
  _is_segmented = ae_common::init_params->get_env_axis_is_segmented();
  _nb_segments  = ae_common::init_params->get_env_axis_nb_segments();
  
  if ( _is_segmented )
  {
    _segments = new ae_env_segment* [_nb_segments];
    
    for ( int16_t i = 0 ; i < _nb_segments; i++ )
    {
      _segments[i] = new ae_env_segment(  ae_common::init_params->get_env_axis_segment_boundaries()[i], 
                                          ae_common::init_params->get_env_axis_segment_boundaries()[i+1], 
                                          ae_common::init_params->get_env_axis_features()[i] );
    }
  
    _area_by_feature = new double [NB_FEATURES];
  }
  else
  {
    _segments = NULL;
  }
  
  
  // ---------------------------------
  // 4) Manage environmental variation
  // ---------------------------------
  _var_method = ae_common::init_params->get_env_var_method();
  _var_sigma        = ae_common::init_params->get_env_var_sigma();
  _var_tau          = ae_common::init_params->get_env_var_tau();
  _var_prng         = new ae_jumping_mt( ae_common::init_params->get_env_seed() );
  
  
  // -------------------------
  // 5) Set sampling and build
  // -------------------------
  _sampling = ae_common::init_params->get_env_sampling();
  build();
}*/

//~ void ae_environment::add_custom_point( double x, double y )
//~ {
  //~ _custom_points->add( new ae_point_2d( x, y ) );
//~ }

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
    ae_list_node<ae_point_2d*>* custom_point_node = _custom_points->get_first();
    ae_point_2d*  custom_point;
    for ( int16_t i = 0 ; i < nb_custom_points ; i++ )
    {
      custom_point = custom_point_node->get_obj();
      
      custom_point->save( backup_file );
      
      custom_point_node = custom_point_node->get_next();
    }
  }
  
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
  

  // ------------------------------
  //  Retrieve sampling
  // ------------------------------
  gzread( backup_file, &_sampling, sizeof(_sampling) );
  
  
  // -------------------------
  //  Retreive custom points
  // -------------------------
  int16_t nb_custom_points;
  gzread( backup_file, &nb_custom_points, sizeof(nb_custom_points) );
  if ( nb_custom_points > 0 ) _custom_points = new ae_list<ae_point_2d*>();
  for ( int16_t i = 0 ; i < nb_custom_points ; i++ )
  {
    _custom_points->add( new ae_point_2d( backup_file ) );
  }
  
  
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

void ae_environment::add_initial_gaussian( double a, double b, double c )
{
  _initial_gaussians->add( new ae_gaussian( a, b, c ) );
}

void ae_environment::build( void )
{
  // ----------------------------------------
  // 1) Generate sample points from gaussians
  // ----------------------------------------
  // NB : Extreme points (at abscissa MIN_X and MAX_X) will be generated, we need to erase the list first
  _points->erase( true );
  
  ae_list_node<ae_gaussian*>* node = NULL;
  
  for ( int16_t i = 0 ; i <= _sampling ; i++ )
  {
    ae_point_2d* new_point = new ae_point_2d( X_MIN + (double)i * (X_MAX - X_MIN) / (double)_sampling, 0.0 );
    node = _gaussians->get_first();
    
    while ( node )
    {
      new_point->y += node->get_obj()->compute_y( new_point->x );
      
      node = node->get_next();
    }
    
    _points->add( new_point );
  }
  
  
  // --------------------
  // 2) Add custom points
  // --------------------
  #warning Environment custom points not yet implemented
  
  
  
  // ---------------------------------------
  // 3) Simplify (get rid of useless points)
  // ---------------------------------------
  add_lower_bound( Y_MIN );
  add_upper_bound( Y_MAX );
  simplify();
  
  
  // ---------------------------------------
  // 4) Compute areas (total and by feature)
  // ---------------------------------------
  _compute_area();
  
  
  // -------------------------------------------------------------------
  //  5) If needed, create a copy of the initial state of the gaussians
  // -------------------------------------------------------------------
  if ( _initial_gaussians == NULL && (_var_method != NO_VAR || is_noise_allowed()) )
  {
    _initial_gaussians = new ae_list<ae_gaussian*>();
    
    int32_t nb_gaussians = _gaussians->get_nb_elts();
    ae_list_node<ae_gaussian*>* gaussian_node = _gaussians->get_first();
    ae_gaussian*  gaussian      = NULL;
    for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
    {
      gaussian = gaussian_node->get_obj();
      
      _initial_gaussians->add( new ae_gaussian( *gaussian ) );
      
      gaussian_node = gaussian_node->get_next();
    }
  }
}

/*!
  Apply a gaussian noise to the set of points with a variance sigma and a coefficient alpha :
 
    alpha = [ -1.0 ------------------ 0.0 -------------------- 1.0 ]
           white noise          uniform fractal             unique draw
                                     noise
 */
void ae_environment::apply_noise( void )
{
  if ( _noise_method != NO_NOISE && _noise_prng->random() < _noise_prob )
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
        _cur_noise->create_interpolated_point( cur_x );
        cur_x += interval;
      }
    }
    else // _cur_noise has already been created -> reinitialize all its points to 0
    {
      ae_list_node<ae_point_2d*>* point_node  = _cur_noise->get_points()->get_first();
      while ( point_node != NULL )
      {
        point_node->get_obj()->y = 0;
        
        point_node = point_node->get_next();
      }
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
      int32_t nb_points_in_each_zone = _cur_noise->get_points()->get_nb_elts() / nb_zones;
      
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
      ae_list_node<ae_point_2d*>* point_node  = _cur_noise->get_points()->get_first();
      ae_point_2d*  point       = NULL;
      int32_t       point_index = 0;
      while ( point_node != NULL )
      {
        point = point_node->get_obj();
        
        num_zone = floor( point_index++ / nb_points_in_each_zone );
        point->y += noise_component[num_zone];
        
        point_node = point_node->get_next();
      }
      delete noise_component;
      
      fractal_step++;
    }
    
    
    // <DEBUG>
    //~ ae_list_node* point_node  = _cur_noise->get_points()->get_first();
    //~ ae_point_2d*  point       = NULL;
    //~ while ( point_node != NULL )
    //~ {
      //~ point = point_node->get_obj();
      
      //~ printf( "  x: %f\ty: %e\n", point->x, point->y );
      
      //~ point_node = point_node->get_next();
    //~ }
    // </DEBUG>
    
    
    // Apply the fractal noise to the environment
    this->add( _cur_noise );
    
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

