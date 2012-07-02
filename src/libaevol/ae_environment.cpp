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
  _gaussians      = NULL;
  _sampling       = 0;
  _custom_points  = NULL;
  
  _total_area = 0.0;
  
  // Environment segmentation
  _is_segmented     = false;
  _nb_segments      = 0;
  _segments         = NULL;
  _area_by_feature  = NULL; // TODO: replace by an "area" attribute in ae_segment
  
  // Variation management
  _alea_var         = NULL;
  _variation_method = NONE;  
  _var_sigma        = 0.0;
  _var_tau          = 0;
  
  // Noise management
  _alea_noise = NULL;
}

// =================================================================
//                             Destructors
// =================================================================
ae_environment::~ae_environment( void )
{
  _gaussians->erase( DELETE_OBJ );
  delete _gaussians;
  
  _custom_points->erase( DELETE_OBJ );
  delete _custom_points;
  
  delete _alea_var;
  delete _alea_noise;
  
  
  if ( _is_segmented )
  {
    for ( int16_t i = 0 ; i < _nb_segments; i++ )
    {
      delete _segments[i];
    }
    delete [] _segments;
    
    delete [] _area_by_feature;
  }
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
  _variation_method = ae_common::init_params->get_env_var_method();
  _var_sigma        = ae_common::init_params->get_env_var_sigma();
  _var_tau          = ae_common::init_params->get_env_var_tau();
  _alea_var         = new ae_rand_mt( ae_common::init_params->get_env_seed() );
  
  
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

void ae_environment::write_to_backup( gzFile* backup_file ) const
{
  // ---------------------
  //  Write gaussians
  // ---------------------
  int16_t nb_gaussians = ( _gaussians == NULL ) ? 0 : _gaussians->get_nb_elts();
  gzwrite( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
  
  if ( _gaussians != NULL )
  {
    ae_list_node* gaussian_node = _gaussians->get_first();
    ae_gaussian*  gaussian;
    for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
    {
      gaussian = ( ae_gaussian* ) gaussian_node->get_obj();
      
      gaussian->write_to_backup( backup_file );
      
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
    ae_list_node* custom_point_node = _custom_points->get_first();
    ae_point_2d*  custom_point;
    for ( int16_t i = 0 ; i < nb_custom_points ; i++ )
    {
      custom_point = ( ae_point_2d* ) custom_point_node->get_obj();
      
      custom_point->write_to_backup( backup_file );
      
      custom_point_node = custom_point_node->get_next();
    }
  }
  
  // ---------------------
  //  Write variation data
  // ---------------------
  int8_t tmp_variation_method = _variation_method;
  gzwrite( backup_file, &tmp_variation_method,  sizeof(tmp_variation_method) );
  
  if ( _variation_method != NONE )
  {
    _alea_var->write_to_backup( backup_file );
    gzwrite( backup_file, &_var_sigma, sizeof(_var_sigma) );
    gzwrite( backup_file, &_var_tau,   sizeof(_var_tau)   );
  }
  
  // ---------------------
  //  Write noise data
  // ---------------------
}

void ae_environment::read_from_backup( gzFile* backup_file )
{
  // ---------------------
  //  Retreive gaussians
  // ---------------------
  int16_t nb_gaussians;
  gzread( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
  if ( nb_gaussians > 0 ) _gaussians = new ae_list();
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
  if ( nb_custom_points > 0 ) _custom_points = new ae_list();
  for ( int16_t i = 0 ; i < nb_custom_points ; i++ )
  {
    _custom_points->add( new ae_point_2d( backup_file ) );
  }
  
  
  // -------------------------------
  //  Retrieve x-axis segmentation
  // -------------------------------
  int8_t tmp_is_segmented;
  gzread( backup_file, &tmp_is_segmented, sizeof(tmp_is_segmented) );
  _is_segmented = (tmp_is_segmented != 0);
  
  if ( _is_segmented )
  {
    gzread( backup_file, &_nb_segments, sizeof(_nb_segments) );
    
    _segments = new ae_env_segment* [_nb_segments];
    
    for ( int16_t i = 0 ; i < _nb_segments; i++ )
    {
      _segments[i] = new ae_env_segment();
      _segments[i]->read_from_backup( backup_file );
      //~ _segments[i] = new ae_env_segment(  ae_common::init_params->get_env_axis_segment_boundaries()[i], 
                                          //~ ae_common::init_params->get_env_axis_segment_boundaries()[i+1], 
                                          //~ ae_common::init_params->get_env_axis_features()[i] );
    }
  
    _area_by_feature = new double [NB_FEATURES];
  }
  
  
  // ----------------------------------------
  //  Retrieve environmental variation data
  // ----------------------------------------
  _alea_var = new ae_rand_mt( backup_file );
  int8_t tmp_variation_method;
  gzread( backup_file, &tmp_variation_method, sizeof(tmp_variation_method) );
  _variation_method = (ae_env_var) tmp_variation_method;
  
  if ( _variation_method != NONE )
  {
    gzread( backup_file, &_var_sigma, sizeof(_var_sigma) );  
    gzread( backup_file, &_var_tau, sizeof(_var_tau) );
  }
  
  
  // ------------------------------------
  //  Retrieve environmental noise data
  // ------------------------------------
  #warning TODO
  
  
  // ---------------------------------------------------------------
  //  If needed, keep a copy of the initial state of the gaussians
  // ---------------------------------------------------------------
  if ( 1/*NOISE OR VAR*/ )
  #warning TODO
  {
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

void ae_environment::build( void )
{
  // ----------------------------------------
  // 1) Generate sample points from gaussians
  // ----------------------------------------
  // NB : Extreme points (at abscissa MIN_X and MAX_X) will be generated, we need to erase the list first
  _points->erase( DELETE_OBJ );
  
  ae_list_node* node = NULL;
  
  for ( int16_t i = 0 ; i <= _sampling ; i++ )
  {
    ae_point_2d* new_point = new ae_point_2d( X_MIN + (double)i * (X_MAX - X_MIN) / (double)_sampling, 0.0 );
    node = _gaussians->get_first();
    
    while ( node )
    {
      new_point->y += ((ae_gaussian*)node->get_obj())->compute_y( new_point->x );
      
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
}

bool ae_environment::fitness_is_composite( void ) const
{
  if ( ! _is_segmented ) return false;
  
  int features = 0;
  for ( int16_t i = 0 ; i < _nb_segments ; i++ )
  {
    features |= _segments[i]->feature;
  }
  
  return ( (features & METABOLISM) && (features & SECRETION) );
}

void ae_environment::_compute_area( void )
{
  if ( _is_segmented )
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
  else
  {
    _total_area = get_geometric_area();
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
  // with the parameters ssd and tau
  
  int16_t nb_gaussians = _gaussians->get_nb_elts();
  
  ae_list_node* gaussian_node = _gaussians->get_first();
  ae_gaussian*  gaussian;
  ae_list_node* ref_gaussian_node = _initial_gaussians->get_first();
  ae_gaussian*  ref_gaussian;
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    gaussian      = ( ae_gaussian* ) gaussian_node->get_obj();
    ref_gaussian  = ( ae_gaussian* ) ref_gaussian_node->get_obj();
    
    // find the current delta_mean = current_mean - ref_mean
    double delta_mean = gaussian->get_mean() - ref_gaussian->get_mean();

    // compute the next value :
    // Dm(t+1) = Dm(t)*(1-1/tau) + ssd/tau*sqrt(2*tau-1)*normal_random()
    delta_mean =  delta_mean * (1.0 - 1.0/_var_tau) + (_var_sigma/_var_tau) * sqrt(2*_var_tau- 1.0) * _alea_var->gaussian_random();

    // deduce the new value of the mean : ref_mean + delta_m
    gaussian->set_mean( ref_gaussian->get_mean() + delta_mean );
    
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
