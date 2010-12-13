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
  // Copy gaussians from ae_common
  _gaussians = new ae_list();
  int16_t nb_gaussians = ae_common::env_gaussians.get_nb_elts();
  ae_list_node* ref_gaussian_node = ae_common::env_gaussians.get_first();
  ae_gaussian*  ref_gaussian;
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    ref_gaussian = ( ae_gaussian* ) ref_gaussian_node->get_obj();

    _gaussians->add( new ae_gaussian( *ref_gaussian ) );

    ref_gaussian_node = ref_gaussian_node->get_next();
  }
  
  // Copy custom points from ae_common
  _custom_points  = new ae_list();
  int16_t nb_points = ae_common::env_custom_points.get_nb_elts();
  ae_list_node* ref_point_node = ae_common::env_custom_points.get_first();
  ae_point_2d*  ref_point;
  for ( int16_t i = 0 ; i < nb_points ; i++ )
  {
    ref_point = ( ae_point_2d* ) ref_point_node->get_obj();

    _custom_points->add( new ae_point_2d( *ref_point ) );

    ref_point_node = ref_point_node->get_next();
  }
  
  _alea = new ae_rand_mt( ae_common::env_seed );

  _sampling         = ae_common::env_sampling;
  _variation_method = ae_common::env_var_method;
  _sigma            = ae_common::env_sigma;
  _tau              = ae_common::env_tau;
  
  build();
}


ae_environment::ae_environment( gzFile* backup_file ) :
#ifdef __NO_X
  ae_fuzzy_set()
#elif defined __X11
  ae_fuzzy_set_X11()
#else
#error You must specify a graphic option
#endif
{
  _area_by_feature  = NULL;
  _segments         = NULL;
  
  // Retreive gaussians
  int16_t nb_gaussians;
  gzread( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
  printf( "    nb_gaussians : %d\n", nb_gaussians );
  
  _gaussians = new ae_list();
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    _gaussians->add( new ae_gaussian( backup_file ) );
  }
 
  
  // Retreive custom points
  int16_t nb_custom_points;
  gzread( backup_file, &nb_custom_points, sizeof(nb_custom_points) );
  _custom_points = new ae_list();
  for ( int16_t i = 0 ; i < nb_custom_points ; i++ )
  {
    _custom_points->add( new ae_point_2d( backup_file ) );
  }
    
  // Retreive variation data
  _alea = new ae_rand_mt( backup_file );
  int8_t tmp_variation_method;
  gzread( backup_file, &tmp_variation_method, sizeof(tmp_variation_method)  );
  _variation_method = (ae_env_var) tmp_variation_method;
  gzread( backup_file, &_sigma,            sizeof(_sigma)             );
  gzread( backup_file, &_tau,              sizeof(_tau)               );


  // Retreive miscellaneous data
  gzread( backup_file, &_sampling,        sizeof(_sampling)       );
  
  // Regenerate environment
  build();
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
  
  delete _alea;
  
  
  
  if ( ae_common::env_axis_is_segmented )
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
//~ void ae_environment::add_custom_point( double x, double y )
//~ {
  //~ _custom_points->add( new ae_point_2d( x, y ) );
//~ }

//~ void ae_environment::add_gaussian( double a, double b, double c )
//~ {
  //~ _gaussians->add( new ae_gaussian( a, b, c ) );
//~ }

void ae_environment::build( void )
{
  // TODO : Check consistency (not gaussians and custom points, not mean variation if points, ...)
  //        Do that in the first constructor and put an if statement here
  
  // ----------------------------------------
  // 1) Generate sample points from gaussians
  // ----------------------------------------
  // NB : Extreme points (at abscissa MIN_X and MAX_X) will be generated, we need to erase the list first
  _points->erase( DELETE_OBJ );
  
  ae_list_node* node = NULL;
  
  for ( int16_t i = 0 ; i <= _sampling ; i++ )
  {
    ae_point_2d* new_point = new ae_point_2d( MIN_X + (double)i * (MAX_X - MIN_X) / (double)_sampling, 0.0 );
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
  // TODO
  
  
  
  
  // ---------------------------------------
  // 3) Simplify (get rid of useless points)
  // ---------------------------------------
  add_lower_bound( MIN_Y );
  add_upper_bound( MAX_Y );
  simplify();
  
  
  
  // TODO : do we really want this here? wouldn't it be better to put it directly in the contructor
  // (or at least to put it in a function whose name tells really what it does)
  // -----------------------------
  // 4) Manage x-axis segmentation
  // -----------------------------
  // Check parameters consistency
  if ( ( ae_common::env_axis_features == NULL && ae_common::env_axis_segment_boundaries != NULL ) ||
       ( ae_common::env_axis_features != NULL && ae_common::env_axis_segment_boundaries == NULL ) )
  {
    printf( "ERROR in param file : you must specifies both the boundaries and the features of the x-axis segments\n" );
    exit( EXIT_FAILURE );
  }
  
  if ( ! ae_common::env_axis_is_segmented )
  {
    _segments = NULL;
  }
  else
  {
    _nb_segments = ae_common::env_axis_nb_segments;
    _segments = new ae_env_segment*[_nb_segments];
    
    for ( int16_t i = 0 ; i < _nb_segments; i++ )
    {
      _segments[i] = new ae_env_segment(  ae_common::env_axis_segment_boundaries[i], 
                                          ae_common::env_axis_segment_boundaries[i+1], 
                                          ae_common::env_axis_features[i] );
    }
  
    _area_by_feature = new double [NB_FEATURES];
  }
  
  //~ // <DEBUG>
  //~ printf( "AXIS SEGMENTATION :\n" );
  //~ for ( int16_t i = 0 ; i < _nb_segments ; i++ )
  //~ {
    //~ printf( "[%lf ; %lf] : ", _segments[i]->start, _segments[i]->stop );
    //~ if ( _segments[i]->feature == NEUTRAL ) printf( "NEUTRAL\n" );
    //~ else if ( _segments[i]->feature == METABOLISM ) printf( "METABOLISM\n" );
    //~ else if ( _segments[i]->feature == SECRETION ) printf( "SECRETION\n" );
    //~ else if ( _segments[i]->feature == TRANSFER ) printf( "TRANSFER\n" );
  //~ }
  //~ // </DEBUG>

  
  
  
  // ---------------------------------------
  // 5) Compute areas (total and by feature)
  // ---------------------------------------
  _compute_area();
}

void ae_environment::_compute_area( void )
{
  if ( ae_common::env_axis_is_segmented )
  {
    _total_area = 0.0;
    for ( int8_t i = 0 ; i < NB_FEATURES ; i++ )
    {
      _area_by_feature[i] = 0.0;
    }

    // TODO : We should take into account that we compute the areas in order (from the leftmost segment, rightwards)
    //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have 
    //      already been through them!)
    for ( int16_t i = 0 ; i < ae_common::env_axis_nb_segments ; i++ )
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

void ae_environment::write_to_backup( gzFile* backup_file )
{
  // Write gaussians
  int16_t nb_gaussians = _gaussians->get_nb_elts();
  gzwrite( backup_file, &nb_gaussians, sizeof(nb_gaussians) );
  
  ae_list_node* gaussian_node = _gaussians->get_first();
  ae_gaussian*  gaussian;
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    gaussian = ( ae_gaussian* ) gaussian_node->get_obj();
    
    gaussian->write_to_backup( backup_file );
    
    gaussian_node = gaussian_node->get_next();
  }
  
  // Write custom points
  int16_t nb_custom_points = _custom_points->get_nb_elts();
  gzwrite( backup_file, &nb_custom_points, sizeof(nb_custom_points) );
  
  ae_list_node* custom_point_node = _custom_points->get_first();
  ae_point_2d*  custom_point;
  for ( int16_t i = 0 ; i < nb_custom_points ; i++ )
  {
    custom_point = ( ae_point_2d* ) custom_point_node->get_obj();
    
    custom_point->write_to_backup( backup_file );
    
    custom_point_node = custom_point_node->get_next();
  }
  
  // Write variation data
  _alea->write_to_backup( backup_file );
  int8_t tmp_variation_method = _variation_method;
  gzwrite( backup_file, &tmp_variation_method,  sizeof(tmp_variation_method) );
  gzwrite( backup_file, &_sigma,                sizeof(_sigma)            );
  gzwrite( backup_file, &_tau,                  sizeof(_tau)              );
  
  // Write miscellaneous data
  gzwrite( backup_file, &_sampling,       sizeof(_sampling)       );
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
  ae_list_node* ref_gaussian_node = ae_common::env_gaussians.get_first();
  ae_gaussian*  ref_gaussian;
  for ( int16_t i = 0 ; i < nb_gaussians ; i++ )
  {
    gaussian      = ( ae_gaussian* ) gaussian_node->get_obj();
    ref_gaussian  = ( ae_gaussian* ) ref_gaussian_node->get_obj();
    
    // find the current delta_mean = current_mean - ref_mean
    double delta_mean = gaussian->get_mean() - ref_gaussian->get_mean();

    // compute the next value :
    // Dm(t+1) = Dm(t)*(1-1/tau) + ssd/tau*sqrt(2*tau-1)*normal_random()
    delta_mean =  delta_mean * (1.0 - 1.0/_tau) + (_sigma/_tau) * sqrt(2*_tau- 1.0) * _alea->gaussian_random();

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
