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


/*! \class
    \brief
*/


#ifndef  __AE_ENVIRONMENT_H__
#define  __AE_ENVIRONMENT_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_env_segment.h>
#include <ae_jumping_mt.h>
#include <ae_fuzzy_set.h>
#include <ae_gaussian.h>
#ifdef __X11
  #include <ae_fuzzy_set_X11.h>
#endif




// =================================================================
//                          Class declarations
// =================================================================







#ifdef __NO_X
class ae_environment : public ae_fuzzy_set
#elif defined __X11
class ae_environment : public ae_fuzzy_set_X11
#else
#error You must specify a graphic option
#endif
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_environment( void );
    ae_environment( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_environment( void );

    // =================================================================
    //                         Accessors: getters
    // =================================================================
    inline ae_list<ae_gaussian*>* get_gaussians( void ) const;
    inline double               get_total_area( void ) const;
    inline int16_t              get_nb_segments( void ) const;
    inline ae_env_segment**     get_segments( void ) const;
    inline double               get_segment_boundaries( int16_t i ) const;
    inline ae_env_axis_feature  get_axis_feature( int16_t i ) const;
    inline double               get_area_by_feature( int8_t feature ) const;
    inline ae_env_var           get_var_method( void ) const;
    inline double               get_var_sigma( void )  const;
    inline int32_t              get_var_tau( void )    const;
    inline ae_env_noise         get_noise_method( void ) const;
    inline bool                 is_noise_allowed( void ) const;
    
    // =================================================================
    //                         Accessors: setters
    // =================================================================
    inline void   set_gaussians( ae_list<ae_gaussian*>* gaussians );
    inline void   set_custom_points( ae_list<ae_point_2d*>* custom_points );
    inline void   set_sampling( int16_t val );
    inline void   set_segmentation( int16_t nb_segments, double* boundaries, ae_env_axis_feature* features, bool separate_segments = false );
    inline void   set_var_method( ae_env_var var_method );
    inline void   set_var_prng( ae_jumping_mt* prng );
    inline void   set_var_sigma( double sigma );
    inline void   set_var_tau( int32_t tau );
    inline void   set_var_sigma_tau( double sigma, int32_t tau );
    inline void   set_noise_method( ae_env_noise noise_method );
    inline void   set_noise_prng( ae_jumping_mt* prng );
    inline void   set_noise_sigma( double sigma );
    inline void   set_noise_alpha( double alpha );
    inline void   set_noise_prob( double prob );
    inline void   set_noise_sampling_log( int32_t sampling_log );


    // =================================================================
    //                            Public Methods
    // =================================================================
    void save( gzFile backup_file ) const;
    void load( gzFile backup_file );

    void add_custom_point( double x, double y );
    void add_gaussian( double a, double b, double c );
    void add_initial_gaussian( double a, double b, double c );
    void build( void );
    
    inline void apply_variation( void );
    void apply_noise( void );


    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_environment( void )
    //~ {
      //~ printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    ae_environment( const ae_environment &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    void _apply_autoregressive_mean_variation( void );
    void _apply_local_gaussian_variation( void );
    void _compute_area( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_list<ae_gaussian*>*  _initial_gaussians; // List containing all the gaussians of the environment in their initial state
    ae_list<ae_gaussian*>*  _gaussians;         // List containing all the gaussians of the environment
    int16_t   _sampling;                        // Number of points to be generated from the gaussians.
    ae_list<ae_point_2d*>*  _custom_points;     // List containing all the custom points of the environment.
                                                // This can not be used in conjunction with gaussians.
    
    int16_t           _nb_segments;
    ae_env_segment**  _segments; // Ordered table of segments.
                                 // Each ae_env_segment knows its boundaries and corresponding feature.
                                 // When the environment is not segmented, this table contains a single
                                 // segment with feature METABOLIC and boundaries MIN_X and MAX_X

    double  _total_area;      // Geometric area of the whole function
    double* _area_by_feature; // Geometric area of each feature
    
    // Variation management
    ae_env_var      _var_method;  // Variation method
    ae_jumping_mt*  _var_prng;    // PRNG used for variation
    double          _var_sigma;   // Autoregressive mean variation sigma parameter
    int32_t         _var_tau;     // Autoregressive mean variation tau parameter
    
    // Noise management
    ae_fuzzy_set*   _cur_noise;           // Current noise (pure noise that is added to the environment fuzzy set)
    ae_jumping_mt*  _noise_prng;          // PRNG used for noise
    ae_env_noise    _noise_method;        // Probability of variation.
    double          _noise_alpha;         // Alpha value (variance coefficient)
    double          _noise_sigma;         // Variance of the noise
    double          _noise_prob;          // Probability of variation.
    int32_t         _noise_sampling_log;  // Log2 of the number of points in the noise fuzzy_set
};


// =====================================================================
//                          Getters' definitions
// =====================================================================
inline int16_t ae_environment::get_nb_segments( void ) const
{
  return _nb_segments;
}

inline ae_list<ae_gaussian*>* ae_environment::get_gaussians( void ) const
{
  return _gaussians;
}

inline ae_env_segment** ae_environment::get_segments( void ) const
{
  return _segments;
}

inline double ae_environment::get_segment_boundaries( int16_t i ) const
{
  assert( i <= _nb_segments );
  
  if ( i == _nb_segments ) return _segments[i-1]->stop;
  else return _segments[i]->start;
}

inline ae_env_axis_feature ae_environment::get_axis_feature( int16_t i ) const
{
  assert( i < _nb_segments );
  
  return _segments[i]->feature;
}

inline double ae_environment::get_area_by_feature( int8_t feature ) const
{
  return _area_by_feature[ feature ];
}

inline double ae_environment::get_total_area( void ) const
{
  return _total_area;
}

inline ae_env_var ae_environment::get_var_method( void ) const
{
  return _var_method;
}

inline double ae_environment::get_var_sigma( void ) const
{
  return _var_sigma;
}

inline int32_t ae_environment::get_var_tau( void ) const
{
  return _var_tau;
}

inline ae_env_noise ae_environment::get_noise_method( void ) const
{
  return _noise_method;
}

inline bool ae_environment::is_noise_allowed( void ) const
{
  return ( _noise_method != NO_NOISE );
}


// =====================================================================
//                          Setters' definitions
// =====================================================================
inline void ae_environment::set_gaussians( ae_list<ae_gaussian*>* gaussians )
{
  _gaussians = gaussians;
}

inline void ae_environment::set_custom_points( ae_list<ae_point_2d*>* custom_points )
{
  _custom_points = custom_points;
}

inline void ae_environment::set_sampling( int16_t val )
{
  _sampling = val;
}


inline void ae_environment::set_segmentation( int16_t nb_segments, double* boundaries, ae_env_axis_feature* features, bool separate_segments /*= false*/ )
{
  // Delete the data to be replaced
  for ( int16_t i = 0 ; i < _nb_segments ; i++ )
  {
    delete _segments[i];
  }
  delete _segments;
  
  
  // Now replace with the new data
  _nb_segments  = nb_segments;
  _segments     = new ae_env_segment* [_nb_segments];
  
  for ( int16_t i = 0 ; i < _nb_segments; i++ )
  {
    _segments[i] = new ae_env_segment( boundaries[i], boundaries[i+1], features[i] );
  }
  
  
  // TODO : Manage separate_segments
}

inline void ae_environment::set_var_method( ae_env_var var_method )
{
  _var_method = var_method;
}

inline void ae_environment::set_var_prng( ae_jumping_mt* prng )
{
  _var_prng = prng;
}

inline void ae_environment::set_var_sigma( double sigma )
{
  _var_sigma = sigma;
}

inline void ae_environment::set_var_tau( int32_t tau )
{
  _var_tau = tau;
}

inline void ae_environment::set_var_sigma_tau( double sigma, int32_t tau )
{
  _var_sigma  = sigma;
  _var_tau    = tau;
}

inline void ae_environment::set_noise_method( ae_env_noise noise_method )
{
  _noise_method = noise_method;
}

inline void ae_environment::set_noise_prng( ae_jumping_mt* prng )
{
  _noise_prng = prng;
}

inline void ae_environment::set_noise_prob( double prob )
{
  _noise_prob = prob;
}

inline void ae_environment::set_noise_alpha( double alpha )
{
  _noise_alpha = alpha;
}

inline void ae_environment::set_noise_sigma( double sigma )
{
  _noise_sigma = sigma;
}

inline void ae_environment::set_noise_sampling_log( int32_t sampling_log )
{
  _noise_sampling_log = sampling_log;
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_environment::apply_variation( void )
{
  switch ( _var_method )
  {
    case NO_VAR :
      return;
    case AUTOREGRESSIVE_MEAN_VAR :
      _apply_autoregressive_mean_variation();
      break;
    case LOCAL_GAUSSIANS_VAR :
      _apply_local_gaussian_variation();
      break;
    default :
      printf( "ERROR : in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
  }
  
  // Environment has changed, recompute its area
  _compute_area();
}

#endif // __AE_ENVIRONMENT_H__
