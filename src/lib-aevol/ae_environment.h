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
#include <ae_rand_mt.h>
#include <ae_fuzzy_set.h>

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
    ae_environment( gzFile* backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_environment( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_list*         get_gaussians( void ) const;
    inline ae_env_segment** get_segments( void ) const;
    inline int16_t          get_nb_segments( void ) const;
    inline double           get_area_by_feature( int feature ) const;
    inline double           get_total_area( void ) const;
    
    inline void   set_sampling( int16_t val );
    inline void   set_variation_method( ae_env_var var_method );
    inline void   set_sigma( double sigma );
    inline void   set_tau( int32_t tau );
    inline void   set_sigma_tau( double sigma, int32_t tau );

    // =================================================================
    //                            Public Methods
    // =================================================================
    void add_custom_point( double x, double y );
    void add_gaussian( double a, double b, double c );
    void build( void );
    inline void apply_variation( void );

    void write_to_backup( gzFile* backup_file );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_environment( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    ae_environment( const ae_environment &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
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
    ae_list*  _gaussians;
    int16_t   _sampling;        // Number of points to be generated from the gaussians.
    ae_list*  _custom_points;
    
    ae_env_segment** _segments; // When the environment is segmented, this is the (ordered) table of segments.
                                // Each ae_env_segment knows its boundaries and corresponding feature.
                                // NULL when the environment is not segmented
    int16_t          _nb_segments;

    double  _total_area;      // Geometric area of the whole function
    double* _area_by_feature; // Geometric area of each feature
    
    // Variation management
    ae_rand_mt* _alea;              // An environment has its own random generator
    ae_env_var  _variation_method;  
    double      _sigma;             // Autoregressive mean variation sigma parameter
    int32_t     _tau;               // Autoregressive mean variation tau parameter
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

inline ae_list* ae_environment::get_gaussians( void ) const
{
  return _gaussians;
}



inline ae_env_segment** ae_environment::get_segments( void ) const
{
  return _segments;
}

inline int16_t ae_environment::get_nb_segments( void ) const
{
  return _nb_segments;
}

inline double ae_environment::get_area_by_feature( int feature ) const
{
  return _area_by_feature[ feature ];
}

inline double ae_environment::get_total_area( void ) const
{
  return _total_area;
}


inline void ae_environment::set_sampling( int16_t val )
{
  _sampling = val;
}

inline void ae_environment::set_variation_method( ae_env_var var_method )
{
  _variation_method = var_method;
}

inline void ae_environment::set_sigma( double sigma )
{
  _sigma = sigma;
}

inline void ae_environment::set_tau( int32_t tau )
{
  _tau = tau;
}

inline void ae_environment::set_sigma_tau( double sigma, int32_t tau )
{
  _sigma = sigma;
  _tau = tau;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_environment::apply_variation( void )
{
  switch ( _variation_method )
  {
    case NONE :
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
  
  _compute_area();
}

#endif // __AE_ENVIRONMENT_H__
