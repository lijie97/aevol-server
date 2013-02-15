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


/*! \class ae_params_init
    \brief Singleton. Contains all the parameters needed only at startup e.g. init_genome_size
*/


#ifndef __AE_PARAMS_INIT_H__
#define __AE_PARAMS_INIT_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_enums.h>
#include <ae_macros.h>
#include <ae_list.h>
#include <ae_gaussian.h>
#include <ae_point_2d.h>




// =================================================================
//                          Class declarations
// =================================================================






class ae_params_init : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_params_init( void );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_params_init( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline int32_t  get_seed( void )                    const;
    inline int32_t  get_env_seed( void )                const;
    inline int32_t  get_initial_genome_length( void )   const;
    inline int8_t   get_init_method( void )             const;
    inline int32_t  get_init_pop_size( void )           const;
    inline bool     get_delete_old_stats( void )        const;
    inline double   get_secretion_init( void )          const;
    inline int32_t  get_plasmid_initial_length( void )  const;
    inline int32_t  get_plasmid_initial_gene( void )    const;
    
    inline ae_list<ae_gaussian*>* get_env_gaussians( void )           const;
    inline ae_list<ae_point_2d*>* get_env_custom_points( void )       const;
    inline int16_t  get_env_sampling( void )            const;
    
    inline bool                 get_env_axis_is_segmented( void )       const;
    inline int16_t              get_env_axis_nb_segments( void )        const;
    inline double*              get_env_axis_segment_boundaries( void ) const;
    inline ae_env_axis_feature* get_env_axis_features( void )           const;
    inline bool                 get_env_axis_separate_segments( void )  const;
    
    inline ae_env_var get_env_var_method( void )  const;
    inline double     get_env_var_sigma( void )   const;
    inline int32_t    get_env_var_tau( void )     const;

    #ifdef __REGUL
      // Binding matrix
      inline double get_binding_zeros_percentage( void ) const;
    #endif
    
    inline void set_seed( int32_t seed );
    inline void set_env_seed( int32_t env_seed );
    inline void set_initial_genome_length( int32_t init_g_len );
    inline void set_init_method( int8_t init_meth );
    inline void set_init_pop_size( int32_t init_pop_size );
    inline void set_delete_old_stats( bool delete_old_stats );
    inline void set_secretion_init( double secretion_init );
    inline void set_plasmid_initial_length( int32_t plasmid_init_len );
    inline void set_plasmid_initial_gene( int32_t plasmid_init_gene );
    
    inline void set_env_sampling( int16_t env_sampling );
    
    inline void set_env_axis_is_segmented( bool env_axis_is_segmented );
    inline void set_env_axis_nb_segments( int16_t env_axis_nb_segments );
    inline void set_env_axis_segment_boundaries( double* env_axis_segment_boundaries ); // set_env_axis_nb_segments first
    inline void set_env_axis_features( ae_env_axis_feature* env_axis_features );        // set_env_axis_nb_segments first
    inline void set_env_axis_separate_segments( bool env_axis_separate_segments );
    
    inline void set_env_var_method( ae_env_var env_var_method );
    inline void set_env_var_sigma( double env_var_sigma );
    inline void set_env_var_tau( int32_t env_var_tau );

    #ifdef __REGUL
      // Binding matrix
      inline void set_binding_zeros_percentage( double zeros_percentage );
    #endif
    
    inline void add_env_gaussian( ae_gaussian* new_gaussian );
    inline void add_env_custom_point( ae_point_2d* new_point );

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    // Backup functions
    void save( gzFile backup_file );
    void load( gzFile backup_file, bool verbose = false );
    void print_to_file( FILE* file );
    
    void clean( void );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_params_init( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_params_init( const ae_params_init &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    // PseudoRandom Number Generator
    int32_t _seed;
    int32_t _env_seed; // Seed for the environment's random generator (environmental variation)

    // Initial conditions
    int32_t  _initial_genome_length;
    int8_t   _init_method;
    int32_t  _init_pop_size;
    
    // Whether to delete the existing statistics file (otherwise kept with the suffix ".old")
    bool _delete_old_stats;
    
    // Environment
    ae_list<ae_gaussian*>* _env_gaussians;
    ae_list<ae_point_2d*>* _env_custom_points;
    int16_t  _env_sampling;
    
    // Environment x-axis segmentation
    bool                 _env_axis_is_segmented;        // true if environment is segmented    
    int16_t              _env_axis_nb_segments;
    double*              _env_axis_segment_boundaries;  // (sorted) including MIN_X and MAX_X
    ae_env_axis_feature* _env_axis_features;            // Corresponding features
    bool                 _env_axis_separate_segments;   // Automatically separate segments
    
    // Environment variation
    ae_env_var  _env_var_method;
    double      _env_var_sigma;
    int32_t     _env_var_tau;
    
    // Environment noise // TODO
    

    // Secretion
    // starting configuration of secretion grid; 0, all are 0; 1, point source of secreted compund
    double _secretion_init;
    
    // Plasmids
    int32_t _plasmid_initial_length;
    int32_t _plasmid_initial_gene;

    #ifdef __REGUL
      // Binding matrix
      double _binding_zeros_percentage;
    #endif
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline int32_t ae_params_init::get_seed( void ) const
{
  return _seed;
}

inline int32_t ae_params_init::get_env_seed( void ) const
{
  return _env_seed;
}

inline int32_t ae_params_init::get_initial_genome_length( void ) const
{
  return _initial_genome_length;
}

inline int8_t ae_params_init::get_init_method( void ) const
{
  return _init_method;
}

inline int32_t ae_params_init::get_init_pop_size( void ) const
{
  return _init_pop_size;
}

inline bool ae_params_init::get_delete_old_stats( void ) const
{
  return _delete_old_stats;
}

inline double ae_params_init::get_secretion_init( void ) const
{
  return _secretion_init;
}

inline int32_t ae_params_init::get_plasmid_initial_length( void ) const
{
  return _plasmid_initial_length;
}

inline int32_t ae_params_init::get_plasmid_initial_gene( void ) const
{
  return _plasmid_initial_gene;
}

inline ae_list<ae_gaussian*>* ae_params_init::get_env_gaussians( void ) const
{
  return _env_gaussians;
}

inline ae_list<ae_point_2d*>* ae_params_init::get_env_custom_points( void ) const
{
  return _env_custom_points;
}

inline int16_t ae_params_init::get_env_sampling( void ) const
{
  return _env_sampling;
}
    
inline bool ae_params_init::get_env_axis_is_segmented( void ) const
{
  return _env_axis_is_segmented;
}

inline int16_t ae_params_init::get_env_axis_nb_segments( void ) const
{
  return _env_axis_nb_segments;
}

inline double* ae_params_init::get_env_axis_segment_boundaries( void ) const
{
  return _env_axis_segment_boundaries;
}

inline ae_env_axis_feature* ae_params_init::get_env_axis_features( void ) const
{
  return _env_axis_features;
}

inline bool ae_params_init::get_env_axis_separate_segments( void ) const
{
  return _env_axis_separate_segments;
}

inline ae_env_var ae_params_init::get_env_var_method( void ) const
{
  return _env_var_method;
}

inline double ae_params_init::get_env_var_sigma( void ) const
{
  return _env_var_sigma;
}

inline int32_t ae_params_init::get_env_var_tau( void ) const
{
  return _env_var_tau;
}


inline void ae_params_init::set_seed( int32_t seed )
{
  _seed = seed;
}

inline void ae_params_init::set_env_seed( int32_t env_seed )
{
  _env_seed = env_seed;
}

inline void ae_params_init::set_initial_genome_length( int32_t init_g_len )
{
  _initial_genome_length = init_g_len;
}

inline void ae_params_init::set_init_method( int8_t init_meth )
{
  _init_method = init_meth;
}

inline void ae_params_init::set_init_pop_size( int32_t init_pop_size )
{
  _init_pop_size = init_pop_size;
}

inline void ae_params_init::set_delete_old_stats( bool delete_old_stats )
{
  _delete_old_stats = delete_old_stats;
}

inline void ae_params_init::set_secretion_init( double secretion_init )
{
  _secretion_init = secretion_init;
}

inline void ae_params_init::set_plasmid_initial_length( int32_t plasmid_init_len )
{
  _plasmid_initial_length = plasmid_init_len;
}

inline void ae_params_init::set_plasmid_initial_gene( int32_t plasmid_init_gene )
{
  _plasmid_initial_gene = plasmid_init_gene;
}

inline void ae_params_init::set_env_sampling( int16_t env_sampling )
{
  _env_sampling = env_sampling;
}
    
inline void ae_params_init::set_env_axis_is_segmented( bool env_axis_is_segmented )
{
  _env_axis_is_segmented = env_axis_is_segmented;
}

inline void ae_params_init::set_env_axis_nb_segments( int16_t env_axis_nb_segments )
{
  _env_axis_nb_segments = env_axis_nb_segments;
}

inline void ae_params_init::set_env_axis_segment_boundaries( double* env_axis_segment_boundaries )
{
  _env_axis_segment_boundaries = env_axis_segment_boundaries;
}

inline void ae_params_init::set_env_axis_features( ae_env_axis_feature* env_axis_features )
{
  _env_axis_features = env_axis_features;
}

inline void ae_params_init::set_env_axis_separate_segments( bool env_axis_separate_segments )
{
  _env_axis_separate_segments = env_axis_separate_segments;
}

inline void ae_params_init::set_env_var_method( ae_env_var env_var_method )
{
  _env_var_method = env_var_method;
}

inline void ae_params_init::set_env_var_sigma( double env_var_sigma )
{
  _env_var_sigma = env_var_sigma;
}

inline void ae_params_init::set_env_var_tau( int32_t env_var_tau )
{
  _env_var_tau = env_var_tau;
}



#ifdef __REGUL
// Binding matrix
inline void ae_params_init::set_binding_zeros_percentage( double zeros_percentage )
{
  _binding_zeros_percentage = zeros_percentage;
}
#endif
    
inline void ae_params_init::add_env_gaussian( ae_gaussian* new_gaussian )
{
  _env_gaussians->add( new_gaussian );
}

inline void ae_params_init::add_env_custom_point( ae_point_2d* new_point )
{
  _env_custom_points->add( new_point );
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_PARAMS_INIT_H__
