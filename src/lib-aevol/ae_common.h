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

#ifndef __AE_COMMON_H__
#define  __AE_COMMON_H__

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
#include <ae_enums.h>
#include <ae_macros.h>
#include <ae_rand_mt.h>
#include <ae_params_init.h>
#include <ae_params_record.h>
#include <ae_params.h>
#ifdef __REGUL
  #include <ae_array_short.h>
  #include <ae_matrix_double.h>
#endif

// =================================================================
//                          Class declarations
// =================================================================
class ae_simulation;
#ifdef __X11
  class ae_simulation_X11;
#endif
class ae_list;






class ae_common
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================

    // =================================================================
    //                             Destructors
    // =================================================================

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    // Backup functions
    static void write_to_backup( gzFile* backup_file );
    static void read_from_backup( gzFile* backup_file, bool verbose = false );
    static void print_to_file( void );

    // Binding functions
    static          void    init_binding_matrix( void );
    static  inline  double  get_binding( int16_t i, int16_t j );
    
    static void clean( void ); // Release memory

    // =================================================================
    //                           Public Attributes
    // =================================================================
    static ae_simulation* sim;

#ifdef __X11
    static ae_simulation_X11* sim_display;
#endif

    // Initialization parameters
    static ae_params_init* init_params;
    
    // Permanent use parameters
    static ae_params* params;
    
    // Recording parameters
    static ae_params_record* rec_params;
    
    // Number of generations to be simulated
    static int32_t  nb_generations;

    // Population properties
    static bool     pop_structure;
    static int16_t  grid_x; 
    static int16_t  grid_y;
    static int32_t  migration_number;
    

#ifdef __REGUL
    // Regulation factors
    static  double  hill_shape_n;
    static  double  hill_shape_theta;
    static  double  hill_shape;

    // Degradation equation
    static  double  degradation_rate;
    static  double  degradation_step;

    // Individual life
    static  int16_t         individual_evaluation_nbr;
    static  ae_array_short* individual_evaluation_dates;

    // Heredity
    static  bool    with_heredity;
    static  double  protein_presence_limit;
#endif

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_common( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_common( const ae_common &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
#ifdef __REGUL
    // Binding matrix
    static  ae_matrix_double* _binding_matrix;
#endif

};

// =====================================================================
//                          Accessors definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
#ifdef __REGUL
double ae_common::get_binding( int16_t i, int16_t j )
{
 return ae_common::_binding_matrix->get_value( i, j );
}
#endif

#endif // __AE_COMMON_H__
