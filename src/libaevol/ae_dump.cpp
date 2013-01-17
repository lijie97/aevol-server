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
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_dump.h>
#include <ae_exp_manager.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_genetic_unit.h>
#ifdef __REGUL
  #include <ae_influence_R.h>
  #include <ae_protein_R.h>
#endif





//##############################################################################
//                                                                             #
//                                Class ae_dump                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

// =================================================================
//                             Destructors
// =================================================================

// =================================================================
//                            Public Methods
// =================================================================

const char* DUMP_FORMAT = "\t%d\t%d\t%f\n";

void ae_dump::write_current_generation_dump( void )
{
  //  printf("Begin dump\n");
  write_fitness_total();
  write_secretion_present();
  write_fitness_metabolic();
  write_secreted_amount();
  //  printf("End dump\n");
}

void ae_dump::write_fitness_total( void )
{
  if ( _exp_m->is_spatially_structured() )
  {
    sprintf( filename_buffer, "dump/fitness_total_%04"PRId32".dat", _exp_m->get_num_gener() ) ;
    current_file = fopen( filename_buffer, "w+" );
    
    double** map = _exp_m->get_spatial_structure()->get_total_fitness_grid();
    fprintf( current_file, "#\tX\tY\tfitness_total(X, Y)\n" );
    
    for( int16_t x = 0 ; x < _exp_m->get_grid_width() ; x++ )
    {
      for( int16_t y = 0 ; y < _exp_m->get_grid_height() ; y++ )
      {
        fprintf( current_file, DUMP_FORMAT, x, y, map [x][y] );
      }
      fprintf( current_file, "\n" );
    }
    
    fflush( current_file );
    fclose( current_file );
    
    // Has been allocated in ae_spatial_structure::get_total_fitness_grid()
    for ( int16_t x = 0 ; x < _exp_m->get_grid_width() ; x++ )
    {
      delete [] map[x];
    }
    delete [] map; 
  }
}

void ae_dump :: write_secreted_amount ( void )
{
  if ( _exp_m->is_spatially_structured() )
  {
    sprintf( filename_buffer, "dump/secreted_amount_%04"PRId32".dat", _exp_m->get_num_gener() ) ;
    current_file = fopen( filename_buffer, "w+" );
    
    double** map = _exp_m->get_spatial_structure()->get_secreted_amount_grid();
    fprintf( current_file, "#\tX\tY\tsecreted_amount(X, Y)\n" );
    for( int16_t x = 0 ; x < _exp_m->get_grid_width() ; x++ )
    {
      for( int16_t y = 0 ; y < _exp_m->get_grid_height() ; y++ )
      {
        fprintf( current_file, DUMP_FORMAT, x, y, map [x][y] );
      }
      fprintf( current_file, "\n" );
    }
    fflush( current_file );
    fclose( current_file );
  }
}


void ae_dump :: write_fitness_metabolic ( void )
{
  if ( _exp_m->is_spatially_structured() )
  {
    sprintf( filename_buffer, "dump/fitness_metabolic_%04"PRId32".dat", _exp_m->get_num_gener() ) ;
    current_file = fopen( filename_buffer, "w+" );
    
    double** map = _exp_m->get_spatial_structure()->get_metabolic_fitness_grid();
    fprintf( current_file, "#\tX\tY\tfitness_metabolic(X, Y)\n" );
    for( int16_t x = 0 ; x < _exp_m->get_grid_width() ; x++ )
    {
      for( int16_t y = 0 ; y < _exp_m->get_grid_height() ; y++ )
      {
        fprintf( current_file, DUMP_FORMAT, x, y, map [x][y] );
      }
      fprintf( current_file, "\n" );
    }
    fflush( current_file );
    fclose( current_file );
  }
}

void ae_dump :: write_secretion_present ( void )
{
  if ( _exp_m->is_spatially_structured() )
  {
    sprintf( filename_buffer, "dump/secretion_present_%04"PRId32".dat", _exp_m->get_num_gener() ) ;
    current_file = fopen( filename_buffer, "w+" );
    
    double** map = _exp_m->get_spatial_structure()->get_secretion_present_grid();
    fprintf( current_file, "#\tX\tY\tsecretion_present(X, Y)\n" );
    for( int16_t x = 0 ; x < _exp_m->get_grid_width() ; x++ )
    {
      for( int16_t y = 0 ; y < _exp_m->get_grid_height() ; y++ )
        fprintf( current_file, DUMP_FORMAT, x, y, map [x][y] );
      fprintf( current_file, "\n" );
    }
    fflush( current_file );
    fclose( current_file );
  }
}


// =================================================================
//                           Protected Methods
// =================================================================
