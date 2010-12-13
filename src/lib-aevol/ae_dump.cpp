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
#include <ae_simulation.h>
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

#define DUMP_FORMAT "\t%d\t%d\t%lf\n"

void ae_dump :: write_current_generation_dump( void )
{
  //  printf("Begin dump\n");
  write_fitness_total();
  write_secretion_present();
  write_fitness_metabolic();
  write_secreted_amount();
  //  printf("End dump\n");
}

void ae_dump :: write_fitness_total ( void )
{
  if( ae_common :: pop_structure )
  {
    sprintf( filename_buffer, "dump/fitness_total_%04"PRId32".dat", ae_common :: sim -> get_num_gener() ) ;
    current_file = fopen(filename_buffer, "w+");
    ae_population* pop = ae_common :: sim ->get_pop() ;
    double** map = pop->get_fitness_total();
    fprintf( current_file, "#\tX\tY\tfitness_total(X, Y)\n" );
    for( int16_t x = 0 ; x < ae_common::grid_x ; x++ )
    {
      for( int16_t y = 0 ; y < ae_common::grid_y ; y++ )
      {
        fprintf( current_file, DUMP_FORMAT, x, y, map [x][y] );
      }
      fprintf( current_file, "\n" );
    }
    
    fflush( current_file );
    fclose( current_file );
  }
}

void ae_dump :: write_secreted_amount ( void )
{
  if( ae_common :: pop_structure )
  {
    sprintf( filename_buffer, "dump/secreted_amount_%04"PRId32".dat", ae_common :: sim -> get_num_gener() ) ;
    current_file = fopen( filename_buffer, "w+" );
    ae_population* pop = ae_common :: sim ->get_pop() ;
    double** map = pop -> get_secreted_amount();
    fprintf( current_file, "#\tX\tY\tsecreted_amount(X, Y)\n" );
    for( int16_t x = 0 ; x < ae_common::grid_x ; x++ )
    {
      for( int16_t y = 0 ; y < ae_common::grid_y ; y++ )
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
  if( ae_common :: pop_structure )
  {
    sprintf( filename_buffer, "dump/fitness_metabolic_%04"PRId32".dat", ae_common :: sim -> get_num_gener() ) ;
    current_file = fopen( filename_buffer, "w+" );
    ae_population* pop = ae_common :: sim ->get_pop() ;
    double** map = pop -> get_fitness_metabolic();
    fprintf( current_file, "#\tX\tY\tfitness_metabolic(X, Y)\n" );
    for( int16_t x = 0 ; x < ae_common::grid_x ; x++ )
    {
      for( int16_t y = 0 ; y < ae_common::grid_y ; y++ )
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
  if( ae_common :: pop_structure )
  {
    sprintf( filename_buffer, "dump/secretion_present_%04"PRId32".dat", ae_common :: sim -> get_num_gener() ) ;
    current_file = fopen( filename_buffer, "w+" );
    ae_population* pop = ae_common :: sim ->get_pop() ;
    double** map = pop -> get_secretion_present();
    fprintf( current_file, "#\tX\tY\tsecretion_present(X, Y)\n" );
    for( int16_t x = 0 ; x < ae_common::grid_x ; x++ )
    {
      for( int16_t y = 0 ; y < ae_common::grid_y ; y++ )
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
