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


/*! \class ae_spatial_structure
    \brief
*/
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_spatial_structure.h>
#include <ae_population.h>




//##############################################################################
//                                                                             #
//                         Class ae_spatial_structure                          #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_spatial_structure::ae_spatial_structure( void )
{
  _pop_grid = NULL;
}

ae_spatial_structure::ae_spatial_structure( gzFile* backup_file )
{
  _prng = new ae_jumping_mt( backup_file );
  
  gzread( backup_file, &_grid_width,  sizeof(_grid_width) );
  gzread( backup_file, &_grid_height, sizeof(_grid_height) );
  
  _pop_grid = new ae_grid_cell** [_grid_width];
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    _pop_grid[x] = new ae_grid_cell* [_grid_height];
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      _pop_grid[x][y] = new ae_grid_cell( backup_file );
    }
  }
  
  gzread( backup_file, &_migration_number,           sizeof(_migration_number) );
  gzread( backup_file, &_secretion_diffusion_prop,   sizeof(_secretion_diffusion_prop) );
  gzread( backup_file, &_secretion_degradation_prop, sizeof(_secretion_degradation_prop) );
}

// =================================================================
//                             Destructors
// =================================================================
ae_spatial_structure::~ae_spatial_structure( void )
{
  if ( _pop_grid != NULL )
  { 
    for ( int16_t i = 0 ; i < _grid_width ; i++ )
    {
      for ( int16_t j = 0 ; j < _grid_height ; j++ )
      {
        delete _pop_grid[i][j];
      }
      
      delete [] _pop_grid[i];
    }
    
    delete [] _pop_grid;
  }
  delete _prng;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_spatial_structure::update_secretion_grid ( void )
{
  int16_t cur_x, cur_y; 
  
  double ** new_secretion = new double*[_grid_width];
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    new_secretion[x] = new double[_grid_height];
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      new_secretion[x][y] = _pop_grid[x][y]->get_compound_amount();
    }
  }

  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      // look at the entire neighborhood
      for ( int8_t i = -1 ; i < 2 ; i++ )
      {
        for ( int8_t j = -1 ; j < 2 ; j ++ )
        {
          cur_x = (x + i + _grid_width)  % _grid_width;
          cur_y = (y + j + _grid_height) % _grid_height;
          
          // add the diffusion from the neighboring cells
          new_secretion[x][y] += _pop_grid[cur_x][cur_y]->get_compound_amount() * _secretion_diffusion_prop;
        }
      }
    }
  }
  
  // substract what has diffused from each cell, and calculate the compound degradation
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      _pop_grid[x][y]->set_compound_amount( new_secretion[x][y] - 9 * _pop_grid[x][y]->get_compound_amount()
                                                                    * _secretion_diffusion_prop );
      _pop_grid[x][y]->set_compound_amount(   _pop_grid[x][y]->get_compound_amount()
                                            * (1 - _secretion_degradation_prop) );
    }
  }
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    delete [] new_secretion[x];
  }  
  delete [] new_secretion; 
}

void ae_spatial_structure::do_random_migrations ( void ) 
{
  ae_individual * tmp_swap;

  int16_t old_x; int16_t old_y; int16_t new_x; int16_t new_y;
  for ( int16_t i = 0 ; i < _migration_number ; i++ )
  {
    old_x = (int16_t) (_prng->random() * _grid_width);
    old_y = (int16_t) (_prng->random() * _grid_height);
    new_x = (int16_t) (_prng->random() * _grid_width);
    new_y = (int16_t) (_prng->random() * _grid_height);


    // swap the individuals in these grid cells...
    tmp_swap = _pop_grid[old_x][old_y]->get_individual();
    _pop_grid[old_x][old_y]->set_individual( _pop_grid[new_x][new_y]->get_individual() );
    _pop_grid[new_x][new_y]->set_individual( tmp_swap );
  }
}

void ae_spatial_structure::save( gzFile* backup_file ) const
{
  if ( _prng == NULL )
  {
    printf( "%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  if ( _pop_grid == NULL )
  {
    printf( "%s:%d: error: grid not initialized.\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  _prng->save( backup_file );
  
  gzwrite( backup_file, &_grid_width,   sizeof(_grid_width) );
  gzwrite( backup_file, &_grid_height,  sizeof(_grid_height) );
  
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      _pop_grid[x][y]->save( backup_file );
    }
  }
  
  gzwrite( backup_file, &_migration_number,           sizeof(_migration_number) );
  gzwrite( backup_file, &_secretion_diffusion_prop,   sizeof(_secretion_diffusion_prop) );
  gzwrite( backup_file, &_secretion_degradation_prop, sizeof(_secretion_degradation_prop) );
}

// =================================================================
//                           Protected Methods
// =================================================================

// =================================================================
//                          Non inline accessors
// =================================================================
