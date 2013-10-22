// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
// 
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 
//*****************************************************************************




// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_grid_cell.h>
#include <ae_individual.h>




//##############################################################################
//                                                                             #
//                              Class ae_grid_cell                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_grid_cell::ae_grid_cell( int16_t x, int16_t y, ae_individual* indiv )
{
  _x = x;
  _y = y;
  
  _compound_amount  = 0.0;
  
  _individual       = indiv;
}

ae_grid_cell::ae_grid_cell( gzFile backup_file )
{
  gzread( backup_file, &_x, sizeof(_x) );
  gzread( backup_file, &_y, sizeof(_y) );
  
  gzread( backup_file, &_compound_amount, sizeof(_compound_amount) );
  
  _individual = NULL;
}

// =================================================================
//                             Destructors
// =================================================================
ae_grid_cell::~ae_grid_cell( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_grid_cell::save( gzFile backup_file ) const
{
  gzwrite( backup_file, &_x, sizeof(_x) );
  gzwrite( backup_file, &_y, sizeof(_y) );
  
  gzwrite( backup_file, &_compound_amount, sizeof(_compound_amount) );
}

// =================================================================
//                           Protected Methods
// =================================================================
