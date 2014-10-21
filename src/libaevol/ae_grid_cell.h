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
// ****************************************************************************


#ifndef __AE_GRID_CELL_H__
#define __AE_GRID_CELL_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ae_individual.h>




// =================================================================
//                          Class declarations
// =================================================================






class ae_grid_cell
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_grid_cell( int16_t x, int16_t y, ae_individual* indiv );
    ae_grid_cell( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_grid_cell( void );


    // =================================================================
    //                        Accessors: getters
    // =================================================================
    inline int16_t get_x( void ) const;
    inline int16_t get_y( void ) const;
    inline double get_compound_amount( void ) const;
    inline ae_individual* get_individual( void ) const;
  
    inline double get_secreted_amount( void ) const;
    inline double get_metabolic_fitness( void ) const;
    inline double get_total_fitness( void ) const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    inline void set_compound_amount( double compound_amount );
    inline void set_individual( ae_individual * indiv );

    // =================================================================
    //                            Public Methods
    // =================================================================
    void save( gzFile backup_file ) const;

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_grid_cell( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    //~ ae_grid_cell( const ae_grid_cell &model )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    // Position on the grid
    int16_t _x;
    int16_t _y;

    // Amount of secreted compound currently present in the grid cell 
    double _compound_amount; 
    
    // pointer to the individual in this cell 
    ae_individual * _individual; 

};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline int16_t ae_grid_cell::get_x( void ) const
{
  return _x;
}

inline int16_t ae_grid_cell::get_y( void ) const
{
  return _y;
}

inline double ae_grid_cell::get_compound_amount( void ) const
{
  return _compound_amount;
}

inline ae_individual* ae_grid_cell::get_individual( void ) const
{
  return _individual;
}

inline double ae_grid_cell::get_secreted_amount( void ) const
{
  return _individual->get_fitness_by_feature( SECRETION );
}

inline double ae_grid_cell::get_metabolic_fitness( void ) const
{
  return _individual->get_fitness_by_feature( METABOLISM );
}

inline double ae_grid_cell::get_total_fitness( void ) const
{
  return _individual->get_fitness();
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void ae_grid_cell::set_compound_amount(double compound_amount)
{
  _compound_amount = compound_amount;
}

inline void ae_grid_cell::set_individual( ae_individual * indiv )
{
  _individual = indiv;
  if ( _individual->get_grid_cell() != this )
  {
    _individual->set_grid_cell( this );
  }
}



// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_GRID_CELL_H__
