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


#ifndef __AE_GRID_CELL_H__
#define __AE_GRID_CELL_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_individual;






class ae_grid_cell : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    //~ ae_grid_cell( int16_t x, int16_t y );
    ae_grid_cell( int16_t x, int16_t y, ae_individual* indiv );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_grid_cell( void );

    // =================================================================
    //                              Accessors
    // =================================================================

    inline int16_t get_x( void );
    inline int16_t get_y( void );
    inline double get_compound_amount( void );
    inline void set_compound_amount( double compound_amount );
    inline ae_individual* get_individual( void );
    inline void set_individual( ae_individual * indiv );

    // =================================================================
    //                            Public Methods
    // =================================================================

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
    
    // grid size
    int16_t _x;
    int16_t _y;

    // the amount of secreted compound currently present in the grid cell 
    double _compound_amount; 
    
    // pointer to the individual in this cell 
    ae_individual * _individual; 

};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

inline int16_t ae_grid_cell::get_x( void )
{
  return _x;
}

inline int16_t ae_grid_cell::get_y( void )
{
  return _y;
}

inline double ae_grid_cell::get_compound_amount( void )
{
  return _compound_amount;
}

inline void ae_grid_cell::set_compound_amount(double compound_amount)
{
  _compound_amount = compound_amount;
}

inline ae_individual* ae_grid_cell::get_individual( void )
{
  return _individual; 
}

inline void ae_grid_cell::set_individual( ae_individual * indiv )
{
  _individual = indiv;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_GRID_CELL_H__
