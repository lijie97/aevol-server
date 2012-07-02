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
#include <ae_individual.h>




// =================================================================
//                          Class declarations
// =================================================================






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
    void write_to_backup( gzFile* backup_file ) const;
    void read_from_backup( gzFile* backup_file );

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
}



// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_GRID_CELL_H__
