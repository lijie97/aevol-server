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


#ifndef __AE_SPATIAL_STRUCTURE_H__
#define __AE_SPATIAL_STRUCTURE_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_grid_cell.h>





// =================================================================
//                          Class declarations
// =================================================================






class ae_spatial_structure : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_spatial_structure( void );
    ae_spatial_structure( gzFile* backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_spatial_structure( void );

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    inline int16_t          get_grid_width( void ) const;
    inline int16_t          get_grid_height( void ) const;
    inline int32_t          get_migration_number( void ) const;
    inline ae_grid_cell***  get_pop_grid( void ) const;
    inline ae_grid_cell*    get_grid_cell( int16_t x, int16_t y ) const;
    inline ae_individual*   get_indiv_at( int16_t x, int16_t y ) const;
  
    inline double** get_secretion_present_grid( void ) const;
    inline double** get_secreted_amount_grid( void ) const;
    inline double** get_metabolic_fitness_grid( void ) const;
    inline double** get_total_fitness_grid( void ) const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    inline void set_grid_size( int16_t grid_x, int16_t grid_y );
    inline void set_migration_number( int32_t migration_number );

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void update_secretion_grid( void ); 
    void do_random_migrations ( void );
    void write_to_backup( gzFile* backup_file ) const;

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*ae_spatial_structure( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_spatial_structure( const ae_spatial_structure &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_rand_mt* _alea;
    
    int16_t _grid_width; 
    int16_t _grid_height;
    
    ae_grid_cell*** _pop_grid;
    
    int32_t _migration_number;
    double  _secretion_diffusion_prop;
    double  _secretion_degradation_prop;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline int16_t ae_spatial_structure::get_grid_width( void ) const
{
  return _grid_width;
}

inline int16_t ae_spatial_structure::get_grid_height( void ) const
{
  return _grid_height;
}

inline int32_t ae_spatial_structure::get_migration_number( void ) const
{
  return _migration_number;
}

inline ae_grid_cell*** ae_spatial_structure::get_pop_grid( void ) const
{
  return _pop_grid;
}

inline ae_grid_cell* ae_spatial_structure::get_grid_cell( int16_t x, int16_t y ) const
{
  return _pop_grid[x][y];
}

inline ae_individual* ae_spatial_structure::get_indiv_at( int16_t x, int16_t y ) const
{
  return _pop_grid[x][y]->get_individual();
}

inline double** ae_spatial_structure::get_secretion_present_grid( void ) const
{
  double** ret = new double*[_grid_width];
  
  for ( int16_t x = 0; x < _grid_width ; x++ )
  {
    ret[x] = new double[_grid_height];
    for ( int16_t y = 0; y < _grid_height ; y++ )
    { 
      ret[x][y] = _pop_grid[x][y]->get_compound_amount(); 
    }
  }
  
  return ret;
}

inline double** ae_spatial_structure::get_secreted_amount_grid( void ) const
{
  double** ret = new double*[_grid_width];
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    ret[x] = new double[_grid_height];
    for ( int16_t y = 0; y < _grid_height ; y++ )
    {
      ret[x][y] = _pop_grid[x][y]->get_secreted_amount();
    }
  }
  
  return ret;
}

inline double** ae_spatial_structure::get_metabolic_fitness_grid( void ) const
{
  double** ret = new double*[_grid_width];
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    ret[x] = new double[_grid_height];
    for ( int16_t y = 0; y < _grid_height ; y++ )
    {
      ret[x][y] = _pop_grid[x][y]->get_metabolic_fitness();
    }
  }
  
  return ret;
}

inline double** ae_spatial_structure::get_total_fitness_grid( void ) const
{
  double** ret = new double*[_grid_width];
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    ret[x] = new double[_grid_height];
    for ( int16_t y = 0; y < _grid_height ; y++ )
    {
      ret[x][y] = _pop_grid[x][y]->get_total_fitness();
    }
  }
  
  return ret;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void ae_spatial_structure::set_grid_size( int16_t grid_width, int16_t grid_height )
{
  _grid_width   = grid_width;
  _grid_height  = grid_height;
}

inline void ae_spatial_structure::set_migration_number( int32_t migration_number )
{
  _migration_number = migration_number;
}

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_SPATIAL_STRUCTURE_H__
