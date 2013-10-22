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
    ae_spatial_structure( gzFile backup_file );

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
    inline void set_prng( ae_jumping_mt* prng );
    inline void set_secretion_degradation_prop( double degradation_prop );
    inline void set_secretion_diffusion_prop( double diffusion_prop );

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void update_secretion_grid( void ); 
    void do_random_migrations ( void );
    void save( gzFile backup_file ) const;

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
    ae_jumping_mt* _prng;
    
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
inline void ae_spatial_structure::set_prng( ae_jumping_mt* prng )
{
  _prng = prng;
}

inline void ae_spatial_structure::set_grid_size( int16_t grid_width, int16_t grid_height )
{
  _grid_width   = grid_width;
  _grid_height  = grid_height;
  
  _pop_grid = new ae_grid_cell** [_grid_width];
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    _pop_grid[x] = new ae_grid_cell* [_grid_height];
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      _pop_grid[x][y] = new ae_grid_cell( x, y, NULL );
    }
  }
}

inline void ae_spatial_structure::set_migration_number( int32_t migration_number )
{
  _migration_number = migration_number;
}

inline void ae_spatial_structure::set_secretion_degradation_prop( double degradation_prop )
{
  _secretion_degradation_prop=degradation_prop;
}
inline void ae_spatial_structure::set_secretion_diffusion_prop( double diffusion_prop )
{
  _secretion_diffusion_prop=diffusion_prop;
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_SPATIAL_STRUCTURE_H__
