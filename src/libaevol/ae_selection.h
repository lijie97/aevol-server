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


#ifndef __AE_SELECTION_H__
#define __AE_SELECTION_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_spatial_structure.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;






class ae_selection : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_selection( ae_exp_manager* exp_m );
    ae_selection( ae_exp_manager* exp_m, gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_selection( void );

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    inline ae_selection_scheme  get_selection_scheme( void ) const;
    inline double               get_selection_pressure( void ) const;
    inline double*              get_prob_reprod(void) const;
    inline ae_jumping_mt*       get_prng(void) const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    // ----------------------------------------- Pseudo-random number generator
    inline void set_prng( ae_jumping_mt* prng );

    // -------------------------------------------------------------- Selection
    inline void set_selection_scheme( ae_selection_scheme sel_scheme );
    inline void set_selection_pressure( double sel_pressure );

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void step_to_next_generation( void );
    void step_to_next_generation_grid( void );
    void write_setup_file( gzFile setup_file ) const;
    void write_setup_file( FILE* setup_file ) const;
    void save( gzFile& backup_file ) const;
    void load( gzFile& exp_setup_file, gzFile& backup_file, bool verbose );
    void load( FILE*&  exp_setup_file, gzFile& backup_file, bool verbose );
    
    ae_individual* do_replication( ae_individual* parent,
                                   int32_t index,
                                   int16_t x = -1,
                                   int16_t y = -1 );
    void compute_prob_reprod( void );
    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_selection( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_selection( const ae_selection &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };


    // =================================================================
    //                           Protected Methods
    // =================================================================
    //void compute_prob_reprod( void );
    void compute_local_prob_reprod( void );
    //ae_individual* do_replication( ae_individual* parent, int32_t index, int16_t x = -1, int16_t y = -1 );
    ae_individual* calculate_local_competition ( int16_t x, int16_t y );

    // =======================================================================
    //                             Protected Attributes
    // =======================================================================
    ae_exp_manager* _exp_m;
    
    // ----------------------------------------- Pseudo-random number generator
    ae_jumping_mt* _prng;

    // -------------------------------------------------------------- Selection
    ae_selection_scheme  _selection_scheme;
    double               _selection_pressure;

    // --------------------------- Probability of reproduction of each organism
    double* _prob_reprod;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline ae_jumping_mt* ae_selection::get_prng(void) const
{
  return _prng;
}

inline ae_selection_scheme ae_selection::get_selection_scheme( void ) const
{
  return _selection_scheme;
}

inline double ae_selection::get_selection_pressure( void ) const
{
  return _selection_pressure;
}

inline double* ae_selection::get_prob_reprod(void) const
{
  if ( _prob_reprod == NULL )
  {
    printf( "ERROR, _prob_reprod has not been computed %s:%d\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  return _prob_reprod;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
// ----------------------------------------- Pseudo-random number generator
inline void ae_selection::set_prng( ae_jumping_mt* prng )
{
  _prng = prng;
}

// -------------------------------------------------------------- Selection
inline void ae_selection::set_selection_scheme( ae_selection_scheme sel_scheme )
{
  _selection_scheme = sel_scheme;
}

inline void ae_selection::set_selection_pressure( double sel_pressure )
{
  _selection_pressure = sel_pressure;
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_SELECTION_H__
