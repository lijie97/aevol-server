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


/*! \class ae_selection
    \brief
*/


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
    ae_selection( ae_exp_manager* exp_m, gzFile* backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_selection( void );

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    // -------------------------------------------------------------- Selection
    inline ae_selection_scheme  get_selection_scheme( void ) const;
    inline double               get_selection_pressure( void ) const;
    
    // --------------------------------------------------------------- Transfer
    inline bool get_with_HT( void ) const;
    //~ inline double get_HT_ins_rate( void ) const;
    //~ inline double get_HT_repl_rate( void ) const;
    inline bool get_with_plasmid_HT( void ) const;
  
    // ------------------------------------------------------ Spatial structure
    inline bool                   is_spatially_structured( void ) const;
    inline ae_spatial_structure*  get_spatial_structure( void ) const;
    inline ae_grid_cell*          get_grid_cell( int16_t x, int16_t y ) const;
    inline int16_t                get_grid_width( void ) const;
    inline int16_t                get_grid_height( void ) const;
  
    // -------------------------------------------------------------- Secretion
    inline bool   get_use_secretion( void ) const;
    inline double get_secretion_contrib_to_fitness( void ) const;
    inline double get_secretion_cost( void ) const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    // ----------------------------------------- Pseudo-random number generator
    inline void set_prng( ae_jumping_mt* prng );

    // -------------------------------------------------------------- Selection
    inline void set_selection_scheme( ae_selection_scheme sel_scheme );
    inline void set_selection_pressure( double sel_pressure );
    
    // --------------------------------------------------------------- Transfer
    inline void set_with_HT( bool with_HT );
    inline void set_HT_ins_rate( double HT_ins_rate );
    inline void set_HT_repl_rate( double HT_repl_rate );
    inline void set_with_plasmid_HT( bool with_p_HT );
    inline void set_nb_plasmid_HT( int16_t nb_p_HT );
    inline void set_prob_plasmid_HT( double prob_p_HT );
    inline void set_swap_GUs( bool swap_GUs );
    
    // ------------------------------------------------------ Spatial structure
    inline void set_spatially_structured( bool structured );
    inline void set_grid_size( int16_t grid_width, int16_t grid_height );
    inline void set_migration_number( int32_t migration_number );
    
    // -------------------------------------------------------------- Secretion
    inline void set_use_secretion( bool use_secretion );
    inline void set_secretion_contrib_to_fitness( double secretion_contrib );
    inline void set_secretion_cost( double secretion_cost );
    inline void set_secretion_degradation_prop( double degradation_prop );
    inline void set_secretion_diffusion_prop( double diffusion_prop );

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void step_to_next_generation( void );
    void step_to_next_generation_grid( void );
    void write_setup_file( gzFile* setup_file ) const;
    void write_setup_file( FILE* setup_file ) const;
    void save( gzFile* sp_struct_file ) const;
    void load( gzFile* exp_setup_file, gzFile* sp_struct_file );
    void load( FILE* exp_setup_file, gzFile* sp_struct_file );

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
    void compute_prob_reprod( void );
    void compute_local_prob_reprod( void );
    ae_individual* do_replication( ae_individual* parent, int32_t index, int16_t x = -1, int16_t y = -1 );
    ae_individual* calculate_local_competition ( int16_t x, int16_t y );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;
    
    // ----------------------------------------- Pseudo-random number generator
    ae_jumping_mt* _prng;

    // -------------------------------------------------------------- Selection
    ae_selection_scheme  _selection_scheme;
    double               _selection_pressure;

    // --------------------------- Probability of reproduction of each organism
    double* _prob_reprod;
    double  _prob_reprod_previous_best;
    
    // --------------------------------------------------------------- Transfer
    bool    _with_HT;
    double  _HT_ins_rate;
    double  _HT_repl_rate;
    int16_t _nb_plasmid_HT;   // TODO: explain
    double  _prob_plasmid_HT; // TODO: explain
    bool    _swap_GUs; // Whether plasmid HT is uni- or bidirectional
    
    // ------------------------------------------------------ Spatial structure
    bool                  _spatially_structured;
    ae_spatial_structure* _spatial_structure;
    
    // -------------------------------------------------------------- Secretion
    bool    _use_secretion;
    double  _secretion_contrib_to_fitness;
    double  _secretion_cost;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline ae_selection_scheme ae_selection::get_selection_scheme( void ) const
{
  return _selection_scheme;
}

inline double ae_selection::get_selection_pressure( void ) const
{
  return _selection_pressure;
}

inline bool ae_selection::get_with_HT( void ) const
{
  return _with_HT;
}

//~ inline double ae_selection::get_HT_ins_rate( void ) const
//~ {
  //~ return _HT_ins_rate;
//~ }

//~ inline double ae_selection::get_HT_repl_rate( void ) const
//~ {
  //~ return _HT_repl_rate;
//~ }

inline bool ae_selection::get_with_plasmid_HT( void ) const
{
  return (_prob_plasmid_HT > 0);
}

inline bool ae_selection::is_spatially_structured( void ) const
{
  return _spatially_structured;
}

inline ae_spatial_structure* ae_selection::get_spatial_structure( void ) const
{
  return _spatial_structure;
}

inline ae_grid_cell* ae_selection::get_grid_cell( int16_t x, int16_t y ) const
{
  return _spatial_structure->get_grid_cell( x, y );
}

inline int16_t ae_selection::get_grid_width( void ) const
{
  return _spatial_structure->get_grid_width();
}

inline int16_t ae_selection::get_grid_height( void ) const
{
  return _spatial_structure->get_grid_height();
}

inline bool ae_selection::get_use_secretion( void ) const
{
  return _use_secretion;
}

inline double ae_selection::get_secretion_contrib_to_fitness( void ) const
{
  return _secretion_contrib_to_fitness;
}

inline double ae_selection::get_secretion_cost( void ) const
{
  return _secretion_cost;
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


// --------------------------------------------------------------- Transfer
inline void ae_selection::set_with_HT( bool with_HT )
{
  _with_HT = with_HT;
}

inline void ae_selection::set_HT_ins_rate( double HT_ins_rate )
{
  _HT_ins_rate = HT_ins_rate;
}

inline void ae_selection::set_HT_repl_rate( double HT_repl_rate )
{
  _HT_repl_rate = HT_repl_rate;
}

inline void ae_selection::set_nb_plasmid_HT( int16_t nb_p_HT )
{
  _nb_plasmid_HT = nb_p_HT;
}

inline void ae_selection::set_prob_plasmid_HT( double prob_p_HT )
{
  _prob_plasmid_HT = prob_p_HT;
}

inline void ae_selection::set_swap_GUs( bool swap_GUs )
{
  _swap_GUs = swap_GUs;
}


// ------------------------------------------------------ Spatial structure
inline void ae_selection::set_spatially_structured( bool structured )
{
  _spatially_structured = structured;
  _spatial_structure    = new ae_spatial_structure();
}

inline void ae_selection::set_grid_size( int16_t grid_width, int16_t grid_height )
{
  assert( _spatially_structured );
  _spatial_structure->set_grid_size( grid_width, grid_height );
}

inline void ae_selection::set_migration_number( int32_t migration_number )
{
  assert( _spatially_structured );
  _spatial_structure->set_migration_number( migration_number );
}

// -------------------------------------------------------------- Secretion
inline void ae_selection::set_use_secretion( bool use_secretion )
{
  _use_secretion = use_secretion;
}

inline void ae_selection::set_secretion_contrib_to_fitness( double secretion_contrib )
{
  _secretion_contrib_to_fitness = secretion_contrib;
}

inline void ae_selection::set_secretion_cost( double secretion_cost )
{
  _secretion_cost = secretion_cost;
}
inline void ae_selection::set_secretion_degradation_prop( double degradation_prop )
{
  assert( _spatially_structured );
  _spatial_structure->set_secretion_degradation_prop(degradation_prop);
}
inline void ae_selection::set_secretion_diffusion_prop( double diffusion_prop )
{
  assert( _spatially_structured );
  _spatial_structure->set_secretion_diffusion_prop(diffusion_prop);
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_SELECTION_H__
