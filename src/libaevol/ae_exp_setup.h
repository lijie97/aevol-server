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


/*! \class ae_exp_setup
    \brief This class contains all the data regarding a simulation (either main run or post-treatment)
*/
 
 
#ifndef  __AE_EXP_SETUP_H__
#define  __AE_EXP_SETUP_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_environment.h>
#include <ae_selection.h>
#include <ae_stats.h>
#include <ae_logs.h>
#include <ae_tree.h>
#include <ae_dump.h>
#include <ae_rand_mt.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_param_loader;
class ae_param_overloader;




 
class ae_exp_setup : public ae_object
{
  friend class ae_exp_manager;
  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_exp_setup( ae_exp_manager* exp_m );
    //~ ae_exp_setup( ae_param_overloader* param_overloader /*= NULL*/ );
    //~ ae_exp_setup( char* backup_file_name, bool to_be_run = true, ae_param_overloader* param_overloader = NULL );
    //~ ae_exp_setup( char* organism_file_name, ae_param_overloader* param_overloader = NULL );
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_exp_setup( void );
  
    // =================================================================
    //                         Accessors: getters
    // =================================================================
    inline ae_environment*  get_env( void ) const;
    
    inline bool     fitness_is_composite( void ) const;
    inline int16_t  get_nb_env_segments( void ) const;
    
    inline ae_selection*        get_sel( void ) const;
    inline ae_selection_scheme  get_selection_scheme( void ) const;
    inline double               get_selection_pressure( void ) const;
    
    // Global settings
    //~ inline bool get_with_alignments( void ) const;
    
    inline bool                   is_spatially_structured( void ) const;
    inline ae_spatial_structure*  get_spatial_structure( void ) const;
    inline ae_grid_cell*          get_grid_cell( int16_t x, int16_t y ) const;
    inline int16_t                get_grid_width( void ) const;
    inline int16_t                get_grid_height( void ) const;
    
    inline bool     get_with_plasmid_HT( void ) const;
    
    inline bool   get_use_secretion( void ) const;
    inline double get_secretion_contrib_to_fitness( void ) const;
    inline double get_secretion_cost( void ) const;
  
    // =================================================================
    //                         Accessors: setters
    // =================================================================
    inline void set_env( ae_environment* env );
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    
    // Initialization
    void create_from_param_file( ae_param_loader* param_loader, ae_param_overloader* param_overloader = NULL );
    
    // Backups
    inline void write_setup_file( gzFile* setup_file ) const;
    inline void write_setup_file( FILE* setup_file ) const;
    void save( gzFile* env_file, gzFile* sp_struct_file ) const;
    void load( gzFile* exp_setup_file, gzFile* env_file, gzFile* sp_struct_file, bool verbose );
    void load( FILE* exp_setup_file, gzFile* env_file, gzFile* sp_struct_file, bool verbose );
    
    inline void step_to_next_generation( void );
    
    void run_evolution( void );
    
    void set_total_generations( int32_t nb_generations );
    
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_exp_setup( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_exp_setup( const ae_exp_setup &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    virtual void display( void ) {};
  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;
      
    // Environment
    ae_environment* _env;
      
    // Selection context
    ae_selection*   _sel;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline ae_environment* ae_exp_setup::get_env( void ) const
{
  return _env;
}

inline bool ae_exp_setup::fitness_is_composite( void ) const
{
  return _env->fitness_is_composite();
}

inline int16_t ae_exp_setup::get_nb_env_segments( void ) const
{
  return _env->get_nb_segments();
}

inline ae_selection* ae_exp_setup::get_sel( void ) const
{
  return _sel;
}

inline ae_selection_scheme ae_exp_setup::get_selection_scheme( void ) const
{
  return _sel->get_selection_scheme();
}

inline double ae_exp_setup::get_selection_pressure( void ) const
{
  return _sel->get_selection_pressure();
}

// Global settings    
//~ inline bool ae_exp_setup::get_with_alignments( void ) const
//~ {
  //~ return ;
//~ }

inline bool ae_exp_setup::get_with_plasmid_HT( void ) const
{
  return _sel->get_with_plasmid_HT();
}

inline bool ae_exp_setup::is_spatially_structured( void ) const
{
  return _sel->is_spatially_structured();
}

inline ae_spatial_structure* ae_exp_setup::get_spatial_structure( void ) const
{
  return _sel->get_spatial_structure();
}

inline ae_grid_cell* ae_exp_setup::get_grid_cell( int16_t x, int16_t y ) const
{
  return _sel->get_grid_cell( x, y );
}

inline int16_t ae_exp_setup::get_grid_width( void ) const
{
  return _sel->get_grid_width();
}

inline int16_t ae_exp_setup::get_grid_height( void ) const
{
  return _sel->get_grid_height();
}

inline bool ae_exp_setup::get_use_secretion( void ) const
{
  return _sel->get_use_secretion();
}

inline double ae_exp_setup::get_secretion_contrib_to_fitness( void ) const
{
  return _sel->get_secretion_contrib_to_fitness();
}

inline double ae_exp_setup::get_secretion_cost( void ) const
{
  return _sel->get_secretion_cost();
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void ae_exp_setup::set_env( ae_environment* env )
{
  _env = env;
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
void ae_exp_setup::write_setup_file( gzFile* setup_file ) const
{
  _sel->write_setup_file( setup_file );
}

void ae_exp_setup::write_setup_file( FILE* setup_file ) const
{
  _sel->write_setup_file( setup_file );
}

inline void ae_exp_setup::step_to_next_generation( void )
{
  // Apply environmental variation
  _env->apply_variation();
  
  // Apply environmental noise
  _env->apply_noise();
  
  // Make the individuals reproduce
  _sel->step_to_next_generation();
}


#endif // __AE_EXP_SETUP_H__
