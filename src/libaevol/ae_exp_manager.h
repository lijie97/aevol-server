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


/*! \class ae_exp_manager
    \brief This class allows to manage an experiment.
    
    An experiment manager allows to... manage an experiment.
    It owns a population and an experimental_setup that can be loaded from a pair of aevol binary files (pop and exp_setup)
*/


#ifndef __AE_EXP_MANAGER_H__
#define __AE_EXP_MANAGER_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_rand_mt.h>
#include <ae_exp_setup.h>
#include <ae_output_manager.h>
#include <ae_population.h>




// =================================================================
//                          Class declarations
// =================================================================






class ae_exp_manager : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_exp_manager( void );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_exp_manager( void );

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    inline ae_population*     get_pop( void ) const;
    inline ae_exp_setup*      get_exp_s( void ) const;
    inline ae_output_manager* get_output_m( void ) const;
  
    inline int32_t  get_first_gener( void ) const;
    inline int32_t  get_num_gener( void ) const;
  
    inline bool quit_signal_received( void ) const;
    
    inline ae_environment* get_env( void ) const;
    
    inline int16_t  get_nb_env_segments( void ) const;
    
    inline ae_selection_scheme get_selection_scheme( void ) const;
    inline double get_selection_pressure( void ) const;
  
    // -------------------------------------------------------- Global settings
    // Global settings accessors
    inline bool                   is_spatially_structured( void ) const;
    inline ae_spatial_structure*  get_spatial_structure( void ) const;
    inline ae_grid_cell*          get_grid_cell( int16_t x, int16_t y ) const;
    inline int16_t                get_grid_width( void ) const;
    inline int16_t                get_grid_height( void ) const;
    
    inline bool     get_with_plasmid_HT( void ) const;
    
    inline bool   fitness_is_composite( void ) const;
    inline bool   get_use_secretion( void ) const;
    inline double get_secretion_contrib_to_fitness( void ) const;
    inline double get_secretion_cost( void ) const;
    
    //~ inline bool   get_with_alignments( void ) const;
    
    // Accessors to population stuff
    inline int32_t        get_nb_indivs( void ) const;
    inline ae_individual* get_best_indiv( void ) const;
    inline ae_list*       get_indivs( void ) const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    inline void set_first_gener( int32_t first_gener );
    inline void set_nb_gener( int32_t nb_gener );
    //~ inline void set_min_genome_length( int32_t min_genome_length );
    //~ inline void set_max_genome_length( int32_t max_genome_length );

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void write_setup_files( void );
    void save_experiment( void ) const;
    void load_experiment( char* exp_setup_file_name,
                          char* out_prof_file_name,
                          char* env_file_name,
                          char* pop_file_name,
                          char* sp_struct_file_name,
                          bool verbose );
    void run_evolution( void );
    virtual void display( void ) {};

    // =================================================================
    //                           Public Attributes
    // =================================================================
    





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*ae_exp_manager( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_exp_manager( const ae_exp_manager &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_rand_mt* _alea;
    
    ae_population*      _pop;
    ae_exp_setup*       _exp_s;
    ae_output_manager*  _output_m;
    
    int32_t _num_gener;
    int32_t _first_gener;
    int32_t _last_gener;
    
    bool _quit_signal_received;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline ae_population* ae_exp_manager::get_pop( void ) const
{
  return _pop;
}

inline ae_exp_setup* ae_exp_manager::get_exp_s( void ) const
{
  return _exp_s;
}

inline ae_output_manager* ae_exp_manager::get_output_m( void ) const
{
  return _output_m;
}

inline int32_t ae_exp_manager::get_first_gener( void ) const
{
  return _first_gener;
}

inline int32_t ae_exp_manager::get_num_gener( void ) const
{
  return _num_gener;
}

inline bool ae_exp_manager::quit_signal_received( void ) const
{
  return _quit_signal_received;
}

inline ae_environment* ae_exp_manager::get_env( void ) const
{
  return _exp_s->get_env();
}

inline int16_t ae_exp_manager::get_nb_env_segments( void ) const
{
  return _exp_s->get_nb_env_segments();
}

inline ae_selection_scheme ae_exp_manager::get_selection_scheme( void ) const
{
  return _exp_s->get_selection_scheme();
}

inline double ae_exp_manager::get_selection_pressure( void ) const
{
  return _exp_s->get_selection_pressure();
}

// Global settings
inline bool ae_exp_manager::is_spatially_structured( void ) const
{
  return _exp_s->is_spatially_structured();
}

inline ae_spatial_structure* ae_exp_manager::get_spatial_structure( void ) const
{
  return _exp_s->get_spatial_structure();
}

inline ae_grid_cell* ae_exp_manager::get_grid_cell( int16_t x, int16_t y ) const
{
  return _exp_s->get_grid_cell( x, y );
}

inline int16_t ae_exp_manager::get_grid_width( void ) const
{
  return _exp_s->get_grid_width();
}

inline int16_t ae_exp_manager::get_grid_height( void ) const
{
  return _exp_s->get_grid_height();
}

inline bool ae_exp_manager::get_with_plasmid_HT( void ) const
{
  return _exp_s->get_with_plasmid_HT();
}

inline bool ae_exp_manager::fitness_is_composite( void ) const
{
  return _exp_s->fitness_is_composite();
}

inline bool ae_exp_manager::get_use_secretion( void ) const
{
  return _exp_s->get_use_secretion();
}

inline double ae_exp_manager::get_secretion_contrib_to_fitness( void ) const
{
  return _exp_s->get_secretion_contrib_to_fitness();
}

inline double ae_exp_manager::get_secretion_cost( void ) const
{
  return _exp_s->get_secretion_cost();
}

//~ inline bool ae_exp_manager::get_with_alignments( void ) const
//~ {
  //~ return _exp_s->get_with_alignments();
//~ }


// Accessors to population stuff
inline int32_t ae_exp_manager::get_nb_indivs( void ) const
{
  return _pop->get_nb_indivs();
}

inline ae_individual* ae_exp_manager::get_best_indiv( void ) const
{
  return _pop->get_best();
}

inline ae_list* ae_exp_manager::get_indivs( void ) const
{
  return _pop->get_indivs();
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void ae_exp_manager::set_first_gener( int32_t first_gener )
{
  _first_gener = first_gener;
}

inline void ae_exp_manager::set_nb_gener( int32_t nb_gener )
{
  _last_gener = _first_gener + nb_gener;
}

// Global constraints
//~ inline void ae_exp_manager::set_min_genome_length( int32_t min_genome_length )
//~ {
  //~ _exp_s->set_min_genome_length( min_genome_length );
//~ }

//~ inline void ae_exp_manager::set_max_genome_length( int32_t max_genome_length )
//~ {
  //~ _exp_s->set_max_genome_length( max_genome_length );
//~ }

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_EXP_MANAGER_H__
