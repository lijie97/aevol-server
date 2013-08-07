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


/*!
  \class ae_exp_manager
  
  \brief This is aevol's top-level class. It allows for
  high-level experiment management
  
  An experiment manager allows to... manage an experiment.
  It owns a population and an experimental_setup that can be loaded from a
  pair of aevol binary files (pop and exp_setup)
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
#include <ae_jumping_mt.h>
#include <ae_exp_setup.h>
#include <ae_output_manager.h>
#include <ae_population.h>




// =================================================================
//                          Class declarations
// =================================================================






class ae_exp_manager : public ae_object
{
  public :
    // =======================================================================
    //                                Constructors
    // =======================================================================
    ae_exp_manager( void );

    // =======================================================================
    //                                Destructors
    // =======================================================================
    virtual ~ae_exp_manager( void );

    // =======================================================================
    //                           Accessors: getters
    // =======================================================================
    inline ae_exp_setup*        get_exp_s( void ) const;
    inline ae_population*       get_pop( void ) const;
    inline ae_environment*      get_env( void ) const;
    inline ae_selection*        get_sel( void ) const;
    inline ae_output_manager*   get_output_m( void ) const;
  
    inline int32_t  get_first_gener( void ) const;
    inline int32_t  get_num_gener( void ) const;
  
    inline bool quit_signal_received( void ) const;
    
    
    inline int16_t  get_nb_env_segments( void ) const;
    
    inline ae_selection_scheme get_selection_scheme( void ) const;
    inline double get_selection_pressure( void ) const;
  
    // ------------------------------------------------------ Spatial structure
    inline bool                   is_spatially_structured( void ) const;
    inline ae_spatial_structure*  get_spatial_structure( void ) const;
    inline ae_grid_cell*          get_grid_cell( int16_t x, int16_t y ) const;
    inline int16_t                get_grid_width( void ) const;
    inline int16_t                get_grid_height( void ) const;
    inline ae_grid_cell***        get_pop_grid( void ) const;
  
    // -------------------------------------------------------- Global settings
    inline bool   get_with_HT( void ) const;
    inline double get_HT_ins_rate( void ) const;
    inline double get_HT_repl_rate( void ) const;
    
    // The ability to own a plasmid is a property of the individuals (_allow_plasmids) because it is used during mutations
    // However there is also a property of the experimental setup (_with_plasmids) that indicates whether plasmids are used because we need this during replication and during loading/writting
    // For now when plasmids are used each individual has one and only one plasmid (so these variables should always be equals), however this may change in the future
    // There is no longer property _with_plasmids_HT because the ability to transfer is evolvable and thus may depend on the plasmid itself
  
    inline bool   get_with_plasmids( void ) const;
    inline double get_prob_plasmid_HT( void ) const;
    inline double get_tune_donor_ability( void ) const;
    inline double get_tune_recipient_ability( void ) const;
    inline double get_donor_cost( void ) const;
    inline double get_recipient_cost( void ) const;
    inline bool   get_swap_GUs( void ) const;
    
    inline bool   get_with_secretion( void ) const;
    inline double get_secretion_contrib_to_fitness( void ) const;
    inline double get_secretion_cost( void ) const;
    
    //~ inline bool   get_with_alignments( void ) const;
    
    // Accessors to population stuff
    inline ae_list<ae_individual*>* get_indivs( void ) const;
    inline int32_t                  get_nb_indivs( void ) const;
    
    inline ae_individual* get_best_indiv( void ) const;
    inline ae_individual*	get_indiv_by_id( int32_t id ) const;
    inline ae_individual* get_indiv_by_rank( int32_t rank ) const;
    
    // Accessors to output manager stuff
    inline int32_t	get_backup_step(void) const;
    inline int32_t	get_big_backup_step(void) const;
    inline bool         get_record_tree( void ) const;
    inline int32_t      get_tree_step( void ) const;
    inline ae_tree_mode get_tree_mode( void ) const;
    inline ae_tree*     get_tree( void ) const;
    
    // =======================================================================
    //                          Accessors: setters
    // =======================================================================
    inline void set_first_gener( int32_t first_gener );
    inline void set_nb_gener( int32_t nb_gener );
    //~ inline void set_min_genome_length( int32_t min_genome_length );
    //~ inline void set_max_genome_length( int32_t max_genome_length );
    inline void set_spatial_structure(  int16_t grid_width,
                                        int16_t grid_height,
                                        ae_jumping_mt* prng );

    // =======================================================================
    //                                 Operators
    // =======================================================================

    // =======================================================================
    //                               Public Methods
    // =======================================================================
    void write_setup_files( void );
    void save( void ) const;
    void save_copy( char* dir, int32_t num_gener = 0 ) const;
    void load( int32_t first_gener, bool use_text_files, bool verbose, bool to_be_run = true);
    void load( int32_t first_gener,
               char* exp_setup_file_name,
               char* out_prof_file_name,
               char* env_file_name,
               char* pop_file_name,
               char* sel_file_name,
               char* sp_struct_file_name,
               bool verbose , 
               bool to_be_run = true);
    void run_evolution( void );
    virtual void display( void ) {};

    // =======================================================================
    //                              Public Attributes
    // =======================================================================
    





  protected :

    // =======================================================================
    //                            Forbidden Constructors
    // =======================================================================
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


    // =======================================================================
    //                              Protected Methods
    // =======================================================================
    inline void step_to_next_generation( void );
  
    void create_missing_directories( const char* dir = "." ) const;
    void open_backup_files( gzFile& env_file,
                            gzFile& pop_file,
                            gzFile& sel_file,
                            gzFile& sp_struct_file,
                            int32_t num_gener,
                            const char mode[3],
                            const char* dir = "." ) const;
    void close_backup_files(  gzFile& env_file,
                              gzFile& pop_file,
                              gzFile& sel_file,
                              gzFile& sp_struct_file ) const;
    void open_setup_files(  gzFile& exp_s_gzfile, FILE*& exp_s_txtfile,
                            gzFile& out_p_gzfile, FILE*& out_p_txtfile,
                            int32_t num_gener,
                            const char mode[3],
                            const char* dir = "." ) const;
    void close_setup_files(
            gzFile& exp_s_gzfile, FILE* exp_s_txtfile,
            gzFile& out_p_gzfile, FILE* out_p_txtfile ) const;

    // =======================================================================
    //                             Protected Attributes
    // =======================================================================
    // ---------------------------------------------------- Experimental setup
    ae_exp_setup* _exp_s;
    
    // ------------------------------------------------------------ Population
    ae_population* _pop;
      
    // ----------------------------------------------------------- Environment
    ae_environment* _env;
    
    // ----------------------------------------------------- Spatial structure
    ae_spatial_structure* _spatial_structure;
    
    // -------------------------------------------------------- Output manager
    ae_output_manager* _output_m;
    
    // -------------------------- Generation numbers (first, last and current)
    int32_t _first_gener;
    int32_t _last_gener;
    int32_t _num_gener;
    
    // Set to true when ctrl-Q is received. Will cause the simulation
    // to be ended after current generation is completed
    bool _quit_signal_received;
};


// ===========================================================================
//                             Getters' definitions
// ===========================================================================
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
  return _env;
}

inline ae_selection* ae_exp_manager::get_sel( void ) const
{
  return get_exp_s()->get_sel();
}

inline int16_t ae_exp_manager::get_nb_env_segments( void ) const
{
  return get_env()->get_nb_segments();
}

inline ae_selection_scheme ae_exp_manager::get_selection_scheme( void ) const
{
  return get_sel()->get_selection_scheme();
}

inline double ae_exp_manager::get_selection_pressure( void ) const
{
  return get_sel()->get_selection_pressure();
}

// Global settings
inline bool ae_exp_manager::is_spatially_structured( void ) const
{
  return ( get_spatial_structure() != NULL );
}

inline ae_spatial_structure* ae_exp_manager::get_spatial_structure( void ) const
{
  return _spatial_structure;
}

inline ae_grid_cell* ae_exp_manager::get_grid_cell( int16_t x, int16_t y ) const
{
  return get_spatial_structure()->get_grid_cell( x, y );
}

inline int16_t ae_exp_manager::get_grid_width( void ) const
{
  return get_spatial_structure()->get_grid_width();
}

inline int16_t ae_exp_manager::get_grid_height( void ) const
{
  return get_spatial_structure()->get_grid_height();
}

inline ae_grid_cell*** ae_exp_manager::get_pop_grid( void ) const
{
  return get_spatial_structure()->get_pop_grid();
}

inline bool ae_exp_manager::get_with_HT( void ) const
{
  return get_exp_s()->get_with_HT();
}

inline double ae_exp_manager::get_HT_ins_rate( void ) const
{
  return get_exp_s()->get_HT_ins_rate();
}

inline double ae_exp_manager::get_HT_repl_rate( void ) const
{
  return get_exp_s()->get_HT_repl_rate();
}

inline bool ae_exp_manager::get_with_plasmids( void ) const
{
  return get_exp_s()->get_with_plasmids();
}

inline double ae_exp_manager::get_prob_plasmid_HT( void ) const
{
  return get_exp_s()->get_prob_plasmid_HT();
}

inline double ae_exp_manager::get_tune_donor_ability( void ) const
{
  return get_exp_s()->get_tune_donor_ability();
}

inline double ae_exp_manager::get_tune_recipient_ability( void ) const
{
  return get_exp_s()->get_tune_recipient_ability();
}

inline double ae_exp_manager::get_donor_cost( void ) const
{
  return get_exp_s()->get_donor_cost();
}

inline double ae_exp_manager::get_recipient_cost( void ) const
{
  return get_exp_s()->get_recipient_cost();
}

inline bool ae_exp_manager::get_swap_GUs( void ) const
{
  return get_exp_s()->get_swap_GUs();
}

inline bool ae_exp_manager::get_with_secretion( void ) const
{
  return get_exp_s()->get_with_secretion();
}

inline double ae_exp_manager::get_secretion_contrib_to_fitness( void ) const
{
  return get_exp_s()->get_secretion_contrib_to_fitness();
}

inline double ae_exp_manager::get_secretion_cost( void ) const
{
  return get_exp_s()->get_secretion_cost();
}

//~ inline bool ae_exp_manager::get_with_alignments( void ) const
//~ {
  //~ return _exp_s->get_with_alignments();
//~ }


// Accessors to population stuff
inline int32_t ae_exp_manager::get_nb_indivs( void ) const
{
  return get_pop()->get_nb_indivs();
}

inline ae_individual* ae_exp_manager::get_best_indiv( void ) const
{
  return get_pop()->get_best();
}

inline ae_list<ae_individual*>* ae_exp_manager::get_indivs( void ) const
{
  return get_pop()->get_indivs();
}

inline ae_individual * ae_exp_manager::get_indiv_by_id( int32_t id ) const
{
  return get_pop()->get_indiv_by_id( id );
}

inline ae_individual * ae_exp_manager::get_indiv_by_rank( int32_t rank ) const
{
  return get_pop()->get_indiv_by_rank( rank );
}


// Accessors to output manager stuff
inline int32_t ae_exp_manager::get_backup_step(void) const
{
	return get_output_m()->get_backup_step();
}

inline int32_t ae_exp_manager::get_big_backup_step(void) const
{
	return get_output_m()->get_big_backup_step();
}

inline bool ae_exp_manager::get_record_tree( void ) const
{
	return get_output_m()->get_record_tree();
}

inline int32_t ae_exp_manager::get_tree_step( void ) const
{
	return get_output_m()->get_tree_step();
}

inline ae_tree_mode ae_exp_manager::get_tree_mode( void ) const
{
	return get_output_m()->get_tree_mode();
}

inline ae_tree* ae_exp_manager::get_tree( void ) const
{
	return get_output_m()->get_tree();
}

// ===========================================================================
//                             Setters' definitions
// ===========================================================================
inline void ae_exp_manager::set_first_gener( int32_t first_gener )
{
  _num_gener = _first_gener = first_gener;
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

inline void ae_exp_manager::set_spatial_structure( int16_t grid_width,
                                                   int16_t grid_height,
                                                   ae_jumping_mt* prng )
{
  _spatial_structure = new ae_spatial_structure();
  _spatial_structure->set_grid_size( grid_width, grid_height );
  _spatial_structure->set_prng( prng );
}

// ===========================================================================
//                            Operators' definitions
// ===========================================================================

// ===========================================================================
//                         Inline methods' definition
// ===========================================================================
inline void ae_exp_manager::step_to_next_generation( void )
{
  // Apply environmental variation
  _env->apply_variation();
  
  // Apply environmental noise
  _env->apply_noise();
  
  _exp_s->step_to_next_generation();
  _num_gener++;
}


#endif // __AE_EXP_MANAGER_H__
