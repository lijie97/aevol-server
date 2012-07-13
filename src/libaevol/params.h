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


/*! \class params
    \brief Singleton. Contains all the parameters needed to create an experiment.
*/


#ifndef __AE_PARAMS_H__
#define __AE_PARAMS_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_macros.h>
#include <ae_enums.h>
#include <ae_list.h>
#include <ae_gaussian.h>
#include <ae_point_2d.h>




// =================================================================
//                          Class declarations
// =================================================================






class params : public ae_object
{
  friend class param_loader;
  
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    params( void );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~params( void );

    // =================================================================
    //                         Accessors: getters
    // =================================================================
  
    // ---------------------------------- Number of generations to be simulated
    inline int32_t get_nb_gener( void ) const;
  
    // ----------------------------------------- PseudoRandom Number Generators
    inline int32_t  get_seed( void )                    const;
    inline int32_t  get_env_seed( void )                const;
  
    // ------------------------------------------------------------ Constraints
    inline int32_t  get_min_genome_length( void ) const;
    inline int32_t  get_max_genome_length( void ) const;
    inline double   get_w_max( void ) const;
    
    // ----------------------------------------------------- Initial conditions
    inline int32_t  get_initial_genome_length( void )   const;
    inline int8_t   get_init_method( void )             const;
    inline int32_t  get_init_pop_size( void )           const;
    
    // ------------------------------------------------------------ Environment
    inline ae_list* get_env_gaussians( void )           const;
    inline ae_list* get_env_custom_points( void )       const;
    inline int16_t  get_env_sampling( void )            const;
    
    // ---------------------------------------- Environment x-axis segmentation
    inline bool                 get_env_axis_is_segmented( void )       const;
    inline int16_t              get_env_axis_nb_segments( void )        const;
    inline double*              get_env_axis_segment_boundaries( void ) const;
    inline ae_env_axis_feature* get_env_axis_features( void )           const;
    inline bool                 get_env_axis_separate_segments( void )  const;
    
    // -------------------------------------------------- Environment variation
    inline ae_env_var get_env_var_method( void )  const;
    inline double     get_env_var_sigma( void )   const;
    inline int32_t    get_env_var_tau( void )     const;
    
    // ------------------------------------------------------ Environment noise

    // --------------------------------------------------------- Mutation rates
    inline double   get_point_mutation_rate( void ) const;
    inline double   get_small_insertion_rate( void ) const;
    inline double   get_small_deletion_rate( void ) const;
    inline int16_t  get_max_indel_size( void ) const;

    // -------------------------------------------- Rearrangements and Transfer
    inline bool   get_with_4pts_trans( void ) const;
    inline bool   get_with_alignments( void ) const;
    inline bool   get_with_HT( void ) const;
    inline bool   get_swap_GUs( void ) const;
    inline double get_HT_ins_rate( void ) const;
    inline double get_HT_repl_rate( void ) const;

    // ------------------------------ Rearrangement rates (without alignements)
    inline double get_duplication_rate( void ) const;
    inline double get_deletion_rate( void ) const;
    inline double get_translocation_rate( void ) const;
    inline double get_inversion_rate( void ) const;

    // --------------------------------- Rearrangement rates (with alignements)
    inline double get_neighbourhood_rate( void ) const;
    inline double get_duplication_proportion( void ) const;
    inline double get_deletion_proportion( void ) const;
    inline double get_translocation_proportion( void ) const;
    inline double get_inversion_proportion( void ) const;

    // ------------------------------------------------------------ Alignements
    inline ae_align_fun_shape get_align_fun_shape( void ) const;
    inline double  get_align_sigm_lambda( void ) const;
    inline int16_t get_align_sigm_mean( void ) const;
    inline int16_t get_align_lin_min( void ) const;
    inline int16_t get_align_lin_max( void ) const;
    
    // Maximum shift of one seq on the other
    inline int16_t get_align_max_shift( void ) const;
    // Work zone half length
    inline int16_t get_align_w_zone_h_len( void ) const;
    // Corresponding residues match bonus
    inline int16_t get_align_match_bonus( void ) const;
    // Corresponding residues mismatch cost
    inline int16_t get_align_mismatch_cost( void ) const;

    // -------------------------------------------------------------- Selection
    inline ae_selection_scheme  get_selection_scheme( void ) const;
    inline double               get_selection_pressure( void ) const;
    
    // ------------------------------------------------------ Spatial structure
    inline bool     is_spatially_structured( void ) const;
    inline int16_t  get_grid_width( void ) const;
    inline int16_t  get_grid_height( void ) const;
    inline int32_t  get_migration_number( void ) const;
    
    // -------------------------------------------------------------- Secretion
    inline bool   get_use_secretion( void ) const;
    inline double get_secretion_contrib_to_fitness( void ) const;      // (0,1)
    // proportion that diffuses into each cell, every generation, 0 for no diffusion
    inline double get_secretion_diffusion_prop( void ) const;
    // proportion of secreted substance that degrades every generation generation
    inline double get_secretion_degradation_prop( void ) const;
    // cost of secreting the compound, as a proportion of the amount secreted
    inline double get_secretion_cost( void ) const;
    inline double get_secretion_init( void )          const;
    
    // --------------------------------------------------------------- Plasmids
    inline bool     get_allow_plasmids( void ) const;
    inline int32_t  get_plasmid_initial_length( void ) const;
    inline int32_t  get_plasmid_initial_gene( void ) const;
    inline int32_t  get_plasmid_minimal_length( void ) const;
    inline bool     get_with_plasmid_HT( void ) const;
    inline double   get_prob_plasmid_HT( void ) const;
    inline int16_t  get_nb_plasmid_HT( void ) const;
    inline bool     get_compute_phen_contrib_by_GU( void ) const;
    
    // ------------------------------------------------------- Translation cost
    inline double get_translation_cost( void ) const;
    
    // ---------------------------------------------------------------- Outputs
    //Stats
    inline int8_t get_stats( void ) const;
    inline bool   get_delete_old_stats( void ) const;
    
    // Backups
    inline int32_t get_backup_step( void ) const;
    inline int32_t get_big_backup_step( void ) const;
    
    // Tree
    inline bool         get_record_tree( void ) const;
    inline int32_t      get_tree_step( void ) const;
    inline ae_tree_mode get_tree_mode( void ) const;
    
    // Logs
    inline int8_t get_logs( void ) const;
    
    // Dumps
    inline bool    get_make_dumps( void ) const;
    inline int32_t get_dump_step( void ) const;
    
    // Other
    inline bool    get_more_stats( void ) const;
    

    #ifdef __REGUL
      // ------------------------------------------------------- Binding matrix
      inline double get_binding_zeros_percentage( void ) const;
    #endif

    // =================================================================
    //                         Accessors: setters
    // =================================================================
  
    // ---------------------------------- Number of generations to be simulated
    inline void set_nb_gener( int32_t nb_gener );
    
    // ----------------------------------------- PseudoRandom Number Generators
    inline void set_seed( int32_t seed );
    inline void set_env_seed( int32_t env_seed );
    
    // ------------------------------------------------------------ Constraints
    inline void set_min_genome_length( int32_t min_genome_length );
    inline void set_max_genome_length( int32_t max_genome_length );
    inline void set_w_max( double w_max );
    
    // ----------------------------------------------------- Initial conditions
    inline void set_initial_genome_length( int32_t init_g_len );
    inline void set_init_method( int8_t init_meth );
    inline void set_init_pop_size( int32_t init_pop_size );
    
    // ------------------------------------------------------------ Environment
    inline void set_env_gaussians( ae_list* env_gaussians );
    inline void set_env_custom_points( ae_list* env_custom_points );
    inline void add_env_gaussian( ae_gaussian* new_gaussian );
    inline void add_env_custom_point( ae_point_2d* new_point );
    inline void set_env_sampling( int16_t env_sampling );
    
    // ---------------------------------------- Environment x-axis segmentation
    inline void set_env_axis_is_segmented( bool env_axis_is_segmented );
    inline void set_env_axis_nb_segments( int16_t env_axis_nb_segments );
    inline void set_env_axis_segment_boundaries( double* env_axis_segment_boundaries ); // set_env_axis_nb_segments first
    inline void set_env_axis_features( ae_env_axis_feature* env_axis_features );        // set_env_axis_nb_segments first
    inline void set_env_axis_separate_segments( bool env_axis_separate_segments );
    
    // -------------------------------------------------- Environment variation
    inline void set_env_var_method( ae_env_var env_var_method );
    inline void set_env_var_sigma( double env_var_sigma );
    inline void set_env_var_tau( int32_t env_var_tau );
    
    // ------------------------------------------------------ Environment noise

    // --------------------------------------------------------- Mutation rates
    inline void set_point_mutation_rate( double point_mutation_rate );
    inline void set_small_insertion_rate( double small_insertion_rate );
    inline void set_small_deletion_rate( double small_deletion_rate );
    inline void set_max_indel_size( int16_t max_indel_size );

    // -------------------------------------------- Rearrangements and Transfer
    inline void set_with_4pts_trans( bool with_4pts_trans );
    inline void set_with_alignments( bool with_alignments );
    inline void set_with_HT( bool with_HT );
    inline void set_swap_GUs( bool swap_GUs );
    inline void set_HT_ins_rate( double HT_ins_rate );
    inline void set_HT_repl_rate( double HT_repl_rate );

    // ------------------------------ Rearrangement rates (without alignements)
    inline void set_duplication_rate( double duplication_rate );
    inline void set_deletion_rate( double deletion_rate );
    inline void set_translocation_rate( double translocation_rate );
    inline void set_inversion_rate( double inversion_rate );

    // --------------------------------- Rearrangement rates (with alignements)
    inline void set_neighbourhood_rate( double neighbourhood_rate );
    inline void set_duplication_proportion( double duplication_proportion );
    inline void set_deletion_proportion( double deletion_proportion );
    inline void set_translocation_proportion( double translocation_proportion );
    inline void set_inversion_proportion( double inversion_proportion );

    // ------------------------------------------------------------ Alignements
    inline void set_align_fun_shape( ae_align_fun_shape align_fun_shape );
    inline void set_align_sigm_lambda( double align_sigm_lambda );
    inline void set_align_sigm_mean( int16_t align_sigm_mean );
    inline void set_align_lin_min( int16_t align_lin_min );
    inline void set_align_lin_max( int16_t align_lin_max );
    
    inline void set_align_max_shift( int16_t align_max_shift );
    inline void set_align_w_zone_h_len( int16_t align_w_zone_h_len );
    inline void set_align_match_bonus( int16_t align_match_bonus );
    inline void set_align_mismatch_cost( int16_t align_mismatch_cost );

    // -------------------------------------------------------------- Selection
    inline void set_selection_scheme( ae_selection_scheme selection_scheme );
    inline void set_selection_pressure( double selection_pressure );
    
    // ------------------------------------------------------ Spatial structure
    inline void set_spatially_structured( bool structured );
    inline void set_grid_width( int16_t grid_width ); 
    inline void set_grid_height( int16_t grid_height );
    inline void set_migration_number( int32_t migration_number );
    
    // -------------------------------------------------------------- Secretion
    inline void set_use_secretion( bool use_secretion );
    inline void set_secretion_contrib_to_fitness( double secretion_contrib_to_fitness );
    inline void set_secretion_diffusion_prop( double secretion_diffusion_prop );
    inline void set_secretion_degradation_prop( double secretion_degradation_prop );
    inline void set_secretion_cost( double secretion_cost );
    inline void set_secretion_init( double secretion_init );
    
    // --------------------------------------------------------------- Plasmids
    inline void set_allow_plasmids( bool allow_plasmids );
    inline void set_plasmid_initial_length( int32_t plasmid_init_len );
    inline void set_plasmid_initial_gene( int32_t plasmid_init_gene );
    inline void set_plasmid_minimal_length( int32_t plasmid_minimal_length );
    inline void set_with_plasmid_HT( bool with_plasmid_HT );
    inline void set_prob_plasmid_HT( double prob_plasmid_HT );
    inline void set_nb_plasmid_HT( int16_t nb_plasmid_HT );
    inline void set_compute_phen_contrib_by_GU( bool compute_phen_contrib_by_GU );
    
    // ------------------------------------------------------- Translation cost
    inline void set_translation_cost( double translation_cost );
    
    // ---------------------------------------------------------------- Outputs
    // Stats
    inline void set_stats( int8_t stats );
    inline void set_delete_old_stats( bool delete_old_stats );
    
    // Backups
    inline void set_backup_step( int32_t backup_step );
    inline void set_big_backup_step( int32_t big_backup_step );
    
    // Tree
    inline void set_record_tree( bool record_tree );
    inline void set_tree_step( int32_t tree_step );
    inline void set_tree_mode( ae_tree_mode tree_mode );
    
    // Dumps
    inline void set_dump_step( int32_t dump_step );
    
    // Logs
    inline void set_logs( int8_t logs );
    
    // Other
    inline void set_more_stats( bool more_stats );

    #ifdef __REGUL
      // ------------------------------------------------------- Binding matrix
      inline void set_binding_zeros_percentage( double zeros_percentage );
    #endif

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    //~ void write_to_backup( gzFile* backup_file );
    //~ void read_from_backup( gzFile* backup_file, bool verbose = false );
    void print_to_file( FILE* file );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*params( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    params( const params &model )
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
  
    // ---------------------------------- Number of generations to be simulated
    int32_t _nb_gener;
    
    // ----------------------------------------- PseudoRandom Number Generators
    // Seed for the individuals' random generators
    int32_t _seed;
    // Seed for the environment's random generator (environmental variation)
    int32_t _env_seed;
    
    // ------------------------------------------------------------ Constraints
    int32_t _min_genome_length;
    int32_t _max_genome_length;
    double  _w_max;
    
    // ----------------------------------------------------- Initial conditions
    int32_t  _initial_genome_length;
    int8_t   _init_method;
    int32_t  _init_pop_size;
    
    // ------------------------------------------------------------ Environment
    ae_list* _env_gaussians;
    ae_list* _env_custom_points;
    int16_t  _env_sampling;
    
    // ---------------------------------------- Environment x-axis segmentation
    // Whether the environment is segmented
    bool _env_axis_is_segmented;
    // Number of x-axis segments
    int16_t _env_axis_nb_segments;
    // x-axis segment boundaries (sorted -- including MIN_X and MAX_X)
    double* _env_axis_segment_boundaries;
    // x-axis segment features
    ae_env_axis_feature* _env_axis_features;
    // Whether to automatically separate segments
    bool _env_axis_separate_segments;
    
    // -------------------------------------------------- Environment variation
    ae_env_var  _env_var_method;
    double      _env_var_sigma;
    int32_t     _env_var_tau;
    
    // ------------------------------------------------------ Environment noise // TODO

    // --------------------------------------------------------- Mutation rates
    double  _point_mutation_rate;
    double  _small_insertion_rate;
    double  _small_deletion_rate;
    int16_t _max_indel_size;

    // -------------------------------------------- Rearrangements and Transfer
    bool    _with_4pts_trans;
    bool    _with_alignments;
    bool    _with_HT;
    bool    _swap_GUs;
    double  _HT_ins_rate;
    double  _HT_repl_rate;

    // ------------------------------ Rearrangement rates (without alignements)
    double _duplication_rate;
    double _deletion_rate;
    double _translocation_rate;
    double _inversion_rate;

    // --------------------------------- Rearrangement rates (with alignements)
    double _neighbourhood_rate;
    double _duplication_proportion;
    double _deletion_proportion;
    double _translocation_proportion;
    double _inversion_proportion;

    // ------------------------------------------------------------ Alignements
    ae_align_fun_shape _align_fun_shape;
    double  _align_sigm_lambda;
    int16_t _align_sigm_mean;
    int16_t _align_lin_min;
    int16_t _align_lin_max;
    
    int16_t _align_max_shift;     // Maximum shift of one seq on the other
    int16_t _align_w_zone_h_len;  // Work zone half length
    int16_t _align_match_bonus;   // Corresponding residues match bonus
    int16_t _align_mismatch_cost; // Corresponding residues mismatch cost

    // -------------------------------------------------------------- Selection
    ae_selection_scheme  _selection_scheme;
    double               _selection_pressure;
    
    // ------------------------------------------------------ Spatial structure
    bool     _spatially_structured;
    int16_t  _grid_width; 
    int16_t  _grid_height;
    int32_t  _migration_number; // TODO : explain
    
    // -------------------------------------------------------------- Secretion
    bool   _use_secretion;
    // Proportion of the fitness contributed by secretion
    double _secretion_contrib_to_fitness;      // (0,1)
    // Proportion that diffuses into each cell, every generation
    // (0 for no diffusion)
    double _secretion_diffusion_prop;
    // Proportion of secreted substance that degrades every generation
    double _secretion_degradation_prop;
    // Cost of secreting the compound, as a proportion of the amount secreted
    double _secretion_cost;
    // Starting configuration of secretion grid
    // 0, all are 0; 1, point source of secreted compund
    double _secretion_init;
    
    // --------------------------------------------------------------- Plasmids
    bool      _allow_plasmids;
    int32_t   _plasmid_initial_length;
    int32_t   _plasmid_initial_gene;
    int32_t   _plasmid_minimal_length;
    bool      _with_plasmid_HT;
    double    _prob_plasmid_HT;
    // Max number of plasmids a cell can send in its neighbourhood per generation
    int16_t   _nb_plasmid_HT;
    bool      _compute_phen_contrib_by_GU;
    
    // ------------------------------------------------------- Translation cost
    double _translation_cost;
    
    // ---------------------------------------------------------------- Outputs
    // Stats
    int8_t _stats;
    // Whether to delete the existing statistics file
    // (otherwise kept with the suffix ".old")
    bool _delete_old_stats;
    
    // Backups
    int32_t _backup_step;
    int32_t _big_backup_step;
    
    // Tree
    bool         _record_tree;
    int32_t      _tree_step;
    ae_tree_mode _tree_mode;
    
    // Dumps // TODO : explain
    bool    _make_dumps;
    int32_t _dump_step;
    
    // Logs
    int8_t _logs;
    
    // Other
    bool _more_stats;  // TODO : explain

    #ifdef __REGUL
      // ------------------------------------------------------- Binding matrix
      double _binding_zeros_percentage;
    #endif
};





// =====================================================================
//                          Getters' definitions
// =====================================================================

// --------------------------------------------- PseudoRandom Number Generators
inline int32_t params::get_seed( void ) const
{
  return _seed;
}

inline int32_t params::get_env_seed( void ) const
{
  return _env_seed;
}

// ---------------------------------------------------------------- Constraints
inline int32_t params::get_min_genome_length( void ) const
{
  return _min_genome_length;
}

inline int32_t params::get_max_genome_length( void ) const
{
  return _max_genome_length;
}

inline double params::get_w_max( void ) const
{
  return _w_max;
}

// --------------------------------------------------------- Initial conditions
inline int32_t params::get_initial_genome_length( void ) const
{
  return _initial_genome_length;
}

inline int8_t params::get_init_method( void ) const
{
  return _init_method;
}

inline int32_t params::get_init_pop_size( void ) const
{
  return _init_pop_size;
}

// ---------------------------------------------------------------- Environment
inline ae_list* params::get_env_gaussians( void ) const
{
  return _env_gaussians;
}

inline ae_list* params::get_env_custom_points( void ) const
{
  return _env_custom_points;
}

inline int16_t params::get_env_sampling( void ) const
{
  return _env_sampling;
}

// -------------------------------------------- Environment x-axis segmentation
    
inline bool params::get_env_axis_is_segmented( void ) const
{
  return _env_axis_is_segmented;
}

inline int16_t params::get_env_axis_nb_segments( void ) const
{
  return _env_axis_nb_segments;
}

inline double* params::get_env_axis_segment_boundaries( void ) const
{
  return _env_axis_segment_boundaries;
}

inline ae_env_axis_feature* params::get_env_axis_features( void ) const
{
  return _env_axis_features;
}

inline bool params::get_env_axis_separate_segments( void ) const
{
  return _env_axis_separate_segments;
}

// ------------------------------------------------------ Environment variation
inline ae_env_var params::get_env_var_method( void ) const
{
  return _env_var_method;
}

inline double params::get_env_var_sigma( void ) const
{
  return _env_var_sigma;
}

inline int32_t params::get_env_var_tau( void ) const
{
  return _env_var_tau;
}

// ---------------------------------------------------------- Environment noise

// ------------------------------------------------------------- Mutation rates
inline double params::get_point_mutation_rate( void ) const
{
  return _point_mutation_rate;
}

inline double params::get_small_insertion_rate( void ) const
{
  return _small_insertion_rate;
}

inline double params::get_small_deletion_rate( void ) const
{
  return _small_deletion_rate;
}

inline int16_t params::get_max_indel_size( void ) const
{
  return _max_indel_size;
}

// ------------------------------------------------ Rearrangements and Transfer
inline bool params::get_with_4pts_trans( void ) const
{
  return _with_4pts_trans;
}

inline bool params::get_with_alignments( void ) const
{
  return _with_alignments;
}

inline bool params::get_with_HT( void ) const
{
  return _with_HT;
}

inline bool params::get_swap_GUs( void ) const
{
  return _swap_GUs;
}

inline double params::get_HT_ins_rate( void ) const
{
  return _HT_ins_rate;
}

inline double params::get_HT_repl_rate( void ) const
{
  return _HT_repl_rate;
}

// ---------------------------------- Rearrangement rates (without alignements)
inline double params::get_duplication_rate( void ) const
{
  return _duplication_rate;
}

inline double params::get_deletion_rate( void ) const
{
  return _deletion_rate;
}

inline double params::get_translocation_rate( void ) const
{
  return _translocation_rate;
}

inline double params::get_inversion_rate( void ) const
{
  return _inversion_rate;
}

// ------------------------------------- Rearrangement rates (with alignements)
inline double params::get_neighbourhood_rate( void ) const
{
  return _neighbourhood_rate;
}

inline double params::get_duplication_proportion( void ) const
{
  return _duplication_proportion;
}

inline double params::get_deletion_proportion( void ) const
{
  return _deletion_proportion;
}

inline double params::get_translocation_proportion( void ) const
{
  return _translocation_proportion;
}

inline double params::get_inversion_proportion( void ) const
{
  return _inversion_proportion;
}

// ---------------------------------------------------------------- Alignements
inline ae_align_fun_shape params::get_align_fun_shape( void ) const
{
  return _align_fun_shape;
}

inline double params::get_align_sigm_lambda( void ) const
{
  return _align_sigm_lambda;
}

inline int16_t params::get_align_sigm_mean( void ) const
{
  return _align_sigm_mean;
}

inline int16_t params::get_align_lin_min( void ) const
{
  return _align_lin_min;
}

inline int16_t params::get_align_lin_max( void ) const
{
  return _align_lin_max;
}

inline int16_t params::get_align_max_shift( void ) const
{
  return _align_max_shift;
}

inline int16_t params::get_align_w_zone_h_len( void ) const
{
  return _align_w_zone_h_len;
}

inline int16_t params::get_align_match_bonus( void ) const
{
  return _align_match_bonus;
}

inline int16_t params::get_align_mismatch_cost( void ) const
{
  return _align_mismatch_cost;
}

// ------------------------------------------------------------------ Selection
inline ae_selection_scheme params::get_selection_scheme( void ) const
{
  return _selection_scheme;
}

inline double params::get_selection_pressure( void ) const
{
  return _selection_pressure;
}

// ---------------------------------------------------------- Spatial structure
inline bool params::is_spatially_structured( void ) const
{
  return _spatially_structured;
}

inline int16_t params::get_grid_width( void ) const
{
  return _grid_width;
}

inline int16_t params::get_grid_height( void ) const
{
  return _grid_height;
}

inline int32_t params::get_migration_number( void ) const
{
  return _migration_number;
}


// ------------------------------------------------------------------ Secretion
inline bool params::get_use_secretion( void ) const
{
  return _use_secretion;
}

inline double params::get_secretion_contrib_to_fitness( void ) const
{
  return _secretion_contrib_to_fitness;
}

inline double params::get_secretion_diffusion_prop( void ) const
{
  return _secretion_diffusion_prop;
}

inline double params::get_secretion_degradation_prop( void ) const
{
  return _secretion_degradation_prop;
}

inline double params::get_secretion_cost( void ) const
{
  return _secretion_cost;
}

inline double params::get_secretion_init( void ) const
{
  return _secretion_init;
}

// ------------------------------------------------------------------- Plasmids
inline bool params::get_allow_plasmids( void ) const
{
  return _allow_plasmids;
}

inline int32_t params::get_plasmid_initial_length( void ) const
{
  return _plasmid_initial_length;
}

inline int32_t params::get_plasmid_initial_gene( void ) const
{
  return _plasmid_initial_gene;
}

inline int32_t params::get_plasmid_minimal_length( void ) const
{
  return _plasmid_minimal_length;
}

inline bool params::get_with_plasmid_HT( void ) const
{
  return _with_plasmid_HT;
}

inline double params::get_prob_plasmid_HT( void ) const
{
  return _prob_plasmid_HT;
}

inline int16_t params::get_nb_plasmid_HT( void ) const
{
  return _nb_plasmid_HT;
}

inline bool params::get_compute_phen_contrib_by_GU( void ) const
{
  return _compute_phen_contrib_by_GU;
}

// ----------------------------------------------------------- Translation cost
inline double params::get_translation_cost( void ) const
{
  return _translation_cost;
}

// -------------------------------------------------------------------- Outputs
// Backups
inline int32_t params::get_backup_step( void ) const
{
  return _backup_step;
}

inline int32_t params::get_big_backup_step( void ) const
{
  return _big_backup_step;
}

// Tree
inline bool params::get_record_tree( void ) const
{
  return _record_tree;
}

inline int32_t params::get_tree_step( void ) const
{
  return _tree_step;
}

inline ae_tree_mode params::get_tree_mode( void ) const
{
  return _tree_mode;
}
    
// Logs
inline int8_t params::get_logs( void ) const
{
  return _logs;
}

// Other
inline bool params::get_more_stats( void ) const
{
  return _more_stats;
}

inline bool params::get_make_dumps( void ) const
{
  return _make_dumps;
}

inline int32_t params::get_dump_step( void ) const
{
  return _dump_step;
}

inline int8_t params::get_stats( void ) const
{
  return _stats;
}

inline bool params::get_delete_old_stats( void ) const
{
  return _delete_old_stats;
}

// =====================================================================
//                          Setters' definitions
// =====================================================================

// ---------------------------------- Number of generations to be simulated
inline void params::set_nb_gener( int32_t nb_gener )
{
  _nb_gener = nb_gener;
}
    
// --------------------------------------------- PseudoRandom Number Generators
inline void params::set_seed( int32_t seed )
{
  _seed = seed;
}

inline void params::set_env_seed( int32_t env_seed )
{
  _env_seed = env_seed;
}
    
// ------------------------------------------------------------ Constraints
inline void params::set_min_genome_length( int32_t min_genome_length )
{
  _min_genome_length = min_genome_length;
}

inline void params::set_max_genome_length( int32_t max_genome_length )
{
  _max_genome_length = max_genome_length;
}

inline void params::set_w_max( double w_max )
{
  _w_max = w_max;
}

// --------------------------------------------------------- Initial conditions
inline void params::set_initial_genome_length( int32_t init_g_len )
{
  _initial_genome_length = init_g_len;
}

inline void params::set_init_method( int8_t init_meth )
{
  _init_method = init_meth;
}

inline void params::set_init_pop_size( int32_t init_pop_size )
{
  _init_pop_size = init_pop_size;
}

// ---------------------------------------------------------------- Environment
inline void params::set_env_gaussians( ae_list* env_gaussians )
{
  _env_gaussians = env_gaussians;
}

inline void params::set_env_custom_points( ae_list* env_custom_points )
{
  _env_custom_points = env_custom_points;
}

inline void params::add_env_gaussian( ae_gaussian* new_gaussian )
{
  _env_gaussians->add( new_gaussian );
}

inline void params::add_env_custom_point( ae_point_2d* new_point )
{
  _env_custom_points->add( new_point );
}

inline void params::set_env_sampling( int16_t env_sampling )
{
  _env_sampling = env_sampling;
}

// -------------------------------------------- Environment x-axis segmentation
inline void params::set_env_axis_is_segmented( bool env_axis_is_segmented )
{
  _env_axis_is_segmented = env_axis_is_segmented;
}

inline void params::set_env_axis_nb_segments( int16_t env_axis_nb_segments )
{
  _env_axis_nb_segments = env_axis_nb_segments;
}

inline void params::set_env_axis_segment_boundaries( double* env_axis_segment_boundaries )
{
  _env_axis_segment_boundaries = env_axis_segment_boundaries;
}

inline void params::set_env_axis_features( ae_env_axis_feature* env_axis_features )
{
  _env_axis_features = env_axis_features;
}

inline void params::set_env_axis_separate_segments( bool env_axis_separate_segments )
{
  _env_axis_separate_segments = env_axis_separate_segments;
}

// ------------------------------------------------------ Environment variation
inline void params::set_env_var_method( ae_env_var env_var_method )
{
  _env_var_method = env_var_method;
}

inline void params::set_env_var_sigma( double env_var_sigma )
{
  _env_var_sigma = env_var_sigma;
}

inline void params::set_env_var_tau( int32_t env_var_tau )
{
  _env_var_tau = env_var_tau;
}

// ---------------------------------------------------------- Environment noise

// ------------------------------------------------------------- Mutation rates
inline void params::set_point_mutation_rate( double point_mutation_rate )
{
  _point_mutation_rate = point_mutation_rate;
}

inline void params::set_small_insertion_rate( double small_insertion_rate )
{
  _small_insertion_rate = small_insertion_rate;
}

inline void params::set_small_deletion_rate( double small_deletion_rate )
{
  _small_deletion_rate = small_deletion_rate;
}

inline void params::set_max_indel_size( int16_t max_indel_size )
{
  _max_indel_size = max_indel_size;
}

// ------------------------------------------------ Rearrangements and Transfer
inline void params::set_with_4pts_trans( bool with_4pts_trans )
{
  _with_4pts_trans = with_4pts_trans;
}

inline void params::set_with_alignments( bool with_alignments )
{
  _with_alignments = with_alignments;
}

inline void params::set_with_HT( bool with_HT )
{
  _with_HT = with_HT;
}

inline void params::set_swap_GUs( bool swap_GUs )
{
  _swap_GUs = swap_GUs;
}

inline void params::set_HT_ins_rate( double HT_ins_rate )
{
  _HT_ins_rate = HT_ins_rate;
}

inline void params::set_HT_repl_rate( double HT_repl_rate )
{
  _HT_repl_rate = HT_repl_rate;
}

// ---------------------------------- Rearrangement rates (without alignements)
inline void params::set_duplication_rate( double duplication_rate )
{
  _duplication_rate = duplication_rate;
}

inline void params::set_deletion_rate( double deletion_rate )
{
  _deletion_rate = deletion_rate;
}

inline void params::set_translocation_rate( double translocation_rate )
{
  _translocation_rate = translocation_rate;
}

inline void params::set_inversion_rate( double inversion_rate )
{
  _inversion_rate = inversion_rate;
}

// ------------------------------------- Rearrangement rates (with alignements)
inline void params::set_neighbourhood_rate( double neighbourhood_rate )
{
  _neighbourhood_rate = neighbourhood_rate;
}

inline void params::set_duplication_proportion( double duplication_proportion )
{
  _duplication_proportion = duplication_proportion;
}

inline void params::set_deletion_proportion( double deletion_proportion )
{
  _deletion_proportion = deletion_proportion;
}

inline void params::set_translocation_proportion( double translocation_proportion )
{
  _translocation_proportion = translocation_proportion;
}

inline void params::set_inversion_proportion( double inversion_proportion )
{
  _inversion_proportion = inversion_proportion;
}

// ---------------------------------------------------------------- Alignements
inline void params::set_align_fun_shape( ae_align_fun_shape align_fun_shape )
{
  _align_fun_shape = align_fun_shape;
}

inline void params::set_align_sigm_lambda( double align_sigm_lambda )
{
  _align_sigm_lambda = align_sigm_lambda;
}

inline void params::set_align_sigm_mean( int16_t align_sigm_mean )
{
  _align_sigm_mean = align_sigm_mean;
}

inline void params::set_align_lin_min( int16_t align_lin_min )
{
  _align_lin_min = align_lin_min;
}

inline void params::set_align_lin_max( int16_t align_lin_max )
{
  _align_lin_max = align_lin_max;
}

inline void params::set_align_max_shift( int16_t align_max_shift )
{
  _align_max_shift = align_max_shift;
}

inline void params::set_align_w_zone_h_len( int16_t align_w_zone_h_len )
{
  _align_w_zone_h_len = align_w_zone_h_len;
}

inline void params::set_align_match_bonus( int16_t align_match_bonus )
{
  _align_match_bonus = align_match_bonus;
}

inline void params::set_align_mismatch_cost( int16_t align_mismatch_cost )
{
  _align_mismatch_cost = align_mismatch_cost;
}

// ------------------------------------------------------------------ Selection
inline void params::set_selection_scheme( ae_selection_scheme selection_scheme )
{
  _selection_scheme = selection_scheme;
}

inline void params::set_selection_pressure( double selection_pressure )
{
  _selection_pressure = selection_pressure;
}
    
// ---------------------------------------------------------- Spatial structure
inline void params::set_spatially_structured( bool structured )
{
  _spatially_structured = structured;
}

inline void params::set_grid_width( int16_t grid_width )
{
  _grid_width = grid_width;
}
 
inline void params::set_grid_height( int16_t grid_height )
{
  _grid_height = grid_height;
}

inline void params::set_migration_number( int32_t migration_number )
{
  _migration_number = migration_number;
}


// ------------------------------------------------------------------ Secretion
inline void params::set_use_secretion( bool use_secretion )
{
  _use_secretion = use_secretion;
}

inline void params::set_secretion_contrib_to_fitness( double secretion_contrib_to_fitness )
{
  _secretion_contrib_to_fitness = secretion_contrib_to_fitness;
}

inline void params::set_secretion_diffusion_prop( double secretion_diffusion_prop )
{
  _secretion_diffusion_prop = secretion_diffusion_prop;
}

inline void params::set_secretion_degradation_prop( double secretion_degradation_prop )
{
  _secretion_degradation_prop = secretion_degradation_prop;
}

inline void params::set_secretion_cost( double secretion_cost )
{
  _secretion_cost = secretion_cost;
}

inline void params::set_secretion_init( double secretion_init )
{
  _secretion_init = secretion_init;
}

// ------------------------------------------------------------------- Plasmids
inline void params::set_allow_plasmids( bool allow_plasmids )
{
  _allow_plasmids = allow_plasmids;
}

inline void params::set_plasmid_initial_length( int32_t plasmid_init_len )
{
  _plasmid_initial_length = plasmid_init_len;
}

inline void params::set_plasmid_initial_gene( int32_t plasmid_init_gene )
{
  _plasmid_initial_gene = plasmid_init_gene;
}

inline void params::set_plasmid_minimal_length( int32_t plasmid_minimal_length )
{
  _plasmid_minimal_length = plasmid_minimal_length;
}

inline void params::set_with_plasmid_HT( bool with_plasmid_HT )
{
  _with_plasmid_HT = with_plasmid_HT;
}

inline void params::set_prob_plasmid_HT( double prob_plasmid_HT )
{
  _prob_plasmid_HT = prob_plasmid_HT;
}

inline void params::set_nb_plasmid_HT( int16_t nb_plasmid_HT )
{
  _nb_plasmid_HT = nb_plasmid_HT;
}

inline void params::set_compute_phen_contrib_by_GU( bool compute_phen_contrib_by_GU )
{
  _compute_phen_contrib_by_GU = compute_phen_contrib_by_GU;
}

// ----------------------------------------------------------- Translation cost
inline void params::set_translation_cost( double translation_cost )
{
  _translation_cost = translation_cost;
}

// -------------------------------------------------------------------- Outputs
// Backups
inline void params::set_backup_step( int32_t backup_step )
{
  _backup_step = backup_step;
}

inline void params::set_big_backup_step( int32_t big_backup_step )
{
  _big_backup_step = big_backup_step;
}

// Tree
inline void params::set_record_tree( bool record_tree )
{
  _record_tree = record_tree;
}

inline void params::set_tree_step( int32_t tree_step )
{
  _tree_step = tree_step;
}

inline void params::set_tree_mode( ae_tree_mode tree_mode )
{
  _tree_mode = tree_mode;
}
    
// Logs
inline void params::set_logs( int8_t logs )
{
  _logs = logs;
}

// Other
inline void params::set_more_stats( bool more_stats )
{
  _more_stats = more_stats;
}

inline void params::set_dump_step( int32_t dump_step )
{
  if ( dump_step > 0 )
  {
    _make_dumps = true;
  }
  
  _dump_step = dump_step;
}

inline void params::set_stats( int8_t stats )
{
  _stats = stats;
}

inline void params::set_delete_old_stats( bool delete_old_stats )
{
  _delete_old_stats = delete_old_stats;
}

#ifdef __REGUL
// ------------------------------------------------------------- Binding matrix
inline void params::set_binding_zeros_percentage( double zeros_percentage )
{
  _binding_zeros_percentage = zeros_percentage;
}
#endif


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_PARAMS_H__
