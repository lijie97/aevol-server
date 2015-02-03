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


#ifndef __PARAM_LOADER_H__
#define __PARAM_LOADER_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>

#include <zlib.h>


// =================================================================
//                            Project Files
// =================================================================
#include <f_line.h>
#include <ae_params_mut.h>
#include <ae_jumping_mt.h>

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ae_macros.h>
#include <ae_enums.h>
#include <ae_list.h>
#include <ae_gaussian.h>
#include <point.h>

namespace aevol {
// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;
class ae_environment;
class ae_individual;


class param_loader
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    param_loader( const char* file_name );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~param_loader( void );

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void read_file( void );
    void load( ae_exp_manager* exp_m, bool verbose = false, char* chromosome = NULL, int32_t lchromosome = 0, char* plasmid = NULL, int32_t lplasmid = 0 );


    f_line* get_line( int32_t* );

    // =================================================================
    //                           Public Attributes
    // =================================================================

    void print_to_file( FILE* file );




  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    param_loader( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    param_loader( const param_loader &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    static void format_line( f_line*, char*, bool* );
    void interpret_line( f_line* line, int32_t cur_line );
    ae_individual* create_random_individual( ae_exp_manager* exp_m, ae_params_mut* param_mut, int32_t id ) const;
    ae_individual* create_random_individual_with_good_gene( ae_exp_manager* exp_m, ae_params_mut* param_mut, int32_t id ) const;
    ae_individual* create_clone( ae_individual* dolly, int32_t id ) const;



    // =================================================================
    //                          Protected Attributes
    // =================================================================
    //~ ae_exp_manager* _exp_m;

    ae_jumping_mt* _prng;

    char*   _param_file_name;
    FILE*   _param_file;

    // ----------------------------------------- PseudoRandom Number Generators
    // Seed for the selection random generator
    int32_t _seed;
    // Seed for the mutations random generator
    int32_t _mut_seed;
    // Seed for the stochasticity random generator
    int32_t _stoch_seed;
    // Seed for the environmental variation random generator
    int32_t _env_var_seed;
    // Seed for the environmental noise random generator
    int32_t _env_noise_seed;

    // ------------------------------------------------------------ Constraints
    int32_t _min_genome_length;
    int32_t _max_genome_length;
    double  _w_max;

    // ----------------------------------------------------- Initial conditions
    int32_t  _chromosome_initial_length;
    int8_t   _init_method;
    int32_t  _init_pop_size;
    char*    _strain_name;

    // ------------------------------------------------------------ Environment
    ae_list<ae_gaussian*>* _env_gaussians;
    ae_list<Point*>* _env_custom_points;
    int16_t  _env_sampling;

    // ---------------------------------------- Environment x-axis segmentation
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

    // ------------------------------------------------------ Environment noise
    ae_env_noise _env_noise_method;   // Method... TODO
    double  _env_noise_alpha;         // Alpha value (variance coefficient)
    double  _env_noise_sigma;         // Variance of the noise
    double  _env_noise_prob;          // Probability of variation.
    int32_t _env_noise_sampling_log;  // Log2 of the number of points in the noise fuzzy_set

    // --------------------------------------------------------- Mutation rates
    double  _point_mutation_rate;
    double  _small_insertion_rate;
    double  _small_deletion_rate;
    int16_t _max_indel_size;

    // -------------------------------------------- Rearrangements and Transfer
    bool    _with_4pts_trans;
    bool    _with_alignments;
    bool    _with_HT;
    bool    _repl_HT_with_close_points;
    double  _HT_ins_rate;
    double  _HT_repl_rate;
    double  _repl_HT_detach_rate;

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

    // ----------------------------------------------- Phenotypic Stochasticity
    bool _with_stochasticity;

    // -------------------------------------------------------------- Selection
    ae_selection_scheme  _selection_scheme;
    double               _selection_pressure;

    // ------------------------------------------------------ Spatial structure
    bool     _spatially_structured;
    int16_t  _grid_width;
    int16_t  _grid_height;
    int32_t  _migration_number; // TODO : explain

    // -------------------------------------------------------------- Secretion
    bool   _with_secretion;
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
    int32_t   _plasmid_maximal_length;
    int32_t   _chromosome_minimal_length;
    int32_t   _chromosome_maximal_length;
    double    _prob_plasmid_HT;
    double    _tune_donor_ability;
    double    _tune_recipient_ability;
    double    _donor_cost;
    double    _recipient_cost;
    bool      _compute_phen_contrib_by_GU;
    bool      _swap_GUs;

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

} // namespace aevol
#endif // __param_loader_H__
