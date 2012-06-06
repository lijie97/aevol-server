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


/*! \class ae_params
    \brief Singleton. Contains all the parameters needed during the whole run.
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




// =================================================================
//                          Class declarations
// =================================================================






class ae_params : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_params( void );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_params( void );

    // =================================================================
    //                         Accessors: getters
    // =================================================================
    // ------------------------------------------------------------ Constraints
    inline int32_t  get_min_genome_length( void ) const;
    inline int32_t  get_max_genome_length( void ) const;
    inline double   get_min_w( void ) const;
    inline double   get_max_w( void ) const;

    // --------------------------------------------------------- Mutation rates
    inline double   get_point_mutation_rate( void ) const;
    inline double   get_small_insertion_rate( void ) const;
    inline double   get_small_deletion_rate( void ) const;
    inline int16_t  get_max_indel_size( void ) const;

    // -------------------------------------------- Rearrangements and Transfer
    inline bool   get_with_4pts_trans( void ) const;
    inline bool   get_with_alignments( void ) const;
    inline bool   get_with_transfer( void ) const;
    inline bool   get_swap_GUs( void ) const;
    inline double get_transfer_ins_rate( void ) const;
    inline double get_transfer_repl_rate( void ) const;

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
    
    inline int16_t get_align_max_shift( void ) const;     // Maximum shift of one seq on the other
    inline int16_t get_align_w_zone_h_len( void ) const;  // Work zone half length
    inline int16_t get_align_match_bonus( void ) const;   // Corresponding residues match bonus
    inline int16_t get_align_mismatch_cost( void ) const; // Corresponding residues mismatch cost

    // -------------------------------------------------------------- Selection
    inline ae_selection_scheme  get_selection_scheme( void ) const;
    inline double               get_selection_pressure( void ) const;
    
    // -------------------------------------------------------------- Secretion
    inline bool   get_use_secretion( void ) const;
    inline double get_secretion_fitness_contrib( void ) const;      // (0,1)
    // proportion that diffuses into each cell, every generation, 0 for no diffusion
    inline double get_secretion_diffusion_prop( void ) const;
    // proportion of secreted substance that degrades every generation generation
    inline double get_secretion_degradation_prop( void ) const;
    // cost of secreting the compound, as a proportion of the amount secreted
    inline double get_secretion_cost( void ) const;
    
    // --------------------------------------------------------------- Plasmids
    inline bool     get_allow_plasmids( void ) const;
    inline int32_t  get_plasmid_minimal_length( void ) const;
    inline double   get_prob_horiz_trans( void ) const;
    inline int16_t  get_nb_horiz_trans( void ) const;
    inline bool     get_compute_phen_contrib_by_GU( void ) const;
    
    // ------------------------------------------------------- Translation cost
    inline double get_translation_cost( void ) const;

    // =================================================================
    //                         Accessors: setters
    // =================================================================
    // ------------------------------------------------------------ Constraints
    inline void set_min_genome_length( int32_t min_genome_length );
    inline void set_max_genome_length( int32_t max_genome_length );
    inline void set_min_w( double min_w );
    inline void set_max_w( double max_w );

    // --------------------------------------------------------- Mutation rates
    inline void set_point_mutation_rate( double point_mutation_rate );
    inline void set_small_insertion_rate( double small_insertion_rate );
    inline void set_small_deletion_rate( double small_deletion_rate );
    inline void set_max_indel_size( int16_t max_indel_size );

    // -------------------------------------------- Rearrangements and Transfer
    inline void set_with_4pts_trans( bool with_4pts_trans );
    inline void set_with_alignments( bool with_alignments );
    inline void set_with_transfer( bool with_transfer );
    inline void set_swap_GUs( bool swap_GUs );
    inline void set_transfer_ins_rate( double transfer_ins_rate );
    inline void set_transfer_repl_rate( double transfer_repl_rate );

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
    
    // -------------------------------------------------------------- Secretion
    inline void set_use_secretion( bool use_secretion );
    inline void set_secretion_fitness_contrib( double secretion_fitness_contrib );
    inline void set_secretion_diffusion_prop( double secretion_diffusion_prop );
    inline void set_secretion_degradation_prop( double secretion_degradation_prop );
    inline void set_secretion_cost( double secretion_cost );
    
    // --------------------------------------------------------------- Plasmids
    inline void set_allow_plasmids( bool allow_plasmids );
    inline void set_plasmid_minimal_length( int32_t plasmid_minimal_length );
    inline void set_prob_horiz_trans( double prob_horiz_trans );
    inline void set_nb_horiz_trans( int16_t nb_horiz_trans );
    inline void set_compute_phen_contrib_by_GU( bool compute_phen_contrib_by_GU );
    
    // ------------------------------------------------------- Translation cost
    inline void set_translation_cost( double translation_cost );

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void write_to_backup( gzFile* backup_file );
    void read_from_backup( gzFile* backup_file, bool verbose = false );
    void print_to_file( FILE* file );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_params( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_params( const ae_params &model )
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
    // ------------------------------------------------------------ Constraints
    int32_t _min_genome_length;
    int32_t _max_genome_length;
    double  _min_w;
    double  _max_w;

    // --------------------------------------------------------- Mutation rates
    double  _point_mutation_rate;
    double  _small_insertion_rate;
    double  _small_deletion_rate;
    int16_t _max_indel_size;

    // -------------------------------------------- Rearrangements and Transfer
    bool    _with_4pts_trans;
    bool    _with_alignments;
    bool    _with_transfer;
    bool    _swap_GUs;
    double  _transfer_ins_rate;
    double  _transfer_repl_rate;

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
    
    // -------------------------------------------------------------- Secretion
    bool   _use_secretion;
    double _secretion_fitness_contrib;      // (0,1)
    // proportion that diffuses into each cell, every generation, 0 for no diffusion
    double _secretion_diffusion_prop;
    // proportion of secreted substance that degrades every generation generation
    double _secretion_degradation_prop;
    // cost of secreting the compound, as a proportion of the amount secreted
    double _secretion_cost;
    
    // --------------------------------------------------------------- Plasmids
    bool     _allow_plasmids;
    int32_t  _plasmid_minimal_length;
    double   _prob_horiz_trans;
    int16_t  _nb_horiz_trans; // How many plasmids can a cell send in its neighbourhood per generation
    bool     _compute_phen_contrib_by_GU;
    
    // ------------------------------------------------------- Translation cost
    double _translation_cost;
    
    // -------------------------------------------------------------- Log files
    int8_t _logs;
};


// =====================================================================
//                          Getters' definitions
// =====================================================================
// ------------------------------------------------------------ Constraints
inline int32_t ae_params::get_min_genome_length( void ) const
{
  return _min_genome_length;
}

inline int32_t ae_params::get_max_genome_length( void ) const
{
  return _max_genome_length;
}

inline double ae_params::get_min_w( void ) const
{
  return _min_w;
}

inline double ae_params::get_max_w( void ) const
{
  return _max_w;
}

// --------------------------------------------------------- Mutation rates
inline double ae_params::get_point_mutation_rate( void ) const
{
  return _point_mutation_rate;
}

inline double ae_params::get_small_insertion_rate( void ) const
{
  return _small_insertion_rate;
}

inline double ae_params::get_small_deletion_rate( void ) const
{
  return _small_deletion_rate;
}

inline int16_t ae_params::get_max_indel_size( void ) const
{
  return _max_indel_size;
}

// -------------------------------------------- Rearrangements and Transfer
inline bool ae_params::get_with_4pts_trans( void ) const
{
  return _with_4pts_trans;
}

inline bool ae_params::get_with_alignments( void ) const
{
  return _with_alignments;
}

inline bool ae_params::get_with_transfer( void ) const
{
  return _with_transfer;
}

inline bool ae_params::get_swap_GUs( void ) const
{
  return _swap_GUs;
}

inline double ae_params::get_transfer_ins_rate( void ) const
{
  return _transfer_ins_rate;
}

inline double ae_params::get_transfer_repl_rate( void ) const
{
  return _transfer_repl_rate;
}

// ------------------------------ Rearrangement rates (without alignements)
inline double ae_params::get_duplication_rate( void ) const
{
  return _duplication_rate;
}

inline double ae_params::get_deletion_rate( void ) const
{
  return _deletion_rate;
}

inline double ae_params::get_translocation_rate( void ) const
{
  return _translocation_rate;
}

inline double ae_params::get_inversion_rate( void ) const
{
  return _inversion_rate;
}

// --------------------------------- Rearrangement rates (with alignements)
inline double ae_params::get_neighbourhood_rate( void ) const
{
  return _neighbourhood_rate;
}

inline double ae_params::get_duplication_proportion( void ) const
{
  return _duplication_proportion;
}

inline double ae_params::get_deletion_proportion( void ) const
{
  return _deletion_proportion;
}

inline double ae_params::get_translocation_proportion( void ) const
{
  return _translocation_proportion;
}

inline double ae_params::get_inversion_proportion( void ) const
{
  return _inversion_proportion;
}

// ------------------------------------------------------------ Alignements
inline ae_align_fun_shape ae_params::get_align_fun_shape( void ) const
{
  return _align_fun_shape;
}

inline double ae_params::get_align_sigm_lambda( void ) const
{
  return _align_sigm_lambda;
}

inline int16_t ae_params::get_align_sigm_mean( void ) const
{
  return _align_sigm_mean;
}

inline int16_t ae_params::get_align_lin_min( void ) const
{
  return _align_lin_min;
}

inline int16_t ae_params::get_align_lin_max( void ) const
{
  return _align_lin_max;
}

inline int16_t ae_params::get_align_max_shift( void ) const
{
  return _align_max_shift;
}

inline int16_t ae_params::get_align_w_zone_h_len( void ) const
{
  return _align_w_zone_h_len;
}

inline int16_t ae_params::get_align_match_bonus( void ) const
{
  return _align_match_bonus;
}

inline int16_t ae_params::get_align_mismatch_cost( void ) const
{
  return _align_mismatch_cost;
}

// -------------------------------------------------------------- Selection
inline ae_selection_scheme ae_params::get_selection_scheme( void ) const
{
  return _selection_scheme;
}

inline double ae_params::get_selection_pressure( void ) const
{
  return _selection_pressure;
}

// -------------------------------------------------------------- Secretion
inline bool ae_params::get_use_secretion( void ) const
{
  return _use_secretion;
}

inline double ae_params::get_secretion_fitness_contrib( void ) const
{
  return _secretion_fitness_contrib;
}

inline double ae_params::get_secretion_diffusion_prop( void ) const
{
  return _secretion_diffusion_prop;
}

inline double ae_params::get_secretion_degradation_prop( void ) const
{
  return _secretion_degradation_prop;
}

inline double ae_params::get_secretion_cost( void ) const
{
  return _secretion_cost;
}

// --------------------------------------------------------------- Plasmids
inline bool ae_params::get_allow_plasmids( void ) const
{
  return _allow_plasmids;
}

inline int32_t ae_params::get_plasmid_minimal_length( void ) const
{
  return _plasmid_minimal_length;
}

inline double ae_params::get_prob_horiz_trans( void ) const
{
  return _prob_horiz_trans;
}

inline int16_t ae_params::get_nb_horiz_trans( void ) const
{
  return _nb_horiz_trans;
}

inline bool ae_params::get_compute_phen_contrib_by_GU( void ) const
{
  return _compute_phen_contrib_by_GU;
}

// ------------------------------------------------------- Translation cost
inline double ae_params::get_translation_cost( void ) const
{
  return _translation_cost;
}

// =====================================================================
//                          Setters' definitions
// =====================================================================
// ------------------------------------------------------------ Constraints
inline void ae_params::set_min_genome_length( int32_t min_genome_length )
{
  _min_genome_length = min_genome_length;
}

inline void ae_params::set_max_genome_length( int32_t max_genome_length )
{
  _max_genome_length = max_genome_length;
}

inline void ae_params::set_min_w( double min_w )
{
  _min_w =  min_w;
}

inline void ae_params::set_max_w( double max_w )
{
  _max_w = max_w;
}

// --------------------------------------------------------- Mutation rates
inline void ae_params::set_point_mutation_rate( double point_mutation_rate )
{
  _point_mutation_rate = point_mutation_rate;
}

inline void ae_params::set_small_insertion_rate( double small_insertion_rate )
{
  _small_insertion_rate = small_insertion_rate;
}

inline void ae_params::set_small_deletion_rate( double small_deletion_rate )
{
  _small_deletion_rate = small_deletion_rate;
}

inline void ae_params::set_max_indel_size( int16_t max_indel_size )
{
  _max_indel_size = max_indel_size;
}

// -------------------------------------------- Rearrangements and Transfer
inline void ae_params::set_with_4pts_trans( bool with_4pts_trans )
{
  _with_4pts_trans = with_4pts_trans;
}

inline void ae_params::set_with_alignments( bool with_alignments )
{
  _with_alignments = with_alignments;
}

inline void ae_params::set_with_transfer( bool with_transfer )
{
  _with_transfer = with_transfer;
}

inline void ae_params::set_swap_GUs( bool swap_GUs )
{
  _swap_GUs = swap_GUs;
}

inline void ae_params::set_transfer_ins_rate( double transfer_ins_rate )
{
  _transfer_ins_rate = transfer_ins_rate;
}

inline void ae_params::set_transfer_repl_rate( double transfer_repl_rate )
{
  _transfer_repl_rate = transfer_repl_rate;
}

// ------------------------------ Rearrangement rates (without alignements)
inline void ae_params::set_duplication_rate( double duplication_rate )
{
  _duplication_rate = duplication_rate;
}

inline void ae_params::set_deletion_rate( double deletion_rate )
{
  _deletion_rate = deletion_rate;
}

inline void ae_params::set_translocation_rate( double translocation_rate )
{
  _translocation_rate = translocation_rate;
}

inline void ae_params::set_inversion_rate( double inversion_rate )
{
  _inversion_rate = inversion_rate;
}

// --------------------------------- Rearrangement rates (with alignements)
inline void ae_params::set_neighbourhood_rate( double neighbourhood_rate )
{
  _neighbourhood_rate = neighbourhood_rate;
}

inline void ae_params::set_duplication_proportion( double duplication_proportion )
{
  _duplication_proportion = duplication_proportion;
}

inline void ae_params::set_deletion_proportion( double deletion_proportion )
{
  _deletion_proportion = deletion_proportion;
}

inline void ae_params::set_translocation_proportion( double translocation_proportion )
{
  _translocation_proportion = translocation_proportion;
}

inline void ae_params::set_inversion_proportion( double inversion_proportion )
{
  _inversion_proportion = inversion_proportion;
}

// ------------------------------------------------------------ Alignements
inline void ae_params::set_align_fun_shape( ae_align_fun_shape align_fun_shape )
{
  _align_fun_shape = align_fun_shape;
}

inline void ae_params::set_align_sigm_lambda( double align_sigm_lambda )
{
  _align_sigm_lambda = align_sigm_lambda;
}

inline void ae_params::set_align_sigm_mean( int16_t align_sigm_mean )
{
  _align_sigm_mean = align_sigm_mean;
}

inline void ae_params::set_align_lin_min( int16_t align_lin_min )
{
  _align_lin_min = align_lin_min;
}

inline void ae_params::set_align_lin_max( int16_t align_lin_max )
{
  _align_lin_max = align_lin_max;
}

inline void ae_params::set_align_max_shift( int16_t align_max_shift )
{
  _align_max_shift = align_max_shift;
}

inline void ae_params::set_align_w_zone_h_len( int16_t align_w_zone_h_len )
{
  _align_w_zone_h_len = align_w_zone_h_len;
}

inline void ae_params::set_align_match_bonus( int16_t align_match_bonus )
{
  _align_match_bonus = align_match_bonus;
}

inline void ae_params::set_align_mismatch_cost( int16_t align_mismatch_cost )
{
  _align_mismatch_cost = align_mismatch_cost;
}

// -------------------------------------------------------------- Selection
inline void ae_params::set_selection_scheme( ae_selection_scheme selection_scheme )
{
  _selection_scheme = selection_scheme;
}

inline void ae_params::set_selection_pressure( double selection_pressure )
{
  _selection_pressure = selection_pressure;
}

// -------------------------------------------------------------- Secretion
inline void ae_params::set_use_secretion( bool use_secretion )
{
  _use_secretion = use_secretion;
}

inline void ae_params::set_secretion_fitness_contrib( double secretion_fitness_contrib )
{
  _secretion_fitness_contrib = secretion_fitness_contrib;
}

inline void ae_params::set_secretion_diffusion_prop( double secretion_diffusion_prop )
{
  _secretion_diffusion_prop = secretion_diffusion_prop;
}

inline void ae_params::set_secretion_degradation_prop( double secretion_degradation_prop )
{
  _secretion_degradation_prop = secretion_degradation_prop;
}

inline void ae_params::set_secretion_cost( double secretion_cost )
{
  _secretion_cost = secretion_cost;
}

// --------------------------------------------------------------- Plasmids
inline void ae_params::set_allow_plasmids( bool allow_plasmids )
{
  _allow_plasmids = allow_plasmids;
}

inline void ae_params::set_plasmid_minimal_length( int32_t plasmid_minimal_length )
{
  _plasmid_minimal_length = plasmid_minimal_length;
}

inline void ae_params::set_prob_horiz_trans( double prob_horiz_trans )
{
  _prob_horiz_trans = prob_horiz_trans;
}

inline void ae_params::set_nb_horiz_trans( int16_t nb_horiz_trans )
{
  _nb_horiz_trans = nb_horiz_trans;
}

inline void ae_params::set_compute_phen_contrib_by_GU( bool compute_phen_contrib_by_GU )
{
  _compute_phen_contrib_by_GU = compute_phen_contrib_by_GU;
}

// ------------------------------------------------------- Translation cost
inline void ae_params::set_translation_cost( double translation_cost )
{
  _translation_cost = translation_cost;
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_PARAMS_H__
