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
// ****************************************************************************


#ifndef __AE_INDIVIDUAL_H__
#define __AE_INDIVIDUAL_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include <environment.h>
#include "genetic_unit.h"
#include "ae_dna.h"
#include "ae_rna.h"
#include "ae_protein.h"
#include "ae_phenotype.h"
#include "ae_replication_report.h"
#include "ae_params_mut.h"
#include "fuzzy.h"


namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;
class ae_grid_cell;


/// Models an individual cell.
///
/// Proteins ans RNAs are shared with descent. Genetic units are an individual's own.
class ae_individual
{
  friend class ae_dna;

  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_individual(  ae_exp_manager* exp_m,
                    ae_jumping_mt* mut_prng,
                    ae_jumping_mt* stoch_prng,
                    ae_params_mut* param_mut,
                    double w_max,
                    int32_t min_genome_length,
                    int32_t max_genome_length,
                    bool allow_plasmids,
                    int32_t id,
                    char* strain_name,
                    int32_t age );

    ae_individual( ae_exp_manager* exp_m, gzFile backup_file );
    ae_individual(const ae_individual &model, bool replication_report_copy = false );
    ae_individual(  ae_individual* const parent, int32_t id,
                    ae_jumping_mt* mut_prng, ae_jumping_mt* stoch_prng );
    ae_individual() = delete; // forbidden constructor


    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_individual( void );

    // =================================================================
    //                        Accessors: Getters
    // =================================================================
    const char*  get_strain_name(void) const;
    int32_t      get_age(void) const;

    ae_exp_manager*  get_exp_m( void ) const;
    int16_t          get_nb_genetic_units( void ) const;
    int32_t          get_nb_plasmids( void )  const;
    int32_t          get_amount_of_dna( void ) const;
    const GeneticUnit& get_genetic_unit( int16_t num_unit ) const;
    GeneticUnit& get_genetic_unit_nonconst( int16_t num_unit );
    double           get_dist_to_target_by_feature( ae_env_axis_feature feature ) const;
    double           get_fitness( void ) const;
    double           get_fitness_by_feature( ae_env_axis_feature feature ) const;
    ae_grid_cell*    get_grid_cell( void ) const;
    bool             get_placed_in_population( void ) const;
    ae_jumping_mt*   get_mut_prng( void ) const;
    ae_jumping_mt*   get_stoch_prng( void ) const;

    const std::list<GeneticUnit>& get_genetic_unit_list_std() const;
    std::list<GeneticUnit>& get_genetic_unit_list_std_nonconst();
    /// Keep only the first (main chromosome) and the last chromosomes from the GU.
    void drop_nested_genetic_units();

    const char* get_genetic_unit_sequence   ( int16_t num_unit ) const;
    int32_t     get_genetic_unit_seq_length ( int16_t num_unit ) const;

    int32_t get_id( void ) const;
    int32_t get_rank( void ) const;


    Fuzzy* get_phenotype_activ( void ) const;
    Fuzzy* get_phenotype_inhib( void ) const;
    ae_phenotype* get_phenotype( void )       const;

    double* get_dist_to_target_by_segment( void ) const;

    ae_replication_report* get_replic_report( void ) const;

    const std::list<ae_protein*>&  get_protein_list( void )  const;
    const std::list<ae_rna*>&      get_rna_list( void )      const;

    double get_w_max( void ) const;

    // Genome size constraints
    int32_t get_min_genome_length( void ) const;
    int32_t get_max_genome_length( void ) const;

    // Plasmids settings
    bool     get_allow_plasmids( void ) const;

    // Mutation rates etc...
    double   get_point_mutation_rate( void ) const;
    double   get_small_insertion_rate( void ) const;
    double   get_small_deletion_rate( void ) const;
    int16_t  get_max_indel_size( void ) const;

    double get_duplication_rate( void ) const;
    double get_deletion_rate( void ) const;
    double get_translocation_rate( void ) const;
    double get_inversion_rate( void ) const;
    double get_neighbourhood_rate( void ) const;
    double get_duplication_proportion( void ) const;
    double get_deletion_proportion( void ) const;
    double get_translocation_proportion( void ) const;
    double get_inversion_proportion( void ) const;

    // Transfer
    bool   get_with_4pts_trans( void ) const;
    bool   get_with_HT( void ) const;
    bool   get_repl_HT_with_close_points( void ) const;
    double get_HT_ins_rate( void ) const;
    double get_HT_repl_rate( void ) const;
    double get_repl_HT_detach_rate( void ) const;

    // ------------------------------------------------------------ Alignements
    bool  get_with_alignments(void) const;
    ae_align_fun_shape get_align_fun_shape( void ) const;
    double  get_align_sigm_lambda( void ) const;
    int16_t get_align_sigm_mean( void ) const;
    int16_t get_align_lin_min( void ) const;
    int16_t get_align_lin_max( void ) const;

    // Maximum shift of one seq on the other
    int16_t get_align_max_shift( void ) const;
    // Work zone half length
    int16_t get_align_w_zone_h_len( void ) const;
    // Corresponding residues match bonus
    int16_t get_align_match_bonus( void ) const;
    // Corresponding residues mismatch cost
    int16_t get_align_mismatch_cost( void ) const;


    // ----------------------------------------------- Phenotypic stochasticity
    bool get_with_stochasticity( void ) const;

    // Statistical data
    int32_t  get_total_genome_size( void ); // TODO: duplicate with get_amount_of_dna?
    int32_t  get_nb_coding_RNAs( void );
    int32_t  get_nb_non_coding_RNAs( void );
    double   get_overall_size_coding_RNAs( void );
    double   get_av_size_coding_RNAs( void );
    double   get_overall_size_non_coding_RNAs( void );
    double   get_av_size_non_coding_RNAs( void );
    int32_t  get_nb_genes_activ( void );
    int32_t  get_nb_genes_inhib( void );
    int32_t  get_nb_functional_genes( void );
    int32_t  get_nb_non_functional_genes( void );
    double   get_overall_size_functional_genes( void );
    double   get_av_size_functional_genes( void );
    double   get_overall_size_non_functional_genes( void );
    double   get_av_size_non_functional_genes( void );

    int32_t  get_nb_bases_in_0_CDS( void );
    int32_t  get_nb_bases_in_0_functional_CDS( void );
    int32_t  get_nb_bases_in_0_non_functional_CDS( void );
    int32_t  get_nb_bases_in_0_RNA( void );
    int32_t  get_nb_bases_in_0_coding_RNA( void );
    int32_t  get_nb_bases_in_0_non_coding_RNA( void );
    int32_t  get_nb_bases_in_neutral_regions( void );
    int32_t  get_nb_neutral_regions( void );

    double get_modularity( void ); // Not const

    int32_t* get_int_probes ( void ) const;
    double* get_double_probes ( void) const;


    // =================================================================
    //                        Accessors: Setters
    // =================================================================
    void set_strain_name(char* name);

    void set_exp_m( ae_exp_manager* exp_m );
    void set_id( int32_t id );
    void set_rank( int32_t rank );

    void set_grid_cell( ae_grid_cell* grid_cell );
    void set_placed_in_population( bool placed_in_population );

    void set_replication_report( ae_replication_report * rep );

    void reset_dist_to_target_by_segment( double * dist_to_target_by_segment );

    void set_w_max( double w_max );

    // Genome size constraints
    void set_min_genome_length( int32_t min_genome_length );
    void set_max_genome_length( int32_t max_genome_length );

    // Plasmids
    void set_allow_plasmids( bool allow_plasmids );

    // Mutation rates etc...
    void set_point_mutation_rate( double point_mutation_rate);
    void set_small_insertion_rate( double small_insertion_rate);
    void set_small_deletion_rate( double small_deletion_rate);
    void set_max_indel_size( int16_t max_indel_size);

    void set_duplication_rate( double duplication_rate);
    void set_deletion_rate( double deletion_rate);
    void set_translocation_rate( double translocation_rate);
    void set_inversion_rate( double inversion_rate);
    void set_neighbourhood_rate( double neighbourhood_rate);
    void set_duplication_proportion( double duplication_proportion);
    void set_deletion_proportion( double deletion_proportion);
    void set_translocation_proportion( double translocation_proportion);
    void set_inversion_proportion( double inversion_proportion);

    // Transfer
    void set_with_4pts_trans( bool with_4pts_trans );
    void set_with_alignments( bool with_alignments );
    void set_with_HT( bool with_HT );
    void set_repl_HT_with_close_points( bool repl_HT_with_close_points);
    void set_HT_ins_rate( double HT_ins_rate );
    void set_HT_repl_rate( double HT_repl_rate );
    void set_repl_HT_detach_rate( double repl_HT_detach_rate );

    // ----------------------------------------------- Phenotypic stochasticity
    void set_with_stochasticity( bool with_stoch );

    void set_mut_prng( ae_jumping_mt* prng );
    void set_stoch_prng( ae_jumping_mt* prng );

    //------------------------------------------------ Generic probes
    void set_int_probes ( int32_t* int_probes);
    void set_double_probes ( double* double_probes);


    // =================================================================
    //                            Public Methods
    // =================================================================


    void add_GU( char * &sequence, int32_t length );  // warning: the individual is left in a totally "cleared" state but not reevaluated
    // void add_GU(ae_genetic_unit&& unit ); // warning: the individual is left in a totally "cleared" state but not reevaluated
  void add_GU(ae_individual* indiv, int32_t chromosome_length, ae_jumping_mt* prng);

    // void add_GU(ae_individual* indiv, int32_t length, ae_jumping_mt* prng);
    void remove_GU ( int16_t num_unit ); // warning: the individual is left in a totally "cleared" state but not reevaluated
    void renew_dist_to_target_by_feature( void );
    void renew_fitness_by_feature( void );

    void inject_GU( ae_individual* donor );
    void inject_2GUs( ae_individual* partner );

    virtual void evaluate( Environment* envir = NULL );
    virtual void reevaluate( Environment* envir = NULL );
    virtual void clear_everything_except_dna_and_promoters(void);
    void do_transcription_translation_folding( void );
    void do_transcription( void );
    void do_translation( void );
    void do_folding( void );
    void compute_phenotype( void );
    void compute_distance_to_target( Environment* envir );

    // Computation of a "proper" fitness value (one that increases when the individual is fitter)
    // The behaviour of this function depends on many parameters and most notably on whether it is
    // a "composite" fitness or not, and on the selection scheme.
    void compute_fitness( Environment* envir );

    void compute_statistical_data();
    void compute_non_coding( void );

    virtual void save( gzFile backup_file ) const;

    int32_t get_nb_terminators( void );

    #ifdef DEBUG
      void assert_promoters( void );
      void assert_promoters_order( void );
    #endif

    void compute_experimental_f_nu( int32_t nb_children, double* reproduction_statistics, double* offsprings_statistics = NULL, FILE* replication_file = NULL);
    double compute_theoritical_f_nu( void );
    // These functions compute the probability of neutral reproduction (F_nu).
    // The first method replicates the individual "nb_children" times and counts how often
    // a child has the same fitness as its parent (and if second argument not NULL, how often
    // a child has same or better fitness and stores the result at the adress contained in the
    // pointer). Results are proportions.
    // The second is an estimate based on genome structure as defined by Carole.
    // They have been implemented on the chromosome only !

    void remove_non_coding_bases(void);
    void double_non_coding_bases(void);


    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================
    virtual void make_protein_list( void );
    virtual void make_rna_list( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;

    // Name and "Age" of the strain
    char*   _strain_name;
    int32_t _age;

    // Random number generator
    ae_jumping_mt* _mut_prng;
    ae_jumping_mt* _stoch_prng;

    // Individual ID and rank of the individual in the population
    // WARNING : The ID is no longer corresponding to the rank of the individual.
    //           The reason for this change is that we now need an identifier for the individuals
    //           as soon as they are created (the rank is only known when all the individuals have been evaluated).
    //           The rank will now be handled in a specific new attribute. (1 for the worst indiv, POP_SIZE for the best)
    int32_t  _id;   // [0 ; POP_SIZE[
    int32_t  _rank; // [1 ; POP_SIZE]

    // Total activation (resp. inhibition) of metabolic functions
    Fuzzy* _phenotype_activ;
    Fuzzy* _phenotype_inhib;

    // The phenotype, roughly corresponding to the sum of activ and inhib
    ae_phenotype* _phenotype;

    // Table containing the partial area of the "gap" (difference between the phenotype and the environment)
    // for each environmental segment.
    // Note: if the environment is not segmented, there will be a single cell in this table
    double* _dist_to_target_by_segment;

    // This table contains the aggregated area of the gap for each type of segment (ae_env_axis_feature).
    // When the environment is not segmented, the only meaningfull value will be held in _dist_to_target[METABOLISM]
    double* _dist_to_target_by_feature;

    // This table contains the fitness contribution for each type of segment (ae_env_axis_feature).
    // For metabolic segment, this is the metabolic fitness contribution;
    // For secretion segment, this is the amount of compound that gets secreted into the environment;
    // Note: total fitness is the combination of metabolic fitness and the amount of compound present in the environment,
    //       not the amount of compound secreted by the individual.
    double* _fitness_by_feature;

    // THE fitness
    double _fitness;

    // When using structured population, this is the cell the individual is in
    ae_grid_cell* _grid_cell;

    // The chromosome and plasmids (if allowed)
    std::list<GeneticUnit> _genetic_unit_list;

    // Report of all the mutational events undergone during the individuals creation,
    // i.e. during the replication that gave birth to this individual
    ae_replication_report* _replic_report;

    // Access lists to all the proteins/RNAs of the individual.
    // Please note that these proteins/RNAs are actually managed (i.e. newed and deleted) via genetic units.
    std::list<ae_protein*> _protein_list;
    std::list<ae_rna*>     _rna_list;

    // Generic probes
    int32_t*  _int_probes;        // Table of 5 int32_t values to be used as one wishes
    double*   _double_probes;     // Table of 5 double values to be used as one wishes

    // Mutation rates etc...
    ae_params_mut* _mut_params;

    // ----------------------------------------------- Phenotypic stochasticity
    bool _with_stochasticity;

    // Artificial chemistry
    double _w_max;

    // Genome size constraints
    int32_t _min_genome_length;
    int32_t _max_genome_length;

    // Plasmids settings
    bool _allow_plasmids;



    // --------------------------------------------------
    // "State" of the individual
    // --------------------------------------------------
    // We keep trace of what we have already computed to avoid double computation (mainly in post-treaments)
    bool _evaluated;
    bool _transcribed;
    bool _translated;
    bool _folded;
    bool _phenotype_computed;
    bool _distance_to_target_computed;
    bool _fitness_computed;
    bool _statistical_data_computed;
    bool _non_coding_computed;
    bool _modularity_computed;

    bool _placed_in_population; // TODO: spatial ?



    // ----------------------------------------
    // Statistical data
    // ----------------------------------------
    // Genome, RNAs and genes size and stuff
    int32_t _total_genome_size;                 // Sum of sizes of the genetic units
    int32_t _nb_coding_RNAs;                    // Number of coding RNAs (at least one gene on RNA)
    int32_t _nb_non_coding_RNAs;                // Number of non-coding-RNAs
    double  _overall_size_coding_RNAs;          // Average size of coding RNAs
    double  _overall_size_non_coding_RNAs;      // Average size of non-coding RNAs
    int32_t _nb_genes_activ;                    // Number of genes realizing a function
    int32_t _nb_genes_inhib;                    // Number of genes inhibitting a function
    int32_t _nb_functional_genes;               // Number of functional genes
    int32_t _nb_non_functional_genes;           // Number of non-functional genes
    double  _overall_size_functional_genes;     // Average size of functional genes
    double  _overall_size_non_functional_genes; // Average size of non-functional genes

    // Mutation/Rearrangement statistics are managed in the replication report

    // Coding / non-coding
    int32_t _nb_bases_in_0_CDS;               // Number of bases that are not included in any gene
    int32_t _nb_bases_in_0_functional_CDS;    // Number of bases that are not included in any functional gene
    int32_t _nb_bases_in_0_non_functional_CDS;// Number of bases that are not included in any degenerated gene
    int32_t _nb_bases_in_0_RNA;               // Number of bases that are not included in any RNA
    int32_t _nb_bases_in_0_coding_RNA;        // Number of bases that are not included in any coding RNA
                                              // (RNAs containing at least one CDS)
    int32_t _nb_bases_in_0_non_coding_RNA;    // Number of bases that are not included in any non coding RNA
    int32_t _nb_bases_in_neutral_regions;     // Number of bases that are in a neutral region
                                              // A base is considered neutral when neither itself NOR its corresponding base on the other
                                              // strand belongs to a coding promoter->terminator region (both included)
    int32_t _nb_neutral_regions;              // Number of neutral regions

    double _modularity; // Ratio between the pairwise distance between genes whose corresponding
                        // phenotypic triangles overlap and the average intergenic distance
                        // (ignoring non-functional genes)void compute_phenotype( void );
};

} // namespace aevol
#endif // __AE_INDIVIDUAL_H__
