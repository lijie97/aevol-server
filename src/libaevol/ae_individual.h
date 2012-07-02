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


/** \class
 *  \brief
 */


#ifndef __AE_INDIVIDUAL_H__
#define __AE_INDIVIDUAL_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_list.h>
#include <ae_environment.h>
#include <ae_genetic_unit.h>
#include <ae_dna.h>
#include <ae_rna.h>
#include <ae_protein.h>
#include <ae_phenotype.h>
#include <ae_replication_report.h>
#include <ae_params_mut.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;
class ae_grid_cell;






class ae_individual : public ae_object
{
  friend class ae_dna;
  
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_individual( ae_exp_manager* exp_m, ae_rand_mt* alea, ae_params_mut* param_mut, int32_t id, int32_t age );
    ae_individual( ae_exp_manager* exp_m, gzFile* backup_file );
    ae_individual( const ae_individual &model );
    //~ ae_individual( ae_individual* const parent, int32_t index );
    //~ ae_individual( char* organism_file_name );
    //~ ae_individual( char* genome, int32_t genome_size, int32_t age = 0 );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_individual( void );

    // =================================================================
    //                        Accessors: Getters
    // =================================================================
    inline int16_t          get_nb_genetic_units( void )                                  const;
    inline int32_t          get_amount_of_dna( void )                                     const;
    inline ae_list*         get_genetic_unit_list( void )                                 const;
    inline ae_genetic_unit* get_genetic_unit( int16_t num_unit )                          const;
    inline double           get_dist_to_target_by_feature( ae_env_axis_feature feature )  const;
    inline double           get_fitness( void )                                           const;
    inline double           get_fitness_by_feature( ae_env_axis_feature feature )         const;
    inline ae_grid_cell*    get_grid_cell( void )                                         const;
    inline bool             get_placed_in_population()                                    const;
    
    inline const char* get_genetic_unit_sequence   ( int16_t num_unit ) const;
    inline int32_t     get_genetic_unit_seq_length ( int16_t num_unit ) const;
    
    inline int32_t get_id( void ) const;
    inline int32_t get_rank( void ) const;
    
    
    inline ae_fuzzy_set* get_phenotype_activ( void ) const;
    inline ae_fuzzy_set* get_phenotype_inhib( void ) const;
    inline ae_phenotype* get_phenotype( void )       const;
    
    inline double* get_dist_to_target_by_segment( void ) const;
    
    inline ae_replication_report* get_replic_report( void ) const;
    
    inline ae_list* get_protein_list( void )  const;
    inline ae_list* get_rna_list( void )      const;
    
    inline int32_t get_nb_plasmids( void )  const;
    inline int32_t get_nb_gen_units( void ) const;
    
    inline double get_w_max( void ) const;
    
    // Mutation rates etc...
    inline double   get_point_mutation_rate( void ) const;
    inline double   get_small_insertion_rate( void ) const;
    inline double   get_small_deletion_rate( void ) const;
    inline int16_t  get_max_indel_size( void ) const;
    
    inline double get_duplication_rate( void ) const;
    inline double get_deletion_rate( void ) const;
    inline double get_translocation_rate( void ) const;
    inline double get_inversion_rate( void ) const;
    inline double get_neighbourhood_rate( void ) const;
    inline double get_duplication_proportion( void ) const;
    inline double get_deletion_proportion( void ) const;
    inline double get_translocation_proportion( void ) const;
    inline double get_inversion_proportion( void ) const;
    
    // Statistical data
    inline int32_t  get_total_genome_size( void ); // TODO: duplicate with get_amount_of_dna?
    inline int32_t  get_nb_coding_RNAs( void );
    inline int32_t  get_nb_non_coding_RNAs( void );
    inline double   get_overall_size_coding_RNAs( void );
    inline double   get_av_size_coding_RNAs( void );
    inline double   get_overall_size_non_coding_RNAs( void );
    inline double   get_av_size_non_coding_RNAs( void );
    inline int32_t  get_nb_genes_activ( void );
    inline int32_t  get_nb_genes_inhib( void );
    inline int32_t  get_nb_functional_genes( void );
    inline int32_t  get_nb_non_functional_genes( void );
    inline double   get_overall_size_functional_genes( void );
    inline double   get_av_size_functional_genes( void );
    inline double   get_overall_size_non_functional_genes( void );
    inline double   get_av_size_non_functional_genes( void );
    
    inline int32_t  get_nb_bases_in_0_CDS( void );
    inline int32_t  get_nb_bases_in_0_functional_CDS( void );
    inline int32_t  get_nb_bases_in_0_non_functional_CDS( void );
    inline int32_t  get_nb_bases_in_0_RNA( void );
    inline int32_t  get_nb_bases_in_0_coding_RNA( void );
    inline int32_t  get_nb_bases_in_0_non_coding_RNA( void );
    inline int32_t  get_nb_bases_in_neutral_regions( void );
    inline int32_t  get_nb_neutral_regions( void );
    
    inline double get_modularity( void ); // Not const
    

    // =================================================================
    //                        Accessors: Setters
    // =================================================================
    inline void set_id( int32_t id );
    inline void set_rank( int32_t rank );
    
    inline void set_grid_cell( ae_grid_cell* grid_cell );
    inline void set_placed_in_population( bool placed_in_population );
    
    inline void set_replication_report( ae_replication_report * rep );
    
    inline void reset_dist_to_target_by_segment( double * dist_to_target_by_segment );
    
    inline void set_w_max( double w_max );
    
    // =================================================================
    //                            Public Methods
    // =================================================================
    inline void renew_dist_to_target_by_feature( void );
    inline void renew_fitness_by_feature( void );
    
    void inject_GU( ae_individual* donor );
    void inject_2GUs( ae_individual* partner );
    
    virtual void evaluate( ae_environment* envir = NULL );
    virtual void reevaluate( ae_environment* envir = NULL );
    inline void do_transcription_translation_folding( void );
    inline void do_transcription( void );
    inline void do_translation( void );
    inline void do_folding( void );
    void compute_phenotype( void );
    void compute_distance_to_target( ae_environment* envir );
    
    // Computation of a "proper" fitness value (one that increases when the individual is fitter)
    // The behaviour of this function depends on many parameters and most notably on whether it is
    // a "composite" fitness or not, and on the selection scheme.
    void compute_fitness( ae_environment* envir );
    
    void compute_statistical_data();
    void compute_non_coding( void );
    /*void compute_modularity( void );*/

    virtual void write_to_backup( gzFile* backup_file ) const;
    
    int32_t get_nb_terminators( void );
    
    #ifdef DEBUG
      inline void print_protein_list( void );
    
      inline void assert_promoters( void );
      inline void assert_promoters_order( void );
    #endif

    double compute_experimental_f_nu( int32_t nb_children, double* neutral_or_better = NULL );
    double compute_theoritical_f_nu( void );
    // These functions compute the probability of neutral reproduction (F_nu).
    // The first method replicates the individual "nb_children" times and counts how often
    // a child has the same fitness as its parent (and if second argument not NULL, how often
    // a child has same or better fitness and stores the result at the adress contained in the
    // pointer). Results are proportions.
    // The second is an estimate based on genome structure as defined by Carole.
    // They have been implemented on the chromosome only !
    

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_individual( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    /*ae_individual( const ae_individual &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/

    // =================================================================
    //                           Protected Methods
    // =================================================================
    virtual void make_protein_list( void );
    virtual void make_rna_list( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;
    
    // "Age" of the strain
    int32_t _age;
    
    // Random number generator
    ae_rand_mt* _alea;
  
    // Individual ID and rank of the individual in the population
    // WARNING : The ID is no longer corresponding to the rank of the individual.
    //           The reason for this change is that we now need an identifier for the individuals
    //           as soon as they are created (the rank is only known when all the individuals have been evaluated).
    //           The rank will now be handled in a specific new attribute. (1 for the worst indiv, POP_SIZE for the best)
    int32_t  _id;   // [0 ; POP_SIZE[
    int32_t  _rank; // [1 ; POP_SIZE]
    
    // Total activation (resp. inhibition) of metabolic functions
    ae_fuzzy_set* _phenotype_activ;
    ae_fuzzy_set* _phenotype_inhib;
    
    // The phenotype, roughly corresponding to the sum of activ and inhib
    ae_phenotype* _phenotype;
    
    // Number of segments the phenotype is to be split into (corresponding to the number of segments in the environment
    int16_t _nb_segments;
    
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
    ae_list* _genetic_unit_list;
    
    // Report of all the mutational events undergone during the individuals creation,
    // i.e. during the replication that gave birth to this individual
    ae_replication_report* _replic_report;
    
    // Access lists to all the proteins/RNAs of the individual.
    // Please note that these proteins/RNAs are actually managed (i.e. newed and deleted) via genetic units.
    ae_list*      _protein_list;
    ae_list*      _rna_list;
    
    // Generic probes
    int32_t*  _int_probes;        // Table of 5 int32_t values to be used as one wishes
    double*   _double_probes;     // Table of 5 double values to be used as one wishes
    
    
    // Mutation rates etc...
    ae_params_mut* _mut_params;
    
    // Artificial chemistry
    double _w_max;
    
    
    
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
    
    bool _placed_in_population;
    
    
    
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


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline int32_t ae_individual::get_id( void ) const
{
  return _id;
}

inline void ae_individual::set_id( int32_t id )
{
  _id = id;
  
  if ( _replic_report != NULL )
  {
    _replic_report->set_id( id );
  }
}

inline double* ae_individual::get_dist_to_target_by_segment( void ) const
{
  return _dist_to_target_by_segment;
}

inline void ae_individual::reset_dist_to_target_by_segment( double* dist_to_target_by_segment )
{
  if ( _dist_to_target_by_segment != NULL ) delete [] _dist_to_target_by_segment;
  _dist_to_target_by_segment = dist_to_target_by_segment;
}

/*!
  Get the individual's rank in the population (1 for the worst indiv, POP_SIZE for the best)

  Warning: be sure you call sort_individuals() before using get_rank_in_population
*/
inline int32_t ae_individual::get_rank( void ) const
{
  return _rank;
}

inline void ae_individual::set_rank( int32_t rank )
{
  _rank = rank;
  
  if ( _replic_report != NULL )
  {
    _replic_report->set_rank( rank );
  }
}

inline void ae_individual::set_grid_cell( ae_grid_cell* grid_cell )
{
  _grid_cell = grid_cell;
  _placed_in_population = true;
}

inline void ae_individual::set_placed_in_population( bool placed_in_population )
{
  _placed_in_population = placed_in_population;
}

/*!
  Set the individual's replication report
  To be used by post-treatment only
*/
inline void ae_individual::set_replication_report( ae_replication_report * rep )
{
  _replic_report = rep;
}

/*!
  Returns the number of genetic units
*/
inline int16_t ae_individual::get_nb_genetic_units( void ) const
{
  return _genetic_unit_list->get_nb_elts();
}

inline int32_t ae_individual::get_amount_of_dna( void ) const
{
  int32_t amount = 0;
  
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
    
    amount += gen_unit->get_dna()->get_length();
    
    gen_unit_node = gen_unit_node->get_next();
  }
  
  return amount;
}

/*!
  Returns the list of genetic units
*/
inline ae_list* ae_individual::get_genetic_unit_list( void ) const
{
  return _genetic_unit_list;
}

/*!
  Returns genetic unit number <num_unit> (0 for main chromosome)
*/
inline ae_genetic_unit* ae_individual::get_genetic_unit( int16_t num_unit ) const
{
  ae_list_node* gen_unit_node = _genetic_unit_list->get_first();
  
  for ( int16_t i = 0 ; i < num_unit ; i++ )
  {
    gen_unit_node = gen_unit_node->get_next();
  }
  
  return (ae_genetic_unit*) gen_unit_node->get_obj();
}

inline double ae_individual::get_dist_to_target_by_feature( ae_env_axis_feature feature ) const
{
  assert( _distance_to_target_computed );
  
  return _dist_to_target_by_feature[feature];
}

inline double ae_individual::get_fitness( void ) const
{
  assert( _fitness_computed );
  
  return _fitness;
}

inline double ae_individual::get_fitness_by_feature( ae_env_axis_feature feature ) const
{
  assert( _fitness_computed );
  
  return _fitness_by_feature[feature];
}

inline ae_grid_cell* ae_individual::get_grid_cell( void ) const
{
  return _grid_cell;
}

inline bool ae_individual::get_placed_in_population ( void ) const
{
  return _placed_in_population;
}

/*!
  Returns the sequence of genetic unit number <num_unit> (0 for main chromosome)
*/
inline const char* ae_individual::get_genetic_unit_sequence( int16_t num_unit ) const
{
  return get_genetic_unit( num_unit )->get_sequence();
}

/*!
  Returns the sequence length of genetic unit number <num_unit> (0 for main chromosome)
*/
inline int32_t ae_individual::get_genetic_unit_seq_length( int16_t num_unit ) const
{
  return get_genetic_unit( num_unit )->get_seq_length();
}

ae_fuzzy_set* ae_individual::get_phenotype_activ( void ) const
{
  return _phenotype_activ;
}

ae_fuzzy_set* ae_individual::get_phenotype_inhib( void ) const
{
  return _phenotype_inhib;
}

ae_phenotype* ae_individual::get_phenotype( void ) const
{
  return _phenotype;
}

ae_replication_report* ae_individual::get_replic_report( void ) const
{
  return _replic_report;
}

inline ae_list* ae_individual::get_protein_list( void ) const
{
  assert( _protein_list );
  
  return _protein_list;
}

inline ae_list* ae_individual::get_rna_list( void ) const
{
  return _rna_list;
}

inline int32_t ae_individual::get_nb_plasmids( void ) const
{
  return _genetic_unit_list->get_nb_elts() - 1;
}

inline int32_t ae_individual::get_nb_gen_units( void ) const
{
  return _genetic_unit_list->get_nb_elts();
}

inline int32_t ae_individual::get_total_genome_size( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _total_genome_size;
}


inline int32_t ae_individual::get_nb_coding_RNAs( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _nb_coding_RNAs;
}

inline int32_t ae_individual::get_nb_non_coding_RNAs( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _nb_non_coding_RNAs;
}

inline double ae_individual::get_overall_size_coding_RNAs( void ) 
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _overall_size_coding_RNAs;
}

inline double ae_individual::get_av_size_coding_RNAs( void ) 
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( _nb_coding_RNAs != 0 )
  {
    return _overall_size_coding_RNAs / _nb_coding_RNAs;
  }
  else return 0.0;
}

inline double ae_individual::get_overall_size_non_coding_RNAs( void ) 
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _overall_size_non_coding_RNAs;
}

inline double ae_individual::get_av_size_non_coding_RNAs( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( _nb_non_coding_RNAs != 0 )
  {
    return _overall_size_non_coding_RNAs / _nb_non_coding_RNAs;
  }
  else return 0.0;
}

inline int32_t ae_individual::get_nb_genes_activ( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _nb_genes_activ;
}

inline int32_t ae_individual::get_nb_genes_inhib( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _nb_genes_inhib;
}

inline int32_t ae_individual::get_nb_functional_genes( void ) 
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _nb_functional_genes;
}

inline int32_t ae_individual::get_nb_non_functional_genes( void ) 
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _nb_non_functional_genes;
}

inline double ae_individual::get_overall_size_functional_genes( void ) 
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _overall_size_functional_genes;
}

inline double ae_individual::get_av_size_functional_genes( void ) 
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( _nb_functional_genes != 0 )
  {
    return _overall_size_functional_genes / _nb_functional_genes;
  }
  else return 0.0;
}

inline double ae_individual::get_overall_size_non_functional_genes( void ) 
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  return _overall_size_non_functional_genes;
}

inline double ae_individual::get_av_size_non_functional_genes( void ) 
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( _nb_non_functional_genes != 0 )
  {
    return _overall_size_non_functional_genes / _nb_non_functional_genes;
  }
  else return 0.0;
}

inline int32_t ae_individual::get_nb_bases_in_0_CDS( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_CDS;
}

inline int32_t ae_individual::get_nb_bases_in_0_functional_CDS( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_functional_CDS;
}

inline int32_t ae_individual::get_nb_bases_in_0_non_functional_CDS( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_non_functional_CDS;
}

inline int32_t ae_individual::get_nb_bases_in_0_RNA( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_RNA;
}

inline int32_t ae_individual::get_nb_bases_in_0_coding_RNA( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_coding_RNA;
}

inline int32_t ae_individual::get_nb_bases_in_0_non_coding_RNA( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_non_coding_RNA;
}

inline int32_t ae_individual::get_nb_bases_in_neutral_regions( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_neutral_regions;
}

inline int32_t ae_individual::get_nb_neutral_regions( void )
{
  if ( ! _statistical_data_computed ) compute_statistical_data();
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_neutral_regions;
}

inline double ae_individual::get_modularity( void )
{
  #warning modularity measure not yet implemented
  //~ if ( _modularity < 0 ) compute_modularity();
  //~ return _modularity;
  return 0;
}

inline double ae_individual::get_w_max( void ) const
{
  return _w_max;
}
    
// ------------------------------------------------------------- Mutation rates
inline double ae_individual::get_point_mutation_rate( void ) const
{
  return _mut_params->get_point_mutation_rate();
}

inline double ae_individual::get_small_insertion_rate( void ) const
{
  return _mut_params->get_small_insertion_rate();
}

inline double ae_individual::get_small_deletion_rate( void ) const
{
  return _mut_params->get_small_deletion_rate();
}

inline int16_t ae_individual::get_max_indel_size( void ) const
{
  return _mut_params->get_max_indel_size();
}

// ---------------------------------- Rearrangement rates (without alignements)
inline double ae_individual::get_duplication_rate( void ) const
{
  return _mut_params->get_duplication_rate();
}

inline double ae_individual::get_deletion_rate( void ) const
{
  return _mut_params->get_deletion_rate();
}

inline double ae_individual::get_translocation_rate( void ) const
{
  return _mut_params->get_translocation_rate();
}

inline double ae_individual::get_inversion_rate( void ) const
{
  return _mut_params->get_inversion_rate();
}

// ------------------------------------- Rearrangement rates (with alignements)
inline double ae_individual::get_neighbourhood_rate( void ) const
{
  return _mut_params->get_neighbourhood_rate();
}

inline double ae_individual::get_duplication_proportion( void ) const
{
  return _mut_params->get_duplication_proportion();
}

inline double ae_individual::get_deletion_proportion( void ) const
{
  return _mut_params->get_deletion_proportion();
}

inline double ae_individual::get_translocation_proportion( void ) const
{
  return _mut_params->get_translocation_proportion();
}

inline double ae_individual::get_inversion_proportion( void ) const
{
  return _mut_params->get_inversion_proportion();
}


inline void ae_individual::set_w_max( double w_max )
{
  _w_max = w_max;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_individual::renew_dist_to_target_by_feature( void )
{
  if ( _dist_to_target_by_feature != NULL ) delete [] _dist_to_target_by_feature;
  _dist_to_target_by_feature = new double [NB_FEATURES];
}


inline void ae_individual::renew_fitness_by_feature( void )
{
  if ( _fitness_by_feature != NULL ) delete [] _fitness_by_feature;
  _fitness_by_feature = new double [NB_FEATURES];
}

inline void ae_individual::do_transcription( void )
{
  if ( _transcribed == true ) return; // Transcription has already been performed, nothing to do.
  _transcribed = true;
  
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*)gen_unit_node->get_obj();
    
    gen_unit->do_transcription();
     
    gen_unit_node = gen_unit_node->get_next();
  }
}

inline void ae_individual::do_translation( void )
{
  if ( _translated == true ) return; // ARNs have already been translated, nothing to do.
  _translated = true;
  if ( _transcribed == false ) do_transcription();
  
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*)gen_unit_node->get_obj();
    
    gen_unit->do_translation();
     
    gen_unit_node = gen_unit_node->get_next();
  }
}

inline void ae_individual::do_folding( void )
{
  if ( _folded == true ) return; // Proteins have already been folded, nothing to do.
  _folded = true;
  if ( _translated == false ) do_translation();
  
  ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
  ae_genetic_unit*  gen_unit;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*)gen_unit_node->get_obj();
    
    gen_unit->compute_phenotypic_contribution();
     
    gen_unit_node = gen_unit_node->get_next();
  }
}

void ae_individual::do_transcription_translation_folding( void )
{
  if ( _transcribed == true && _translated == true && _folded == true ) return;
  
  do_transcription();
  do_translation();
  do_folding();
  
  make_protein_list();
}

#ifdef DEBUG
  inline void ae_individual::print_protein_list( void )
  {
    ae_list_node* prot_node = _protein_list->get_first();
    ae_protein*   prot = NULL;
    
    while ( prot_node != NULL )
    {
      prot = (ae_protein*) prot_node->get_obj();
      
      char* prot_sequence = prot->get_AA_sequence();
      printf( "prot at pos : %"PRId32"      length : %"PRId32" AAs\n", prot->get_first_translated_pos(), prot->get_length() );
      printf( "  strand : %s    concentration : %f\n", (prot->get_strand() == LEADING)?"LEADING":"LAGGING", prot->get_concentration() );
      printf( "  sequence : %s\n", prot_sequence );
      delete prot_sequence;
      
      prot_node = prot_node->get_next();
    }
  }
  
  inline void ae_individual::assert_promoters( void )
  {
    // Perform assertion for each genetic unit
    ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
    ae_genetic_unit*  gen_unit      = NULL;
    
    while ( gen_unit_node != NULL )
    {
      gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
      
      gen_unit->assert_promoters();

      gen_unit_node = gen_unit_node->get_next();
    }
  }
  
  inline void ae_individual::assert_promoters_order( void )
  {
    // Perform assertion for each genetic unit
    ae_list_node*     gen_unit_node = _genetic_unit_list->get_first();
    ae_genetic_unit*  gen_unit      = NULL;
    
    while ( gen_unit_node != NULL )
    {
      gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
      
      gen_unit->assert_promoters_order();

      gen_unit_node = gen_unit_node->get_next();
    }
  }
#endif


#endif // __AE_INDIVIDUAL_H__
