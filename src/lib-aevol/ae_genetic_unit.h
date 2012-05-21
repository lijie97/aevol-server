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


#ifndef __AE_GENETIC_UNIT_H__
#define __AE_GENETIC_UNIT_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_dna.h>
#include <ae_rna.h>
#include <ae_protein.h>
#include <ae_fuzzy_set.h>
#include <ae_common.h>
#include <ae_environment.h>
#include <ae_utils.h>



// =================================================================
//                          Class declarations
// =================================================================






class ae_genetic_unit : public ae_object
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_genetic_unit( ae_individual* indiv, int32_t length );
    ae_genetic_unit( ae_individual* indiv, char* seq, int32_t length, ae_list** prom_list = NULL );
    ae_genetic_unit( ae_individual* indiv, const ae_genetic_unit &model );
    ae_genetic_unit( ae_individual* indiv, ae_genetic_unit* const parent );
    ae_genetic_unit( ae_individual* indiv, gzFile* backup_file );
    ae_genetic_unit( ae_individual* indiv, char* organism_file_name );


    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_genetic_unit( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_individual* get_indiv( void )                   const;
    inline ae_dna*        get_dna( void )                     const;
    inline ae_list**      get_rna_list( void )                const;
    inline void           set_rna_list( ae_list** new_list )       ;
    inline ae_list**      get_protein_list( void )            const;
    inline ae_fuzzy_set*  get_activ_contribution( void )      const;
    inline ae_fuzzy_set*  get_inhib_contribution( void )      const;
    inline ae_fuzzy_set*  get_phenotypic_contribution( void ) const;
    
    
    // Direct DNA access
    inline const char*  get_sequence( void ) const;
    inline int32_t      get_seq_length( void ) const;

    
    // Statistical data
    inline int32_t  get_nb_coding_RNAs( void )                    const;
    inline int32_t  get_nb_non_coding_RNAs( void )                const;
    inline double   get_overall_size_coding_RNAs( void )          const;
    inline double   get_av_size_coding_RNAs( void )               const;
    inline double   get_overall_size_non_coding_RNAs( void )      const;
    inline double   get_av_size_non_coding_RNAs( void )           const;
    inline int32_t  get_nb_genes_activ( void )                    const;
    inline int32_t  get_nb_genes_inhib( void )                    const;
    inline int32_t  get_nb_functional_genes( void )               const;
    inline int32_t  get_nb_non_functional_genes( void )           const;
    inline double   get_overall_size_functional_genes( void )     const;
    inline double   get_av_size_functional_genes( void )          const;
    inline double   get_overall_size_non_functional_genes( void ) const;
    inline double   get_av_size_non_functional_genes( void )      const;
    
    inline int32_t  get_nb_bases_in_0_CDS( void );
    inline int32_t  get_nb_bases_in_0_functional_CDS( void );
    inline int32_t  get_nb_bases_in_0_non_functional_CDS( void );
    inline int32_t  get_nb_bases_in_0_RNA( void );
    inline int32_t  get_nb_bases_in_0_coding_RNA( void );
    inline int32_t  get_nb_bases_in_0_non_coding_RNA( void );
    
    inline int32_t  get_nb_bases_non_essential( void );
    inline int32_t  get_nb_bases_non_essential_including_nf_genes( void );
    
    inline int32_t  get_nb_bases_in_neutral_regions( void );
    inline int32_t  get_nb_neutral_regions( void );
    inline int32_t* get_beginning_neutral_regions( void );
    inline int32_t* get_end_neutral_regions( void );
    
    inline double   get_modularity( void )                        const;

    inline double   get_dist_to_target_by_feature( ae_env_axis_feature feature )  const;
    inline double   get_fitness( void )                                           const;
    inline double   get_fitness_by_feature( ae_env_axis_feature feature )         const;
    
    
    
    // =================================================================
    //                            Public Methods
    // =================================================================
    void locate_promoters( void );
    void do_transcription( void );
    void do_translation( void );
    void compute_phenotypic_contribution( void );
    
    // DM: these two are identical to functions from ae_individual 
    void compute_distance_to_target( ae_environment* envir );
    void compute_fitness( ae_environment* envir );
    
    void reset_expression( void ); // useful for post-treatment programs
    
    inline void print_rnas( void ) const;
    inline static void print_rnas( ae_list ** rnas );
    static void print_rnas( ae_list * rnas, ae_strand strand );
    void print_proteins( void ) const;

    bool        is_promoter( ae_strand strand, int32_t pos, int8_t& dist ) const;
    bool        is_terminator( ae_strand strand, int32_t pos ) const;
    bool        is_shine_dalgarno( ae_strand strand, int32_t pos ) const;
    inline bool is_start( ae_strand strand, int32_t pos ) const;
    inline bool is_stop( ae_strand strand, int32_t pos ) const;
    int8_t      get_codon( ae_strand strand, int32_t pos ) const;
    
    void compute_non_coding( void );
  
    
    void duplicate_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list** duplicated_promoters );

    void get_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list** promoters );
    void get_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list* leading_promoters );
    void get_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list* lagging_promoters );
    void get_leading_promoters_starting_after( int32_t pos, ae_list* leading_promoters );
    void get_leading_promoters_starting_before( int32_t pos, ae_list* leading_promoters );
    void get_lagging_promoters_starting_before( int32_t pos, ae_list* lagging_promoters );
    void get_lagging_promoters_starting_after( int32_t pos, ae_list* lagging_promoters );
    
    void invert_promoters_included_in( int32_t pos_1, int32_t pos_2 );
    static void invert_promoters( ae_list** promoter_lists, int32_t seq_length );
    static void invert_promoters( ae_list** promoter_lists, int32_t pos_1, int32_t pos_2 ); // See WARNING
    
    inline void extract_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list** extracted_promoters );
    inline void extract_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list** extracted_promoters );
    void extract_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list* extracted_promoters );
    void extract_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list* extracted_promoters );
    void extract_leading_promoters_starting_after( int32_t pos, ae_list* extracted_promoters );
    void extract_leading_promoters_starting_before( int32_t pos, ae_list* extracted_promoters );
    void extract_lagging_promoters_starting_before( int32_t pos, ae_list* extracted_promoters );
    void extract_lagging_promoters_starting_after( int32_t pos, ae_list* extracted_promoters );
    
    static void shift_promoters( ae_list** promoters_to_shift, int32_t delta_pos, int32_t seq_length );
    void insert_promoters( ae_list** promoters_to_insert );
    void insert_promoters_at( ae_list** promoters_to_insert, int32_t pos );
    
    inline void remove_promoters_around( int32_t pos );
    inline void remove_promoters_around( int32_t pos_1, int32_t pos_2 );
    inline void remove_all_promoters( void );
    
    inline void look_for_new_promoters_around( int32_t pos );
    inline void look_for_new_promoters_around( int32_t pos_1, int32_t pos_2 );
    
    inline void move_all_promoters_after( int32_t pos, int32_t delta_pos );
    void move_all_leading_promoters_after( int32_t pos, int32_t delta_pos );
    void move_all_lagging_promoters_after( int32_t pos, int32_t delta_pos );
    
    // Not implemented!!!
    //~ inline void duplicate_promoters_starting_between( int32_t pos_1, int32_t pos_2, int32_t delta_pos );
    //~ void duplicate_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, int32_t delta_pos );
    //~ void duplicate_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, int32_t delta_pos );
    
    inline void copy_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list** new_promoter_lists );
    inline void copy_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list** new_promoter_lists );
    void copy_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list* new_promoter_list );
    void copy_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list* new_promoter_list );
    //~ inline void copy_all_promoters( ae_list** new_promoter_lists );

    void write_to_backup( gzFile* backup_file );
    
    int32_t get_nb_terminators( void );
    
    //~ // set the genetic unit of all promoters to the appropriate 
    //~ void reasign_promoter_genetic_unit (void);
    
    #ifdef DEBUG
      void assert_promoters( void );
      void assert_promoters_order( void );
    #endif

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_genetic_unit( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_genetic_unit( const ae_genetic_unit &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };


    // =================================================================
    //                           Protected Methods
    // =================================================================
    void init_statistical_data( void );
    
    void remove_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2 );
    void remove_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2 );
    void remove_leading_promoters_starting_before( int32_t pos );
    void remove_lagging_promoters_starting_before( int32_t pos );
    void remove_leading_promoters_starting_after( int32_t pos );
    void remove_lagging_promoters_starting_after( int32_t pos );
    
    void look_for_new_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2 );
    void look_for_new_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2 );
    void look_for_new_leading_promoters_starting_after( int32_t pos );
    void look_for_new_lagging_promoters_starting_after( int32_t pos );
    void look_for_new_leading_promoters_starting_before( int32_t pos );
    void look_for_new_lagging_promoters_starting_before( int32_t pos );


    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_individual*  _indiv;
    ae_dna*         _dna;
    ae_list**       _rna_list;
    ae_list**       _protein_list;
    ae_fuzzy_set*   _activ_contribution;
    ae_fuzzy_set*   _inhib_contribution;
    ae_fuzzy_set*   _phenotypic_contribution;
    // NB : _phenotypic_contribution is only an indicative value, not used for the whole phenotype computation

    // DM: For plasmid work, we sometimes *need* all the data (e.g. fitness, secretion) calculated for each GU
    double* _dist_to_target_per_segment;
    double* _dist_to_target_by_feature;
    double* _fitness_by_feature;
    double  _fitness;    
    
    // Statistical data (intrinsic for this genetic unit)
    int32_t _nb_coding_RNAs;                // Number of coding RNAs (at least one gene on RNA)
    int32_t _nb_non_coding_RNAs;            // Number of non-coding-RNAs
    double  _overall_size_coding_RNAs;      // Average size of coding RNAs
    double  _overall_size_non_coding_RNAs;  // Average size of non-coding RNAs
    int32_t _nb_genes_activ;                // Number of genes realizing a function
    int32_t _nb_genes_inhib;                // Number of genes inhibitting a function
    int32_t _nb_fun_genes;                  // Number of functional genes
    int32_t _nb_non_fun_genes;              // Number of non-functional genes
    double  _overall_size_fun_genes;        // Average size of functional genes
    double  _overall_size_non_fun_genes;    // Average size of non-functional genes
    
    
    int32_t _nb_bases_in_0_CDS;                // Number of bases that are not included in any gene
    int32_t _nb_bases_in_0_functional_CDS;     // Number of bases that are not included in any non-degenerated gene
    int32_t _nb_bases_in_0_non_functional_CDS; // Number of bases that are not included in any degenerated gene
    
    int32_t _nb_bases_in_0_RNA;                // Number of bases that are not included in any RNA
    int32_t _nb_bases_in_0_coding_RNA;         // Number of bases that are not included in any coding RNA
                                               // (RNAs containing at least one CDS)
    int32_t _nb_bases_in_0_non_coding_RNA;     // Number of bases that are not included in any non coding RNA
    
    int32_t _nb_bases_non_essential;                    // Number of bases that are not included in any non-degenerated gene or involved in its expression
    int32_t _nb_bases_non_essential_including_nf_genes; // Number of bases that are not included in any gene or involved in its expression
    int32_t _nb_bases_in_neutral_regions;   // Number of bases that are not included in a neutral region
                                            // A base is considered neutral when neither itself NOR its corresponding base on the other
                                            // strand belongs to a coding promoter->terminator region (both included)
    int32_t _nb_neutral_regions;            // Number of neutral regions
    int32_t* _beginning_neutral_regions;    // Beginning of neutral regions
    int32_t* _end_neutral_regions;          // Corresponding ends of neutral regions

    double  _modularity;  // Ratio between the pairwise distance between genes whose corresponding
                          // phenotypic triangles overlap and the average intergenic distance 
                          // (ignoring non-functional genes)
    
    
    // TODO : check and comment what it is for
    //~ int32_t _nb_promoters[2];
    //~ int32_t _nb_genes[2];
    //~ int32_t _nb_functional_genes[2]; 
    //~ double  _average_gene_size;
    //~ double  _average_functional_gene_size;
    //~ int32_t _nb_coding_bp;
    //~ double  _clustering;
    
    // We keep trace of what we have already computed to avoid double computation (mainly in post-treaments)
    bool _transcribed;
    bool _translated;
    bool _phenotypic_contributions_computed;
    bool _non_coding_computed;
    bool _distance_to_target_computed;
    bool _fitness_computed; 
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline ae_individual* ae_genetic_unit::get_indiv( void ) const
{
  return _indiv;
}

inline ae_dna* ae_genetic_unit::get_dna( void ) const
{
  return _dna;
}

inline ae_list** ae_genetic_unit::get_rna_list( void ) const
{
  return _rna_list;
}

inline void ae_genetic_unit::set_rna_list( ae_list** new_list )
{
  _rna_list = new_list;
} // TODO : erase that! // NOTE : Why?

inline ae_list** ae_genetic_unit::get_protein_list( void ) const
{
  assert( _protein_list );
  assert( _protein_list[LEADING] );
  assert( _protein_list[LAGGING] );

  return _protein_list;
}

inline ae_fuzzy_set* ae_genetic_unit::get_activ_contribution( void ) const
{
  return _activ_contribution;
}

inline ae_fuzzy_set* ae_genetic_unit::get_inhib_contribution( void ) const
{
  return _inhib_contribution;
}

inline ae_fuzzy_set* ae_genetic_unit::get_phenotypic_contribution( void ) const
{
  return _phenotypic_contribution;
}

/*!
  Returns the DNA sequence
*/
inline const char* ae_genetic_unit::get_sequence( void ) const
{
  return _dna->get_data();
}

/*!
  Returns the DNA sequence length
*/
inline int32_t ae_genetic_unit::get_seq_length( void ) const
{
  return _dna->get_length();
}

inline int32_t ae_genetic_unit::get_nb_coding_RNAs( void ) const
{
  return _nb_coding_RNAs;
}

inline int32_t ae_genetic_unit::get_nb_non_coding_RNAs( void ) const
{
  return _nb_non_coding_RNAs;
}

inline double ae_genetic_unit::get_overall_size_coding_RNAs( void ) const
{
  return _overall_size_coding_RNAs;
}

inline double ae_genetic_unit::get_av_size_coding_RNAs( void ) const
{
  if ( _nb_coding_RNAs != 0 )
  {
    return _overall_size_coding_RNAs / _nb_coding_RNAs;
  }
  else return 0.0;
}

inline double ae_genetic_unit::get_overall_size_non_coding_RNAs( void ) const
{
  return _overall_size_non_coding_RNAs;
}

inline double ae_genetic_unit::get_av_size_non_coding_RNAs( void ) const
{
  if ( _nb_non_coding_RNAs != 0 )
  {
    return _overall_size_non_coding_RNAs / _nb_non_coding_RNAs;
  }
  else return 0.0;
}

inline int32_t ae_genetic_unit::get_nb_genes_activ( void ) const
{
  return _nb_genes_activ;
}

inline int32_t ae_genetic_unit::get_nb_genes_inhib( void ) const
{
  return _nb_genes_inhib;
}

inline int32_t ae_genetic_unit::get_nb_functional_genes( void ) const
{
  return _nb_fun_genes;
}

inline int32_t ae_genetic_unit::get_nb_non_functional_genes( void ) const
{
  return _nb_non_fun_genes;
}

inline double ae_genetic_unit::get_overall_size_functional_genes( void ) const
{
  return _overall_size_fun_genes;
}

inline double ae_genetic_unit::get_av_size_functional_genes( void ) const
{
  if ( _nb_fun_genes != 0 )
  {
    return _overall_size_fun_genes / _nb_fun_genes;
  }
  else return 0.0;
}

inline double ae_genetic_unit::get_overall_size_non_functional_genes( void ) const
{
  return _overall_size_non_fun_genes;
}

inline double ae_genetic_unit::get_av_size_non_functional_genes( void ) const
{
  if ( _nb_non_fun_genes != 0 )
  {
    return _overall_size_non_fun_genes / _nb_non_fun_genes;
  }
  else return 0.0;
}

inline int32_t ae_genetic_unit::get_nb_bases_in_0_CDS( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_CDS;
}

inline int32_t ae_genetic_unit::get_nb_bases_in_0_functional_CDS( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_functional_CDS;
}

inline int32_t ae_genetic_unit::get_nb_bases_in_0_non_functional_CDS( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_non_functional_CDS;
}

inline int32_t ae_genetic_unit::get_nb_bases_in_0_RNA( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_RNA;
}

inline int32_t ae_genetic_unit::get_nb_bases_in_0_coding_RNA( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_coding_RNA;
}

inline int32_t ae_genetic_unit::get_nb_bases_in_0_non_coding_RNA( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_0_non_coding_RNA;
}

inline int32_t ae_genetic_unit::get_nb_bases_non_essential( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_non_essential;
}

inline int32_t ae_genetic_unit::get_nb_bases_non_essential_including_nf_genes( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_non_essential_including_nf_genes;
}

inline int32_t ae_genetic_unit::get_nb_bases_in_neutral_regions( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_bases_in_neutral_regions;
}

inline int32_t ae_genetic_unit::get_nb_neutral_regions( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _nb_neutral_regions;
}

inline int32_t* ae_genetic_unit::get_beginning_neutral_regions( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _beginning_neutral_regions;
}

inline int32_t* ae_genetic_unit::get_end_neutral_regions( void )
{
  if ( ! _non_coding_computed ) compute_non_coding();
  return _end_neutral_regions;
}

inline double ae_genetic_unit::get_modularity( void ) const
{
  return _modularity;
}

inline double ae_genetic_unit::get_dist_to_target_by_feature( ae_env_axis_feature feature ) const
{
   assert( _distance_to_target_computed );
  
  return _dist_to_target_by_feature[feature];
}

inline double ae_genetic_unit::get_fitness( void ) const
{
  assert( _fitness_computed );
  
  return _fitness;
}

inline double ae_genetic_unit::get_fitness_by_feature( ae_env_axis_feature feature ) const
{
  assert( _fitness_computed );
  
  return _fitness_by_feature[feature];
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_genetic_unit::print_rnas( void ) const
{
  print_rnas( _rna_list );
}

inline /* static */ void ae_genetic_unit::print_rnas( ae_list ** rnas )
{
  print_rnas( rnas[LEADING], LEADING );
  print_rnas( rnas[LAGGING], LAGGING );
}

inline bool ae_genetic_unit::is_start( ae_strand strand, int32_t index ) const
{
  return ( get_codon( strand, index ) == CODON_START );
}

inline bool ae_genetic_unit::is_stop( ae_strand strand, int32_t index ) const
{
  return ( get_codon( strand, index ) == CODON_STOP );
}

inline void ae_genetic_unit::remove_all_promoters( void )
{
  _rna_list[LEADING]->erase( DELETE_OBJ );
  _rna_list[LAGGING]->erase( DELETE_OBJ );
}

inline void ae_genetic_unit::move_all_promoters_after( int32_t pos, int32_t delta_pos )
{
  move_all_leading_promoters_after( pos, delta_pos );
  move_all_lagging_promoters_after( pos, delta_pos );
}

inline void ae_genetic_unit::extract_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list** extracted_promoters )
{
  assert( pos_1 >= 0 && pos_1 < pos_2 && pos_2 <= _dna->get_length() );
  if ( pos_2 - pos_1 >= PROM_SIZE )
  {
    extract_leading_promoters_starting_between( pos_1, pos_2 - PROM_SIZE + 1, extracted_promoters[LEADING] );
    extract_lagging_promoters_starting_between( pos_1 + PROM_SIZE - 1, pos_2, extracted_promoters[LAGGING] );
  }
}

inline void ae_genetic_unit::extract_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list** extracted_promoters )
{
  extract_leading_promoters_starting_between( pos_1, pos_2, extracted_promoters[LEADING] );
  extract_lagging_promoters_starting_between( pos_1, pos_2, extracted_promoters[LAGGING] );
}


/*!
  \brief  Remove those promoters that would be broken if the chromosome was cut at pos.

  Remove promoters that include BOTH the base before AND after pos (marked X in the cartoon below).
  If the genome is smaller than the size of a promoter, all the promoters will be removed.
  
  \verbatim
     -------------------------------------------------------
    |   |   |   |   | X | X |   |   |   |   |   |   |   |   |
     -------------------------------------------------------
    ^                   ^
    0                  pos
  \endverbatim
*/
inline void ae_genetic_unit::remove_promoters_around( int32_t pos )
{
  if ( _dna->get_length() >= PROM_SIZE )
  {
    remove_leading_promoters_starting_between( ae_utils::mod(pos - PROM_SIZE + 1, _dna->get_length()), pos );
    remove_lagging_promoters_starting_between( pos, ae_utils::mod(pos + PROM_SIZE - 1, _dna->get_length()) );
  }
  else
  {
    remove_all_promoters();
  }
}


/*!
  \brief  Remove those promoters that would be broken if the sequence [pos_1 ; pos_2[ was deleted.

  Remove promoters that     * include BOTH the base before AND after pos_1 (marked X in the cartoon below).
                            * include BOTH the base before AND after pos_2 (marked Y in the cartoon below).
                            * are completely contained between pos_1 and pos_2.
  If the remaining sequence, i.e. [pos_2 ; pos_1[ is smaller than the size of a promoter, all the promoters will be removed.

  \verbatim
     -------------------------------------------------------
    |   |   |   |   | X | X |   |   |   | Y | Y |   |   |   |
     -------------------------------------------------------
    ^                   ^                   ^
    0                 pos_1               pos_2
  \endverbatim
*/
inline void ae_genetic_unit::remove_promoters_around( int32_t pos_1, int32_t pos_2 )
{
  if ( ae_utils::mod(pos_1 - pos_2, _dna->get_length()) >= PROM_SIZE )
  {
    remove_leading_promoters_starting_between( ae_utils::mod(pos_1 - PROM_SIZE + 1, _dna->get_length()), pos_2 );
    remove_lagging_promoters_starting_between( pos_1, ae_utils::mod(pos_2 + PROM_SIZE - 1, _dna->get_length()) );
  }
  else
  {
    remove_all_promoters();
  }
}


/*!
  \brief  Look for promoters that are astride pos and add them to the list of promoters (_rna_list).

  Look for promoters that include BOTH the base before AND after pos (marked X in the cartoon below).
  If the genome is smaller than the size of a promoter, no search is performed.
  
  \verbatim
     -------------------------------------------------------
    |   |   |   |   | X | X |   |   |   |   |   |   |   |   |
     -------------------------------------------------------
    ^                   ^
    0                  pos
  \endverbatim
*/
inline void ae_genetic_unit::look_for_new_promoters_around( int32_t pos )
{
  assert( pos >= 0 && pos <= _dna->get_length() );
  
  if ( _dna->get_length() >= PROM_SIZE )
  {
    look_for_new_leading_promoters_starting_between(  ae_utils::mod(pos - PROM_SIZE + 1, _dna->get_length()),
                                                      ae_utils::mod(pos                , _dna->get_length()) );
    look_for_new_lagging_promoters_starting_between(  ae_utils::mod(pos                , _dna->get_length()),
                                                      ae_utils::mod(pos + PROM_SIZE - 1, _dna->get_length()) );
  }
}


/*!
  \brief  Look for promoters that contain at least 1 base lying in [pos_1 ; pos_2[ and add them to the list of promoters (_rna_list).

  Look for promoters that   * include BOTH the base before AND after pos_1 (marked X in the cartoon below).
                            * include BOTH the base before AND after pos_2 (marked Y in the cartoon below).
                            * are completely contained between pos_1 and pos_2.
  If the genome is smaller than the size of a promoter, no search is performed.
  
  \verbatim
     -------------------------------------------------------
    |   |   |   |   | X | X |   |   |   | Y | Y |   |   |   |
     -------------------------------------------------------
    ^                   ^                   ^
    0                 pos_1               pos_2
  \endverbatim
*/
inline void ae_genetic_unit::look_for_new_promoters_around( int32_t pos_1, int32_t pos_2 )
{
  //~ if ( ae_utils::mod( pos_1 - pos_2, _dna->get_length()) == PROM_SIZE - 1 )
  //~ {
    //~ // We have to look at every possible position on the genome.
    //~ locate_promoters();
  //~ }
  /*else*/ if ( _dna->get_length() >= PROM_SIZE )
  {
    look_for_new_leading_promoters_starting_between( ae_utils::mod(pos_1 - PROM_SIZE + 1, _dna->get_length()), pos_2 );
    look_for_new_lagging_promoters_starting_between( pos_1, ae_utils::mod(pos_2 + PROM_SIZE - 1, _dna->get_length()) );
  }
}

//~ inline void ae_genetic_unit::duplicate_promoters_starting_between( int32_t pos_1, int32_t pos_2, int32_t delta_pos )
//~ {
  //~ duplicate_leading_promoters_starting_between( pos_1, pos_2, delta_pos );
  //~ duplicate_lagging_promoters_starting_between( pos_1, pos_2, delta_pos );
//~ }

inline void ae_genetic_unit::copy_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list** new_promoter_lists )
{
  if ( ae_utils::mod( pos_2 - pos_1 - 1, _dna->get_length() ) + 1 >= PROM_SIZE )
  {
    copy_leading_promoters_starting_between( pos_1, ae_utils::mod( pos_2 - PROM_SIZE + 1, _dna->get_length() ), new_promoter_lists[LEADING] );
    copy_lagging_promoters_starting_between( ae_utils::mod( pos_1 + PROM_SIZE - 1, _dna->get_length() ), pos_2, new_promoter_lists[LAGGING] );
  }
}

inline void ae_genetic_unit::copy_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list** new_promoter_lists )
{
  copy_leading_promoters_starting_between( pos_1, pos_2, new_promoter_lists[LEADING] );
  copy_lagging_promoters_starting_between( pos_1, pos_2, new_promoter_lists[LAGGING] );
}

//~ inline void ae_genetic_unit::copy_all_promoters( ae_list** new_promoter_lists )
//~ {
  //~ ae_list_node* rna_node = NULL;
  
  //~ for ( int8_t strand = LEADING ; strand <= LAGGING ; strand++ )
  //~ {
    //~ rna_node = _rna_list[strand]->get_first();
    
    //~ while ( rna_node != NULL )
    //~ {
      //~ new_promoter_lists[strand]->add( new ae_rna( this, *((ae_rna*)rna_node->get_obj()) ) );
      
      //~ rna_node = rna_node->get_next();
    //~ }
  //~ }
//~ }

#endif // __ae_genetic_unit_H__
