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
#include "ae_list.h"
#include "ae_dna.h"
#include "ae_rna.h"
#include "ae_protein.h"
#include "fuzzy.h"
#include "environment.h"
#include "ae_jumping_mt.h"
#include "ae_utils.h"



using std::vector;
using std::list;

namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;






class ae_genetic_unit
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_genetic_unit( ae_individual* indiv, int32_t length, ae_jumping_mt * prng );
    ae_genetic_unit( ae_individual* indiv, char* seq, int32_t length, ae_list<ae_rna*>** prom_list = NULL );
    ae_genetic_unit( ae_individual* indiv, const ae_genetic_unit& model );
    ae_genetic_unit( ae_individual* indiv, const ae_genetic_unit* parent );
    ae_genetic_unit( ae_individual* indiv, gzFile backup_file );
    ae_genetic_unit( ae_individual* indiv, char* organism_file_name );
    ae_genetic_unit() = delete;
    ae_genetic_unit(const ae_genetic_unit &) = delete;

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_genetic_unit( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    ae_exp_manager*  get_exp_m( void ) const;
    ae_individual*   get_indiv( void )                   const;
    ae_dna*          get_dna( void )                     const;
    Fuzzy*    get_activ_contribution( void )      const;
    Fuzzy*    get_inhib_contribution( void )      const;
    Fuzzy*    get_phenotypic_contribution( void ) const;

    // TODO: re constify
    // TODO return as (rvalue?) reference
    /*const*/ std::vector<std::list<ae_rna*>> get_rna_list() const;
    // TODO return as (rvalue?) reference
    const list<ae_protein*> get_protein_list(ae_strand strand) const;



    // Direct DNA access
    const char*  get_sequence( void ) const;
    int32_t      get_seq_length( void ) const;


    // Statistical data
    int32_t  get_nb_coding_RNAs( void )                    const;
    int32_t  get_nb_non_coding_RNAs( void )                const;
    double   get_overall_size_coding_RNAs( void )          const;
    double   get_av_size_coding_RNAs( void )               const;
    double   get_overall_size_non_coding_RNAs( void )      const;
    double   get_av_size_non_coding_RNAs( void )           const;
    int32_t  get_nb_genes_activ( void )                    const;
    int32_t  get_nb_genes_inhib( void )                    const;
    int32_t  get_nb_functional_genes( void )               const;
    int32_t  get_nb_non_functional_genes( void )           const;
    double   get_overall_size_functional_genes( void )     const;
    double   get_av_size_functional_genes( void )          const;
    double   get_overall_size_non_functional_genes( void ) const;
    double   get_av_size_non_functional_genes( void )      const;

    int32_t  get_nb_bases_in_0_CDS( void );
    int32_t  get_nb_bases_in_0_functional_CDS( void );
    int32_t  get_nb_bases_in_0_non_functional_CDS( void );
    int32_t  get_nb_bases_in_0_RNA( void );
    int32_t  get_nb_bases_in_0_coding_RNA( void );
    int32_t  get_nb_bases_in_0_non_coding_RNA( void );

    int32_t  get_nb_bases_non_essential( void );
    int32_t  get_nb_bases_non_essential_including_nf_genes( void );

    int32_t  get_nb_bases_in_neutral_regions( void );
    int32_t  get_nb_neutral_regions( void );
    int32_t* get_beginning_neutral_regions( void );
    int32_t* get_end_neutral_regions( void );

    double   get_modularity( void )                        const;

    double   get_dist_to_target_by_feature( ae_env_axis_feature feature )  const;
    double   get_fitness( void )                                           const;
    double   get_fitness_by_feature( ae_env_axis_feature feature )         const;

    int32_t  get_min_gu_length( void ) const;
    int32_t  get_max_gu_length( void ) const;


    void set_min_gu_length( int32_t min_gu_length );
    void set_max_gu_length( int32_t max_gu_length );

    void set_exp_m( ae_exp_manager* exp_m );


    // =================================================================
    //                            Public Methods
    // =================================================================
    void locate_promoters( void );
    void do_transcription( void );
    void do_translation( void );
    void compute_phenotypic_contribution( void );

    void take_ownership_of_all_rnas(void) { ae_dna::set_GU(get_rna_list(), this); };


  

    // DM: these two are identical to functions from ae_individual
    void compute_distance_to_target( Environment* envir );
    void compute_fitness( Environment* envir );

    void reset_expression( void ); // useful for post-treatment programs

    void print_rnas( void ) const;
    void print_coding_rnas( void );
    static void print_rnas( ae_list<ae_rna*>** rnas );
    static void print_rnas( ae_list<ae_rna*>* rnas, ae_strand strand );
    void print_proteins( void ) const;

    bool        is_promoter( ae_strand strand, int32_t pos, int8_t& dist ) const;
    bool        is_terminator( ae_strand strand, int32_t pos ) const;
    bool        is_shine_dalgarno( ae_strand strand, int32_t pos ) const;
    bool is_start( ae_strand strand, int32_t pos ) const;
    bool is_stop( ae_strand strand, int32_t pos ) const;
    int8_t      get_codon( ae_strand strand, int32_t pos ) const;

    void compute_non_coding( void );


    void duplicate_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** duplicated_promoters );

    void get_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** promoters );
    void get_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* leading_promoters );
    void get_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* lagging_promoters );
    void get_leading_promoters_starting_after( int32_t pos, ae_list<ae_rna*>* leading_promoters );
    void get_leading_promoters_starting_before( int32_t pos, ae_list<ae_rna*>* leading_promoters );
    void get_lagging_promoters_starting_before( int32_t pos, ae_list<ae_rna*>* lagging_promoters );
    void get_lagging_promoters_starting_after( int32_t pos, ae_list<ae_rna*>* lagging_promoters );

    void invert_promoters_included_in( int32_t pos_1, int32_t pos_2 );
    static void invert_promoters( ae_list<ae_rna*>** promoter_lists, int32_t seq_length );
    static void invert_promoters( ae_list<ae_rna*>** promoter_lists, int32_t pos_1, int32_t pos_2 ); // See WARNING

    void extract_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** extracted_promoters );
    void extract_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** extracted_promoters );
    void extract_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* extracted_promoters );
    void extract_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* extracted_promoters );
    void extract_leading_promoters_starting_after( int32_t pos, ae_list<ae_rna*>* extracted_promoters );
    void extract_leading_promoters_starting_before( int32_t pos, ae_list<ae_rna*>* extracted_promoters );
    void extract_lagging_promoters_starting_before( int32_t pos, ae_list<ae_rna*>* extracted_promoters );
    void extract_lagging_promoters_starting_after( int32_t pos, ae_list<ae_rna*>* extracted_promoters );

    static void shift_promoters( ae_list<ae_rna*>** promoters_to_shift, int32_t delta_pos, int32_t seq_length );
    void insert_promoters( ae_list<ae_rna*>** promoters_to_insert );
    void insert_promoters_at( ae_list<ae_rna*>** promoters_to_insert, int32_t pos );

    void remove_promoters_around( int32_t pos );
    void remove_promoters_around( int32_t pos_1, int32_t pos_2 );
    void remove_all_promoters( void );

    void look_for_new_promoters_around( int32_t pos );
    void look_for_new_promoters_around( int32_t pos_1, int32_t pos_2 );

    void move_all_promoters_after( int32_t pos, int32_t delta_pos );
    void move_all_leading_promoters_after( int32_t pos, int32_t delta_pos );
    void move_all_lagging_promoters_after( int32_t pos, int32_t delta_pos );

    // Not implemented!!!
    //~ void duplicate_promoters_starting_between( int32_t pos_1, int32_t pos_2, int32_t delta_pos );
    //~ void duplicate_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, int32_t delta_pos );
    //~ void duplicate_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, int32_t delta_pos );

    void copy_promoters_included_in( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** new_promoter_lists );
    void copy_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>** new_promoter_lists );
    void copy_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* new_promoter_list );
    void copy_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, ae_list<ae_rna*>* new_promoter_list );
    //~ void copy_all_promoters( ae_list** new_promoter_lists );

    void save( gzFile backup_file ) const;

    int32_t get_nb_terminators( void );

    //~ // set the genetic unit of all promoters to the appropriate
    //~ void reasign_promoter_genetic_unit (void);

    #ifdef DEBUG
      void assert_promoters( void );
      void assert_promoters_order( void );
    #endif

    bool* is_belonging_to_coding_RNA(void);
    void remove_non_coding_bases( void);
    void double_non_coding_bases(void);

    // WARNING: The method below works properly only in the case of a single genetic unit (no plasmid).
    // Translocations between different genetic units is not managed.
    void compute_nb_of_affected_genes(const ae_mutation * mut, int & nb_genes_at_breakpoints, int & nb_genes_in_segment, int & nb_genes_in_replaced_segment);

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
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
    ae_exp_manager* _exp_m;

    ae_individual*  _indiv;
    ae_dna*         _dna;
    Fuzzy*   _activ_contribution;
    Fuzzy*   _inhib_contribution;
    Fuzzy*   _phenotypic_contribution;
    // NB : _phenotypic_contribution is only an indicative value, not used for the whole phenotype computation

    ae_list<ae_rna*>**     _rna_list;
    ae_list<ae_protein*>** _protein_list;

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

    int32_t _min_gu_length;
    int32_t _max_gu_length;

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

} // namespace aevol
#endif // __ae_genetic_unit_H__
