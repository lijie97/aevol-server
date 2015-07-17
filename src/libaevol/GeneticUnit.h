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

#ifndef AEVOL_GENETIC_UNIT_H__
#define AEVOL_GENETIC_UNIT_H__



// =================================================================
//                              Includes
// =================================================================
#include <inttypes.h>
#include <stdio.h>

#include <memory>

#include "Dna.h"
#include "Rna.h"
#include "Protein.h"
#include "Fuzzy.h"
#include "JumpingMT.h"
#include "Utils.h"
#include "PhenotypicTarget.h"

#ifdef __REGUL
#include "raevol/Protein_R.h"
#endif

using std::vector;
using std::list;

namespace aevol {

using Promoters1Strand = std::list<Rna>;
using Promoters2Strands = std::vector<Promoters1Strand>;

class ExpManager;

class GeneticUnit
{
 public :

  // =================================================================
  //                             Constructors
  // =================================================================
  GeneticUnit(Individual * indiv,
              int32_t length,
              std::shared_ptr<JumpingMT> prng);
  GeneticUnit(Individual * indiv,
              char* seq,
              int32_t length,
              const Promoters2Strands& prom_list = {{},{}});
  GeneticUnit( Individual * indiv, const GeneticUnit& model );
  GeneticUnit( Individual * indiv, const GeneticUnit* parent );
  GeneticUnit( Individual * indiv, gzFile backup_file );
  GeneticUnit( Individual * indiv, char* organism_file_name );
  GeneticUnit() = delete;
  GeneticUnit(const GeneticUnit &) = delete;

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~GeneticUnit( void );

  // =================================================================
  //                              Accessors
  // =================================================================
  ExpManager* get_exp_m() const;
  Individual* get_indiv() const;
  Dna* get_dna() const;
  Fuzzy* get_activ_contribution() const;
  Fuzzy* get_inhib_contribution() const;
  Fuzzy* get_phenotypic_contribution() const;

  const Promoters2Strands& get_rna_list() const;
  // TODO return as reference
  #ifndef __REGUL
  std::list<Protein>& get_protein_list(Strand strand);
  #else
  std::list<Protein_R>& get_protein_list(Strand strand);
  #endif
  void clear_protein_list(Strand strand);

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

  int32_t  get_nb_bases_in_0_CDS() const;
  int32_t  get_nb_bases_in_0_functional_CDS() const;
  int32_t  get_nb_bases_in_0_non_functional_CDS() const;
  int32_t  get_nb_bases_in_0_RNA() const;
  int32_t  get_nb_bases_in_0_coding_RNA() const;
  int32_t  get_nb_bases_in_0_non_coding_RNA() const;

  int32_t  get_nb_bases_non_essential() const;
  int32_t  get_nb_bases_non_essential_including_nf_genes() const;

  int32_t  get_nb_bases_in_neutral_regions() const;
  int32_t  get_nb_neutral_regions() const;
  int32_t* get_beginning_neutral_regions() const;
  int32_t* get_end_neutral_regions() const;

  double   get_modularity() const;

  double   get_dist_to_target_by_feature(PhenotypicFeature feature) const;
  double   get_fitness() const;
  double   get_fitness_by_feature(PhenotypicFeature feature) const;

  int32_t  get_min_gu_length( void ) const;
  int32_t  get_max_gu_length( void ) const;


  void set_min_gu_length( int32_t min_gu_length );
  void set_max_gu_length( int32_t max_gu_length );

  void set_exp_m( ExpManager * exp_m );


  // =================================================================
  //                            Public Methods
  // =================================================================
  void locate_promoters( void );
  void do_transcription( void );
  void do_translation( void );
  void compute_phenotypic_contribution( void );

  void take_ownership_of_all_rnas(void) { Dna::set_GU(get_rna_list(), this); };




  // DM: these two are identical to functions from Individual
  void compute_distance_to_target(const PhenotypicTarget& target);
  void compute_fitness(const PhenotypicTarget& target);

  void reset_expression( void ); // useful for post-treatment programs

  void print_rnas() const;
  void print_coding_rnas();
  static void print_rnas(const Promoters2Strands& rnas);
  static void print_rnas(const Promoters1Strand& rnas, Strand strand);
  void print_proteins() const;

  bool        is_promoter( Strand strand, int32_t pos, int8_t& dist ) const;
  bool        is_terminator( Strand strand, int32_t pos ) const;
  bool        is_shine_dalgarno( Strand strand, int32_t pos ) const;
  bool is_start( Strand strand, int32_t pos ) const;
  bool is_stop( Strand strand, int32_t pos ) const;
  int8_t      get_codon( Strand strand, int32_t pos ) const;

  void compute_non_coding( void );


  // these functions are called once, they should likely not be public methods
  void duplicate_promoters_included_in(int32_t pos_1,
                                       int32_t pos_2,
                                       Promoters2Strands& duplicated_promoters);

  void get_promoters_included_in(int32_t pos_1,
                                 int32_t pos_2,
                                 Promoters2Strands& promoters);
  void get_promoters(Strand strand_id,
                     Position start,
                     int32_t pos1,
                     int32_t pos2,
                     Promoters1Strand& promoters);

  void invert_promoters_included_in(int32_t pos_1, int32_t pos_2);
  static void invert_promoters(Promoters2Strands& promoter_lists, int32_t seq_length );
  static void invert_promoters(Promoters2Strands& promoter_lists, int32_t pos_1, int32_t pos_2 ); // See WARNING

  void extract_promoters_included_in(int32_t pos_1,
                                     int32_t pos_2,
                                     Promoters2Strands& extracted_promoters);
  void extract_promoters_starting_between(int32_t pos_1,
                                          int32_t pos_2,
                                          Promoters2Strands& extracted_promoters);
  void extract_leading_promoters_starting_between(int32_t pos_1,
                                                  int32_t pos_2,
                                                  Promoters1Strand& extracted_promoters);
  void extract_lagging_promoters_starting_between(int32_t pos_1,
                                                  int32_t pos_2,
                                                  Promoters1Strand& extracted_promoters);





  // end comment

  static void shift_promoters(Promoters2Strands& promoters_to_shift,
                              int32_t delta_pos,
                              int32_t seq_length);
  void insert_promoters(Promoters2Strands& promoters_to_insert);
  void insert_promoters_at(Promoters2Strands& promoters_to_insert,
                           int32_t pos );

  void remove_promoters_around( int32_t pos );
  void remove_promoters_around( int32_t pos_1, int32_t pos_2 );
  void remove_all_promoters( void );

  void look_for_new_promoters_around( int32_t pos );
  void look_for_new_promoters_around( int32_t pos_1, int32_t pos_2 );

  void move_all_promoters_after( int32_t pos, int32_t delta_pos );
  void move_all_leading_promoters_after( int32_t pos, int32_t delta_pos );
  void move_all_lagging_promoters_after( int32_t pos, int32_t delta_pos );

  void copy_promoters_included_in( int32_t pos_1, int32_t pos_2, Promoters2Strands& new_promoter_lists );
  void copy_promoters_starting_between( int32_t pos_1, int32_t pos_2, Promoters2Strands& new_promoter_lists );
  void copy_leading_promoters_starting_between( int32_t pos_1, int32_t pos_2, Promoters1Strand& new_promoter_list );
  void copy_lagging_promoters_starting_between( int32_t pos_1, int32_t pos_2, Promoters1Strand& new_promoter_list );

  void save( gzFile backup_file ) const;

  int32_t get_nb_terminators();

  //~ // set the genetic unit of all promoters to the appropriate
  //~ void reasign_promoter_genetic_unit (void);

  void assert_promoters( void );
  void assert_promoters_order( void );

  bool* is_belonging_to_coding_RNA(void);
  void remove_non_coding_bases( void);
  void double_non_coding_bases(void);

  // WARNING: The method below works properly only in the case of a single genetic unit (no plasmid).
  // Translocations between different genetic units is not managed.
  void compute_nb_of_affected_genes(const Mutation * mut, int & nb_genes_at_breakpoints, int & nb_genes_in_segment, int & nb_genes_in_replaced_segment);

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
  ExpManager* _exp_m;

  Individual*  _indiv;
  Dna*         _dna;
  Fuzzy*   _activ_contribution;
  Fuzzy*   _inhib_contribution;
  Fuzzy*   _phenotypic_contribution;
  // NB : _phenotypic_contribution is only an indicative value, not used for the whole phenotype computation

  // _rna_list always has 2 elements: make it an std::array
  Promoters2Strands _rna_list = {{},{}};
  // _protein_list always has 2 elements: make it an std::array
  #ifndef __REGUL
  std::array<std::list<Protein>, 2> _protein_list; // = {{},{}};
  #else
  std::array<std::list<Protein_R>, 2> _protein_list; // = {{},{}};
  #endif

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
#endif // AEVOL_GENETIC_UNIT_H__
