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


#ifndef AEVOL_STAT_RECORD_H_
#define AEVOL_STAT_RECORD_H_


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <list>



// =================================================================
//                            Project Files
// =================================================================


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;
class Stats;
class Individual;
class ae_population;

enum indiv_or_pop
{
  INDIV           = 0,
  POP             = 1,
  STDEVS          = 2,
  SKEWNESS        = 3,
  NB_INDIV_OR_POP = 4
};

enum chrom_or_gen_unit
{
  ALL_GU          = 0,
  CHROM           = 1,
  PLASMIDS        = 2,
  NB_CHROM_OR_GU  = 3
};

enum best_or_glob
{
  BEST            = 0,
  GLOB            = 1,
  SDEV            = 2,
  SKEW            = 3,
  NB_BEST_OR_GLOB = 4
};

enum stats_type
{
  FITNESS_STATS   = 0,
  MUTATION_STATS  = 1,
  GENES_STATS     = 2,
  BP_STATS        = 3,
  REAR_STATS      = 4,
  NB_STATS_TYPES  = 5
};






class StatRecord
{
  friend class Stats;
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    StatRecord() = delete;
    StatRecord(ExpManager * exp_m);
    StatRecord(const StatRecord &model);
    StatRecord(ExpManager * exp_m,
                   Individual * indiv,
                   chrom_or_gen_unit chrom_or_gu = CHROM,
                   bool compute_non_coding = true);
    StatRecord(ExpManager * exp_m,
                   const std::list<Individual *> indivs,
                   chrom_or_gen_unit chrom_or_gu = CHROM);
    StatRecord(ExpManager * exp_m,
                   const std::list<Individual *> indivs,
                   const StatRecord * means,
                   chrom_or_gen_unit chrom_or_gu = CHROM);
    StatRecord(ExpManager * exp_m,
                   const std::list<Individual *> indivs,
                   const StatRecord * means,
                   const StatRecord * stdevs,
                   chrom_or_gen_unit chrom_or_gu = CHROM);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~StatRecord();

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void initialize_data();
    void write_to_file( FILE* stat_file, stats_type stat_type_to_print ) const;
    
    void divide( double divisor );
    void divide_record( StatRecord const * means, double power );

    void add( StatRecord * to_add, int32_t index );
    void substract_power( StatRecord const * means, StatRecord const * to_substract, double power );
    
    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ExpManager * exp_m_;
    
    // NB : All the attributes are doubles because they will be used to
    //      compute averages over the population.
    indiv_or_pop record_type_;
    
    int32_t pop_size_;

    double  fitness_;

    double  metabolic_error_;
    double  parent_metabolic_error_;
    double  metabolic_fitness_;

    double  secretion_error_;
    double  parent_secretion_error_;
    double  secretion_fitness_;
  
    double  compound_amount_;
    
    int32_t  amount_of_dna_;
    int32_t  nb_coding_rnas_;
    int32_t  nb_non_coding_rnas_;
    double  av_size_coding_rnas_;      // NOT including promoter but including terminator
    double  av_size_non_coding_rnas_;  // NOT including promoter but including terminator
    int32_t  nb_functional_genes_;
    int32_t  nb_non_functional_genes_;
    double  av_size_functional_gene_;     // NOT including START and STOP codons
    double  av_size_non_functional_gene_; // NOT including START and STOP codons

    int32_t  nb_mut_;
    int32_t  nb_rear_;
    int32_t  nb_switch_;
    int32_t  nb_indels_;
    int32_t  nb_dupl_;
    int32_t  nb_del_;
    int32_t  nb_trans_;
    int32_t  nb_inv_;
    
    double  dupl_rate_;
    double  del_rate_;
    double  trans_rate_;
    double  inv_rate_;
    double  mean_align_score_;
    
    int32_t  nb_bases_in_0_CDS_;
    int32_t  nb_bases_in_0_functional_CDS_;
    int32_t  nb_bases_in_0_non_functional_CDS_;
    int32_t  nb_bases_in_0_RNA_;
    int32_t  nb_bases_in_0_coding_RNA_;
    int32_t  nb_bases_in_0_non_coding_RNA_;
    
    int32_t  nb_bases_non_essential_;
    int32_t  nb_bases_non_essential_including_nf_genes_;
    
    #ifdef __REGUL
      int32_t  nb_influences_;
      int32_t  nb_enhancing_influences_;
      int32_t  nb_operating_influences_;
      double  av_value_influences_;
      double  av_value_enhancing_influences_;
      double  av_value_operating_influences_;
    #endif
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_STAT_RECORD_H_
