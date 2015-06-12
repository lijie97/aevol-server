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


#ifndef AEVOL_STAT_RECORD_H__
#define AEVOL_STAT_RECORD_H__


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
    StatRecord(void) = delete;
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
    virtual ~StatRecord(void);

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void initialize_data( void );
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
    ExpManager * _exp_m;
    
    // NB : All the attributes are doubles because they will be used to
    //      compute averages over the population.
    indiv_or_pop _record_type;
    
    int32_t _pop_size;

    double  _fitness;

    double  _metabolic_error;
    double  _parent_metabolic_error;
    double  _metabolic_fitness;

    double  _secretion_error;
    double  _parent_secretion_error;
    double  _secretion_fitness;
  
    double  _compound_amount;
    
    int32_t  _amount_of_dna;
    int32_t  _nb_coding_rnas;
    int32_t  _nb_non_coding_rnas;
    double  _av_size_coding_rnas;      // NOT including promoter but including terminator
    double  _av_size_non_coding_rnas;  // NOT including promoter but including terminator
    int32_t  _nb_functional_genes;
    int32_t  _nb_non_functional_genes;
    double  _av_size_functional_gene;     // NOT including START and STOP codons
    double  _av_size_non_functional_gene; // NOT including START and STOP codons

    int32_t  _nb_mut;
    int32_t  _nb_rear;
    int32_t  _nb_switch;
    int32_t  _nb_indels;
    int32_t  _nb_dupl;
    int32_t  _nb_del;
    int32_t  _nb_trans;
    int32_t  _nb_inv;
    
    double  _dupl_rate;
    double  _del_rate;
    double  _trans_rate;
    double  _inv_rate;
    double  _mean_align_score;
    
    int32_t  _nb_bases_in_0_CDS;
    int32_t  _nb_bases_in_0_functional_CDS;
    int32_t  _nb_bases_in_0_non_functional_CDS;
    int32_t  _nb_bases_in_0_RNA;
    int32_t  _nb_bases_in_0_coding_RNA;
    int32_t  _nb_bases_in_0_non_coding_RNA;
    
    int32_t  _nb_bases_non_essential;
    int32_t  _nb_bases_non_essential_including_nf_genes;
    
    #ifdef __REGUL
      int32_t  _nb_influences;
      int32_t  _nb_enhancing_influences;
      int32_t  _nb_operating_influences;
      double  _av_value_influences;
      double  _av_value_enhancing_influences;
      double  _av_value_operating_influences;
    #endif
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_STAT_RECORD_H__
