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


#ifndef __AE_STAT_RECORD_H__
#define __AE_STAT_RECORD_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;
class ae_stats;
class ae_individual;
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






class ae_stat_record : public ae_object
{
  friend class ae_stats;
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_stat_record( ae_exp_manager* exp_m );
    ae_stat_record( const ae_stat_record &model );
    ae_stat_record( ae_exp_manager* exp_m, ae_individual const * indiv, chrom_or_gen_unit chrom_or_gu = CHROM, bool compute_non_coding = true, int32_t num_gener = -1 );
    ae_stat_record( ae_exp_manager* exp_m, ae_population const * pop, chrom_or_gen_unit chrom_or_gu = CHROM );
    ae_stat_record( ae_exp_manager* exp_m, ae_population const * pop, ae_stat_record const * means, chrom_or_gen_unit chrom_or_gu = CHROM );
    ae_stat_record( ae_exp_manager* exp_m, ae_population const * pop, ae_stat_record const * means, ae_stat_record const * stdevs, chrom_or_gen_unit chrom_or_gu = CHROM );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_stat_record( void );

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void initialize_data( void );
    void write_to_file( FILE* stat_file, stats_type stat_type_to_print ) const;
    
    void divide( double divisor );
    void divide_record( ae_stat_record const * means, double power );

    void add( ae_stat_record* to_add );
    void substract_power( ae_stat_record const * means, ae_stat_record const * to_substract, double power );
    
    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_stat_record( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    /*ae_stat_record( const ae_stat_record &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;
    
    // NB : All the attributes are doubles because they will be used to
    //      compute averages over the population.
    indiv_or_pop _record_type;
    
    double  _num_gener;
    double  _pop_size;

    double  _fitness;

    double  _metabolic_error;
    double  _parent_metabolic_error;
    double  _metabolic_fitness;

    double  _secretion_error;
    double  _parent_secretion_error;
    double  _secretion_fitness;
  
    double  _compound_amount;
    
    double  _amount_of_dna;
    double  _nb_coding_rnas;
    double  _nb_non_coding_rnas;
    double  _av_size_coding_rnas;      // NOT including promoter but including terminator
    double  _av_size_non_coding_rnas;  // NOT including promoter but including terminator
    double  _nb_functional_genes;
    double  _nb_non_functional_genes;
    double  _av_size_functional_gene;     // NOT including START and STOP codons
    double  _av_size_non_functional_gene; // NOT including START and STOP codons

    double  _nb_mut;
    double  _nb_rear;
    double  _nb_switch;
    double  _nb_indels;
    double  _nb_dupl;
    double  _nb_del;
    double  _nb_trans;
    double  _nb_inv;
    
    double  _dupl_rate;
    double  _del_rate;
    double  _trans_rate;
    double  _inv_rate;
    double  _mean_align_score;
    
    double  _nb_bases_in_0_CDS;
    double  _nb_bases_in_0_functional_CDS;
    double  _nb_bases_in_0_non_functional_CDS;
    double  _nb_bases_in_0_RNA;
    double  _nb_bases_in_0_coding_RNA;
    double  _nb_bases_in_0_non_coding_RNA;
    
    double  _nb_bases_non_essential;
    double  _nb_bases_non_essential_including_nf_genes;
    
    #ifdef __REGUL
      double  _nb_influences;
      double  _nb_enhancing_influences;
      double  _nb_operating_influences;
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


#endif // __AE_STAT_RECORD_H__
