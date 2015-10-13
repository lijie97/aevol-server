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




// =================================================================
//                              Libraries
// =================================================================
#include <list>



// =================================================================
//                            Project Files
// =================================================================
#include "StatRecord.h"

#include "ExpManager.h"
#include "ExpSetup.h"
#include "Individual.h"
#include "GeneticUnit.h"
#include "ReplicationReport.h"
#include "DnaReplicationReport.h"


using std::list;


namespace aevol {



//##############################################################################
//                                                                             #
//                             Class StatRecord                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
StatRecord::StatRecord(ExpManager * exp_m)
{
  _exp_m = exp_m;
  initialize_data();
}

StatRecord::StatRecord(const StatRecord &model)
{
  _exp_m = model._exp_m;
  
  _pop_size  = model._pop_size;
  
  _metabolic_error         = model._metabolic_error;
  _metabolic_fitness       = model._metabolic_fitness;
  _parent_metabolic_error  = model._parent_metabolic_error;
  
  _secretion_error         = model._secretion_error;
  _secretion_fitness       = model._secretion_fitness;
  _parent_secretion_error  = model._parent_secretion_error;
  
  _compound_amount   = model._compound_amount;
  
  _fitness = model._fitness;
  
  _amount_of_dna                = model._amount_of_dna;
  _nb_coding_rnas               = model._nb_coding_rnas;
  _nb_non_coding_rnas           = model._nb_non_coding_rnas;
  _av_size_coding_rnas          = model._av_size_coding_rnas;
  _av_size_non_coding_rnas      = model._av_size_non_coding_rnas;
  _nb_functional_genes          = model._nb_functional_genes;
  _nb_non_functional_genes      = model._nb_non_functional_genes;
  _av_size_functional_gene      = model._av_size_functional_gene;
  _av_size_non_functional_gene  = model._av_size_non_functional_gene;

  _nb_mut    = model._nb_mut;
  _nb_rear   = model._nb_rear;
  _nb_switch = model._nb_switch;
  _nb_indels = model._nb_indels;
  _nb_dupl   = model._nb_dupl;
  _nb_del    = model._nb_del;
  _nb_trans  = model._nb_trans;
  _nb_inv    = model._nb_inv;
  
  _dupl_rate        = model._dupl_rate;
  _del_rate         = model._del_rate;
  _trans_rate       = model._trans_rate;
  _inv_rate         = model._inv_rate;
  _mean_align_score = model._mean_align_score;
  
  _nb_bases_in_0_CDS                = model._nb_bases_in_0_CDS;
  _nb_bases_in_0_functional_CDS     = model._nb_bases_in_0_functional_CDS;
  _nb_bases_in_0_non_functional_CDS = model._nb_bases_in_0_non_functional_CDS;
  _nb_bases_in_0_RNA                = model._nb_bases_in_0_RNA;
  _nb_bases_in_0_coding_RNA         = model._nb_bases_in_0_coding_RNA;
  _nb_bases_in_0_non_coding_RNA     = model._nb_bases_in_0_non_coding_RNA;
      
  _nb_bases_non_essential                     = model._nb_bases_non_essential;
  _nb_bases_non_essential_including_nf_genes  = model._nb_bases_non_essential_including_nf_genes;
    
  #ifdef __REGUL
    _nb_influences                 = model._nb_influences;
    _nb_enhancing_influences       = model._nb_enhancing_influences;
    _nb_operating_influences       = model._nb_operating_influences;
    _av_value_influences           = model._av_value_influences;
    _av_value_enhancing_influences = model._av_value_enhancing_influences;
    _av_value_operating_influences = model._av_value_operating_influences;
  #endif
}

StatRecord::StatRecord(ExpManager* exp_m,
                       Individual* indiv,
                       chrom_or_gen_unit chrom_or_gu,
                       bool compute_non_coding)
{
  _exp_m = exp_m;
  initialize_data();
  _record_type = INDIV;
    
  // ---------------
  // Simulation data
  // ---------------
  _pop_size = 0; // The pop_size value is irrelevant when dealing with a single
                 // individual. It is present for column alignment.
  
  #ifdef __REGUL
    // TODO
  #endif  
    
  // TODO : These conditions are not well managed!!!
  if (indiv->get_nb_genetic_units() == 1)
  { // One single Genetic Unit
    // -------------------------------------------------
    // Compute statistical data for the given individual
    // -------------------------------------------------
    ReplicationReport* replic_report = nullptr;
    if (_exp_m->get_tree() != nullptr)
      replic_report = _exp_m->get_tree()->get_report_by_index(AeTime::get_time(),
                                                              indiv->get_id());

    if (compute_non_coding)
      indiv->compute_non_coding();

    const GeneticUnit& gen_unit = *indiv->get_genetic_unit_list().begin();

    // Metabolic error stats
    _metabolic_error = indiv->get_dist_to_target_by_feature(METABOLISM);
    _metabolic_fitness = indiv->get_fitness_by_feature(METABOLISM);
    _parent_metabolic_error = (replic_report != NULL) ?
                              replic_report->get_parent_metabolic_error() :
                              0.0;

    // Fitness
    _fitness = indiv->get_fitness();

    // Secretion stats
    if (_exp_m->get_with_secretion())
    {
       _secretion_error   = indiv->get_dist_to_target_by_feature(SECRETION);
       _secretion_fitness = indiv->get_fitness_by_feature(SECRETION);
       _compound_amount   = indiv->get_grid_cell()->compound_amount();
       _parent_secretion_error = 0.0;

      if (replic_report != NULL)
      {
        _parent_secretion_error = replic_report->get_parent_secretion_error();
      }
    }
    else
    {
      _secretion_error   = 0.0;
      _secretion_fitness = 0.0;
      _compound_amount   = 0.0;
      _parent_secretion_error = 0.0;
    }

    // Genes and RNA stats
    _amount_of_dna               = gen_unit.get_dna()->length();
    _nb_coding_rnas              = gen_unit.get_nb_coding_RNAs();
    _nb_non_coding_rnas          = gen_unit.get_nb_non_coding_RNAs();
    _av_size_coding_rnas         = gen_unit.get_av_size_coding_RNAs();
    _av_size_non_coding_rnas     = gen_unit.get_av_size_non_coding_RNAs();
    _nb_functional_genes         = gen_unit.get_nb_functional_genes();
    _nb_non_functional_genes     = gen_unit.get_nb_non_functional_genes();
    _av_size_functional_gene     = gen_unit.get_av_size_functional_genes();
    _av_size_non_functional_gene = gen_unit.get_av_size_non_functional_genes();

    // Non coding stats
    if (compute_non_coding) {
      _nb_bases_in_0_CDS                = gen_unit.get_nb_bases_in_0_CDS();
      _nb_bases_in_0_functional_CDS     = gen_unit.get_nb_bases_in_0_functional_CDS();
      _nb_bases_in_0_non_functional_CDS = gen_unit.get_nb_bases_in_0_non_functional_CDS();
      _nb_bases_in_0_RNA                = gen_unit.get_nb_bases_in_0_RNA();
      _nb_bases_in_0_coding_RNA         = gen_unit.get_nb_bases_in_0_coding_RNA();
      _nb_bases_in_0_non_coding_RNA     = gen_unit.get_nb_bases_in_0_non_coding_RNA();

      _nb_bases_non_essential                     = gen_unit.get_nb_bases_non_essential();
      _nb_bases_non_essential_including_nf_genes  = gen_unit.get_nb_bases_non_essential_including_nf_genes();
    }

    // Mutation stats
    if (replic_report != NULL)
    {
      _nb_mut    = replic_report->get_nb(S_MUT);
      _nb_rear   = replic_report->get_nb(REARR);
      _nb_switch = replic_report->get_nb(SWITCH);
      _nb_indels = replic_report->get_nb(INDEL);
      _nb_dupl   = replic_report->get_nb(DUPL);
      _nb_del    = replic_report->get_nb(DEL);
      _nb_trans  = replic_report->get_nb(TRANS);
      _nb_inv    = replic_report->get_nb(INV);

      // Rearrangement rate stats
      int32_t parent_genome_size = replic_report->get_parent_genome_size();
      _dupl_rate  = _nb_dupl  / parent_genome_size;
      _del_rate   = _nb_del   / parent_genome_size;
      _trans_rate = _nb_trans / parent_genome_size;
      _inv_rate   = _nb_inv   / parent_genome_size;

      //~ // <DEBUG>
      //~ if (_nb_dupl + _nb_del + _nb_trans + _nb_inv != 0)
      //~ {
        //~ printf("_nb_dupl : %"PRId32"\n_nb_del : %"PRId32"\n_nb_trans : %"PRId32"\n_nb_inv : %"PRId32"\n",
                //~ (int32_t) _nb_dupl, (int32_t) _nb_del, (int32_t) _nb_trans, (int32_t) _nb_inv);
        //~ printf("parent genome size : %"PRId32"\n", parent_genome_size);
        //~ printf("_dupl_rate : %f\n_del_rate : %f\n_trans_rate : %f\n_inv_rate : %f\n",
                //~ _dupl_rate, _del_rate, _trans_rate, _inv_rate);
        //~ getchar();
      //~ }
      //~ // </DEBUG>

      _mean_align_score = replic_report->get_mean_align_score();
    }
  }
  else if (chrom_or_gu == ALL_GU)
  {
    // -------------------------------------------------
    // Compute statistical data for the given individual
    // -------------------------------------------------
    ReplicationReport* replic_report = nullptr;
    if (_exp_m->get_tree() != nullptr)
      replic_report = _exp_m->get_tree()->get_report_by_index(AeTime::get_time(),
                                                              indiv->get_id());

    // Metabolic error stats
    _metabolic_error = (double) indiv->get_dist_to_target_by_feature(METABOLISM);
    _metabolic_fitness = (double) indiv->get_fitness_by_feature(METABOLISM);
    _parent_metabolic_error = (replic_report != NULL) ? replic_report->get_parent_metabolic_error() : 0.0;

    // Fitness
    _fitness = indiv->get_fitness();

    // Secretion stats
    if (_exp_m->get_with_secretion())
    {
       _secretion_error = (double) indiv->get_dist_to_target_by_feature(SECRETION);
       _secretion_fitness = (double) indiv->get_fitness_by_feature(SECRETION);
       _compound_amount   = (double) indiv->get_grid_cell()->compound_amount();
       _parent_secretion_error = 0.0;

      if (replic_report != NULL)
      {
        _parent_secretion_error = replic_report->get_parent_secretion_error();
      }
    }
    else
    {
      _secretion_error   = 0.0;
      _secretion_fitness = 0.0;
      _compound_amount   = 0.0;
      _parent_secretion_error = 0.0;
    }

    for (auto& gen_unit: indiv->get_genetic_unit_list_nonconst()) {
      // Genes and RNA stats
      _amount_of_dna               += gen_unit.get_dna()->length();
      _nb_coding_rnas              += gen_unit.get_nb_coding_RNAs();
      _nb_non_coding_rnas          += gen_unit.get_nb_non_coding_RNAs();
      _av_size_coding_rnas         += gen_unit.get_av_size_coding_RNAs();
      _av_size_non_coding_rnas     += gen_unit.get_av_size_non_coding_RNAs();
      _nb_functional_genes         += gen_unit.get_nb_functional_genes();
      _nb_non_functional_genes     += gen_unit.get_nb_non_functional_genes();
      _av_size_functional_gene     += gen_unit.get_av_size_functional_genes();
      _av_size_non_functional_gene += gen_unit.get_av_size_non_functional_genes();

      // Non coding stats
      if (compute_non_coding)
      {
        _nb_bases_in_0_CDS                += gen_unit.get_nb_bases_in_0_CDS();
        _nb_bases_in_0_functional_CDS     += gen_unit.get_nb_bases_in_0_functional_CDS();
        _nb_bases_in_0_non_functional_CDS += gen_unit.get_nb_bases_in_0_non_functional_CDS();
        _nb_bases_in_0_RNA                += gen_unit.get_nb_bases_in_0_RNA();
        _nb_bases_in_0_coding_RNA         += gen_unit.get_nb_bases_in_0_coding_RNA();
        _nb_bases_in_0_non_coding_RNA     += gen_unit.get_nb_bases_in_0_non_coding_RNA();
        
        _nb_bases_non_essential                     += gen_unit.get_nb_bases_non_essential();
        _nb_bases_non_essential_including_nf_genes  += gen_unit.get_nb_bases_non_essential_including_nf_genes();
      }

      // Mutation stats
      if (replic_report != NULL)
      {
        _nb_mut    += replic_report->get_nb(S_MUT);
        _nb_rear   += replic_report->get_nb(REARR);
        _nb_switch += replic_report->get_nb(SWITCH);
        _nb_indels += replic_report->get_nb(INDEL);
        _nb_dupl   += replic_report->get_nb(DUPL);
        _nb_del    += replic_report->get_nb(DEL);
        _nb_trans  += replic_report->get_nb(TRANS);
        _nb_inv    += replic_report->get_nb(INV);
      }
    }
    
    // Rearrangement rate stats
    if (replic_report != NULL)
    {
      int32_t parent_genome_size = replic_report->get_parent_genome_size();
      _dupl_rate  = _nb_dupl  / parent_genome_size;
      _del_rate   = _nb_del   / parent_genome_size;
      _trans_rate = _nb_trans / parent_genome_size;
      _inv_rate   = _nb_inv   / parent_genome_size;
      _mean_align_score = replic_report->get_mean_align_score();
    }
  }
  else // => We have a multi-GU individual and we want only the main chromosome or only the plasmids
  // WARNING (TODO) As it is coded, this will work only if there is ONE SINGLE PLASMID!
  {
    if (chrom_or_gu != PLASMIDS and chrom_or_gu != CHROM) {
      printf("%s: error: StatRecord called with inappropriate `chrom_or_gu`\n", __FILE__);
      exit(EXIT_FAILURE);
    }

    GeneticUnit& gen_unit = (chrom_or_gu == PLASMIDS) ?
        *std::next(indiv->get_genetic_unit_list_nonconst().begin()) :
        indiv->get_genetic_unit_list_nonconst().front();

    // -------------------------------------------------
    // Compute statistical data for the given individual
    // -------------------------------------------------
    ReplicationReport* replic_report = nullptr;
    if (_exp_m->get_tree() != nullptr)
      replic_report = _exp_m->get_tree()->get_report_by_index(AeTime::get_time(),
                                                              indiv->get_id());
    
    // Metabolic error stats
    _metabolic_error = (double) gen_unit.get_dist_to_target_by_feature(METABOLISM);
    _metabolic_fitness = (double) gen_unit.get_fitness_by_feature(METABOLISM);
    _parent_metabolic_error = (replic_report != NULL) ? replic_report->get_parent_metabolic_error() : 0.0;
    
    // Fitness
    _fitness = indiv->get_fitness();
  
    // Secretion stats
    if (_exp_m->get_with_secretion())
    {
       _secretion_error = (double) gen_unit.get_dist_to_target_by_feature(SECRETION);
       _secretion_fitness = (double) gen_unit.get_fitness_by_feature(SECRETION);
       _compound_amount   = (double) indiv->get_grid_cell()->compound_amount();
       _parent_secretion_error = 0.0;
  
      if (replic_report != NULL)
      {
        _parent_secretion_error = replic_report->get_parent_secretion_error();
      }
    }
    else
    {
      _secretion_error   = 0.0;
      _secretion_fitness = 0.0;
      _compound_amount   = 0.0;
      _parent_secretion_error = 0.0;
    }
    
      // Genes and RNA stats
    _amount_of_dna               = gen_unit.get_dna()->length();
    _nb_coding_rnas              = gen_unit.get_nb_coding_RNAs();
    _nb_non_coding_rnas          = gen_unit.get_nb_non_coding_RNAs();
    _av_size_coding_rnas         = gen_unit.get_av_size_coding_RNAs();
    _av_size_non_coding_rnas     = gen_unit.get_av_size_non_coding_RNAs();
    _nb_functional_genes         = gen_unit.get_nb_functional_genes();
    _nb_non_functional_genes     = gen_unit.get_nb_non_functional_genes();
    _av_size_functional_gene     = gen_unit.get_av_size_functional_genes();
    _av_size_non_functional_gene = gen_unit.get_av_size_non_functional_genes();
    
      // Non coding stats
    if (compute_non_coding)
    {
      _nb_bases_in_0_CDS                  = gen_unit.get_nb_bases_in_0_CDS();
      _nb_bases_in_0_functional_CDS       = gen_unit.get_nb_bases_in_0_functional_CDS();
      _nb_bases_in_0_non_functional_CDS   = gen_unit.get_nb_bases_in_0_non_functional_CDS();
      _nb_bases_in_0_RNA                  = gen_unit.get_nb_bases_in_0_RNA();
      _nb_bases_in_0_coding_RNA           = gen_unit.get_nb_bases_in_0_coding_RNA();
      _nb_bases_in_0_non_coding_RNA       = gen_unit.get_nb_bases_in_0_non_coding_RNA();
      
      _nb_bases_non_essential                     = gen_unit.get_nb_bases_non_essential();
      _nb_bases_non_essential_including_nf_genes  = gen_unit.get_nb_bases_non_essential_including_nf_genes();
    }
    
    // Mutation stats
    // TODO <david.parsons@inria.fr> Disabled
//    if (gen_unit.get_dna()->get_replication_report() != NULL)
//    {
//      _nb_mut    = gen_unit.get_dna()->get_replication_report()->get_nb(S_MUT);
//      _nb_rear   = gen_unit.get_dna()->get_replication_report()->get_nb(REARR);
//      _nb_switch = gen_unit.get_dna()->get_replication_report()->get_nb(SWITCH);
//      _nb_indels = gen_unit.get_dna()->get_replication_report()->get_nb(INDEL);
//      _nb_dupl   = gen_unit.get_dna()->get_replication_report()->get_nb(DUPL);
//      _nb_del    = gen_unit.get_dna()->get_replication_report()->get_nb(DEL);
//      _nb_trans  = gen_unit.get_dna()->get_replication_report()->get_nb(TRANS);
//      _nb_inv    = gen_unit.get_dna()->get_replication_report()->get_nb(INV);
//    }
    
    // Rearrangement rate stats
    if (replic_report != NULL)
    {
      int32_t parent_genome_size = replic_report->get_parent_genome_size();
      _dupl_rate  = _nb_dupl  / parent_genome_size;
      _del_rate   = _nb_del   / parent_genome_size;
      _trans_rate = _nb_trans / parent_genome_size;
      _inv_rate   = _nb_inv   / parent_genome_size;
      _mean_align_score = replic_report->get_mean_align_score();
    }
  }
}

// Calculate average statistics for all the recorded values 
StatRecord::StatRecord(ExpManager * exp_m,
                               const list<Individual *> indivs,
                               chrom_or_gen_unit chrom_or_gu)
{
  _exp_m = exp_m;
  initialize_data();
  
  _record_type = POP;
  
  // ---------------
  // Simulation data
  // ---------------
  _pop_size = (double) indivs.size();

  // ------------------------------------------------------------------
  // Compute statistical data for the each individual in the population
  // ------------------------------------------------------------------
  for (const auto& indiv: indivs)
  {
    StatRecord * indiv_stat_record = new StatRecord(_exp_m, indiv, chrom_or_gu, false);
    this->add(indiv_stat_record, indiv->get_id());
    delete indiv_stat_record;
  }
  
  
  // ------------------------------------------------------------------
  // Divide every accumulator by the number of indivs in the population
  // ------------------------------------------------------------------
  this->divide(_pop_size);
}

// Calculate standard deviation for all the recorded values 
StatRecord::StatRecord(ExpManager * exp_m,
                               const list<Individual *> indivs,
                               const StatRecord * means,
                               chrom_or_gen_unit chrom_or_gu)
{
  _exp_m = exp_m;
  initialize_data();
  
  _record_type = STDEVS;
  
  // ---------------
  // Simulation data
  // ---------------
  _pop_size = indivs.size();

  // ------------------------------------------------------------------
  // Compute statistical data for the each individual in the population
  // ------------------------------------------------------------------
  for (const auto& indiv : indivs) {
    StatRecord * indiv_stat_record = new StatRecord(_exp_m, indiv, chrom_or_gu, false);
    this->substract_power(means, indiv_stat_record, 2);
    delete indiv_stat_record;
  }
  
  // ---------------------------------------------------------------------------------
  // Divide every accumulator by the square root of number of indivs in the population
  // ---------------------------------------------------------------------------------
  this->divide(pow((_pop_size-1), 0.5));
}

 // Calculate skewness for all the recorded values 
StatRecord::StatRecord(ExpManager * exp_m,
                               const list<Individual *> indivs,
                               const StatRecord * means,
                               const StatRecord * stdevs,
                               chrom_or_gen_unit chrom_or_gu)
{
  _exp_m = exp_m;
  initialize_data();
  
  _record_type = SKEWNESS;
  
  // ---------------
  // Simulation data
  // ---------------
  _pop_size = (double) indivs.size();

  // ------------------------------------------------------------------
  // Compute statistical data for the each individual in the population
  // ------------------------------------------------------------------
  for (const auto& indiv : indivs)
  {
    StatRecord * indiv_stat_record = new StatRecord(_exp_m, indiv, chrom_or_gu, false);
    this->substract_power(means, indiv_stat_record, 3);
    delete indiv_stat_record;
  }

  this->divide(-_pop_size);
  
  this->divide_record(stdevs, 3/2);
}

// =================================================================
//                             Destructors
// =================================================================
StatRecord::~StatRecord(void)
{
}

// =================================================================
//                            Public Methods
// =================================================================
void StatRecord::initialize_data(void)
{
  _pop_size  = 0.0;
  
  _metabolic_error         = 0.0;
  _metabolic_fitness       = 0.0;
  _parent_metabolic_error  = 0.0;

  _secretion_error         = 0.0;
  _parent_secretion_error  = 0.0;
  
  _secretion_fitness = 0.0;
  _compound_amount   = 0.0;
  
  _fitness = 0.0;
  
  _amount_of_dna               = 0.0;
  _nb_coding_rnas              = 0.0;
  _nb_non_coding_rnas          = 0.0;
  _av_size_coding_rnas         = 0.0;
  _av_size_non_coding_rnas     = 0.0;
  _nb_functional_genes         = 0.0;
  _nb_non_functional_genes     = 0.0;
  _av_size_functional_gene     = 0.0;
  _av_size_non_functional_gene = 0.0;

  _nb_mut    = 0.0;
  _nb_rear   = 0.0;
  _nb_switch = 0.0;
  _nb_indels = 0.0;
  _nb_dupl   = 0.0;
  _nb_del    = 0.0;
  _nb_trans  = 0.0;
  _nb_inv    = 0.0;
  
  _dupl_rate  = 0.0;
  _del_rate   = 0.0;
  _trans_rate = 0.0;
  _inv_rate   = 0.0;
  _mean_align_score = 0.0;
  
  _nb_bases_in_0_CDS                = 0.0;
  _nb_bases_in_0_functional_CDS     = 0.0;
  _nb_bases_in_0_non_functional_CDS = 0.0;
  _nb_bases_in_0_RNA                = 0.0;
  _nb_bases_in_0_coding_RNA         = 0.0;
  _nb_bases_in_0_non_coding_RNA     = 0.0;
    
  _nb_bases_non_essential                     = 0.0;
  _nb_bases_non_essential_including_nf_genes  = 0.0;
    
  #ifdef __REGUL
    _nb_influences                 = 0.0;
    _nb_enhancing_influences       = 0.0;
    _nb_operating_influences       = 0.0;
    _av_value_influences           = 0.0;
    _av_value_enhancing_influences = 0.0;
    _av_value_operating_influences = 0.0;
  #endif
}

void StatRecord::write_to_file(FILE* stat_file, stats_type stat_type_to_print) const
{
  if (_record_type == INDIV)
  {
    if (stat_type_to_print == FITNESS_STATS)
    {
      fprintf(stat_file,
          "%" PRId64 " %" PRId32 " %e %" PRId32 " %e %e %e %e %e %e %e",
          AeTime::get_time(),
          _pop_size,
          _fitness,
          _amount_of_dna,
          _metabolic_error,
          _parent_metabolic_error,
          _metabolic_fitness,
          _secretion_error,
          _parent_secretion_error,
          _secretion_fitness,
          _compound_amount);
      #ifdef __REGUL
        fprintf(stat_file,
            " %" PRId32 " %" PRId32 " %" PRId32 " %f %f %f",
            _nb_influences,
            _nb_enhancing_influences,
            _nb_operating_influences,
            _av_value_influences,
            _av_value_enhancing_influences,
            _av_value_operating_influences);
      #endif
    }
    if (stat_type_to_print == MUTATION_STATS)
    {
      fprintf(stat_file,
          "%" PRId64 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 "",
          AeTime::get_time(),
          _nb_mut,
          _nb_rear,
          _nb_switch,
          _nb_indels,
          _nb_dupl,
          _nb_del,
          _nb_trans,
          _nb_inv);
    }
    if (stat_type_to_print == GENES_STATS)
    {
      fprintf(stat_file,
          "%" PRId64 " %" PRId32 " %" PRId32 " %f %f %" PRId32 " %" PRId32 " %f %f ",
          AeTime::get_time(),
          _nb_coding_rnas,
          _nb_non_coding_rnas,
          _av_size_coding_rnas,
          _av_size_non_coding_rnas,
          _nb_functional_genes,
          _nb_non_functional_genes,
          _av_size_functional_gene,
          _av_size_non_functional_gene);
    }
    if (stat_type_to_print == BP_STATS)
    {
      fprintf(stat_file,
          "%" PRId64 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 "",
          AeTime::get_time(),
          _nb_bases_in_0_CDS,
          _nb_bases_in_0_functional_CDS,
          _nb_bases_in_0_non_functional_CDS,
          _nb_bases_in_0_RNA,
          _nb_bases_in_0_coding_RNA,
          _nb_bases_in_0_non_coding_RNA,
          _nb_bases_non_essential,
          _nb_bases_non_essential_including_nf_genes);
    }
    if (stat_type_to_print == REAR_STATS)
    {
      fprintf(stat_file,
          "%" PRId64 " %e %e %e %e %f",
          AeTime::get_time(),
          _dupl_rate,
          _del_rate,
          _trans_rate,
          _inv_rate,
          _mean_align_score);
    }
  }
  else // if _record_type == POP
  {
   if (stat_type_to_print == FITNESS_STATS)
    {
      fprintf(stat_file,
          "%" PRId64 " %" PRId32 " %e %" PRId32 " %e %e %e %e %e %e %e",
          AeTime::get_time(),
          _pop_size,
          _fitness,              
          _amount_of_dna, 
          _metabolic_error,
          _parent_metabolic_error,
          _metabolic_fitness,
          _secretion_error,
          _parent_secretion_error,
          _secretion_fitness,
          _compound_amount);

      #ifdef __REGUL
        fprintf(stat_file,
            " %" PRId32 " %" PRId32 " %" PRId32 " %f %f %f",
            _nb_influences,
            _nb_enhancing_influences,
            _nb_operating_influences,
            _av_value_influences,
            _av_value_enhancing_influences,
            _av_value_operating_influences);
      #endif
    }
    if (stat_type_to_print == MUTATION_STATS)
    {        
      fprintf(stat_file,
          "%" PRId64 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
          " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 "",
          AeTime::get_time(),
          _nb_mut,
          _nb_rear,
          _nb_switch,
          _nb_indels,
          _nb_dupl,
          _nb_del,
          _nb_trans,
          _nb_inv);

    }
    if (stat_type_to_print == GENES_STATS)
    {
      fprintf(stat_file,
          "%" PRId64 " %" PRId32 " %" PRId32 " %f %f %" PRId32 " %" PRId32
          " %f %f",
          AeTime::get_time(),
          _nb_coding_rnas,
          _nb_non_coding_rnas,
          _av_size_coding_rnas,
          _av_size_non_coding_rnas,
          _nb_functional_genes,
          _nb_non_functional_genes,
          _av_size_functional_gene,
          _av_size_non_functional_gene);
    }
    if (stat_type_to_print == BP_STATS)
    {
     // TO DO (if needed) : base-pair stats for all individuals, not just for the best one. 
     //
     // fprintf( stat_file, "%" PRId64 " %f %f %f %f %f %f %f %f",
     //         AeTime::get_time(),
     //         nb_bases_in_0_CDS_,
     //         nb_bases_in_0_functional_CDS_,
     //         nb_bases_in_0_non_functional_CDS_,
     //         nb_bases_in_0_RNA_,
     //         nb_bases_in_0_coding_RNA_,
     //         nb_bases_in_0_non_coding_RNA_,
     //         nb_bases_non_essential_,
     //         nb_bases_non_essential_including_nf_genes_);
    }
    if (stat_type_to_print == REAR_STATS)
    {
      fprintf(stat_file,
          "%" PRId64 " %e %e %e %e %f",
          AeTime::get_time(),
          _dupl_rate,
          _del_rate,
          _trans_rate,
          _inv_rate,
          _mean_align_score);
    }
  }
  
  fprintf(stat_file, "\n");
}

void StatRecord::divide(double divisor)
{
  // NB : pop_size is a "global" value and must not be divided.
  
  _fitness                 /= divisor;
  
  _metabolic_error         /= divisor;
  _parent_metabolic_error  /= divisor;
  _metabolic_fitness       /= divisor;

  _secretion_error         /= divisor;
  _parent_secretion_error  /= divisor;
  
  _secretion_fitness       /= divisor;
  _compound_amount         /= divisor;
  
  _amount_of_dna               /= divisor;
  _nb_coding_rnas              /= divisor;
  _nb_non_coding_rnas          /= divisor;
  _av_size_coding_rnas         /= divisor;
  _av_size_non_coding_rnas     /= divisor;
  _nb_functional_genes          /= divisor;
  _nb_non_functional_genes      /= divisor;
  _av_size_functional_gene      /= divisor;
  _av_size_non_functional_gene  /= divisor;

  _nb_mut    /= divisor;
  _nb_rear   /= divisor;
  _nb_switch /= divisor;
  _nb_indels /= divisor;
  _nb_dupl   /= divisor;
  _nb_del    /= divisor;
  _nb_trans  /= divisor;
  _nb_inv    /= divisor;
  
  //~ printf("PREFINAL %f %f %f %f\n", _dupl_rate, _del_rate, _trans_rate, _inv_rate);
  _dupl_rate  /= divisor;
  _del_rate   /= divisor;
  _trans_rate /= divisor;
  _inv_rate   /= divisor;
  //~ printf("FINAL %f %f %f %f\n", _dupl_rate, _del_rate, _trans_rate, _inv_rate);
  //~ getchar();
  _mean_align_score /= divisor;
  
  _nb_bases_in_0_CDS                /= divisor;
  _nb_bases_in_0_functional_CDS     /= divisor;
  _nb_bases_in_0_non_functional_CDS /= divisor;
  _nb_bases_in_0_RNA                /= divisor;
  _nb_bases_in_0_coding_RNA         /= divisor;
  _nb_bases_in_0_non_coding_RNA     /= divisor;
    
  _nb_bases_non_essential                     /= divisor;
  _nb_bases_non_essential_including_nf_genes  /= divisor;
    
  #ifdef __REGUL
    _nb_influences                 /= divisor;
    _nb_enhancing_influences       /= divisor;
    _nb_operating_influences       /= divisor;
    _av_value_influences           /= divisor;
    _av_value_enhancing_influences /= divisor;
    _av_value_operating_influences /= divisor;
  #endif
}


void StatRecord::divide_record(StatRecord const * to_divide, double power)
{
  // NB : pop_size is a "global" value and must not be divided.
  
  if (to_divide->_fitness != 0) { _fitness    /= pow(to_divide->_fitness, power); }
  
  if (to_divide->_metabolic_error != 0)        { _metabolic_error         /= pow(to_divide->_metabolic_error, power); }
  if (to_divide->_parent_metabolic_error != 0) { _parent_metabolic_error  /= pow(to_divide->_parent_metabolic_error, power); }
  if (to_divide->_metabolic_fitness != 0)        { _metabolic_fitness         /= pow(to_divide->_metabolic_fitness, power); }
  
  if (to_divide->_secretion_error != 0)        { _secretion_error         /= pow(to_divide->_secretion_error, power); }
  if (to_divide->_parent_secretion_error != 0) { _parent_secretion_error  /= pow(to_divide->_parent_secretion_error, power); }
  
  if (to_divide->_secretion_fitness != 0)       { _secretion_fitness       /= pow(to_divide->_secretion_fitness, power); }
  if (to_divide->_compound_amount != 0)        { _compound_amount         /= pow(to_divide->_compound_amount, power); }
  
  if (to_divide->_amount_of_dna != 0)               { _amount_of_dna               /= pow(to_divide->_amount_of_dna, power); }
  if (to_divide->_nb_coding_rnas != 0)              { _nb_coding_rnas              /= pow(to_divide->_nb_coding_rnas, power); }
  if (to_divide->_nb_non_coding_rnas != 0)          { _nb_non_coding_rnas          /= pow(to_divide->_nb_non_coding_rnas, power); }
  if (to_divide->_av_size_coding_rnas != 0)         { _av_size_coding_rnas         /= pow(to_divide->_av_size_coding_rnas, power); }
  if (to_divide->_av_size_non_coding_rnas != 0)     { _av_size_non_coding_rnas     /= pow(to_divide->_av_size_non_coding_rnas, power); }
  if (to_divide->_nb_functional_genes != 0)         { _nb_functional_genes         /= pow(to_divide->_nb_functional_genes, power); }
  if (to_divide->_nb_non_functional_genes != 0)     { _nb_non_functional_genes     /= pow(to_divide->_nb_non_functional_genes, power); }
  if (to_divide->_av_size_functional_gene != 0)     { _av_size_functional_gene     /= pow(to_divide->_av_size_functional_gene, power); }
  if (to_divide->_av_size_non_functional_gene != 0) { _av_size_non_functional_gene /= pow(to_divide->_av_size_non_functional_gene, power); }

  if (to_divide->_nb_mut != 0)     { _nb_mut    /= pow(to_divide->_nb_mut, power); }
  if (to_divide->_nb_rear != 0)    { _nb_rear   /= pow(to_divide->_nb_rear, power); }
  if (to_divide->_nb_switch != 0)  { _nb_switch /= pow(to_divide->_nb_switch, power); }
  if (to_divide->_nb_indels != 0)  { _nb_indels /= pow(to_divide->_nb_indels, power); }
  if (to_divide->_nb_dupl != 0)    { _nb_dupl   /= pow(to_divide->_nb_dupl, power); }
  if (to_divide->_nb_del != 0)     { _nb_del    /= pow(to_divide->_nb_del, power); }
  if (to_divide->_nb_trans != 0)   { _nb_trans  /= pow(to_divide->_nb_trans, power); }
  if (to_divide->_nb_inv != 0)     { _nb_inv    /= pow(to_divide->_nb_inv, power); }
  
  if (to_divide->_dupl_rate != 0)        { _dupl_rate  /= pow(to_divide->_dupl_rate, power); }
  if (to_divide->_del_rate != 0)         { _del_rate   /= pow(to_divide->_del_rate, power); }
  if (to_divide->_trans_rate != 0)       { _trans_rate /= pow(to_divide->_trans_rate, power); }
  if (to_divide->_inv_rate != 0)         { _inv_rate   /= pow(to_divide->_inv_rate, power); }
  if (to_divide->_mean_align_score != 0) { _mean_align_score /= pow(to_divide->_mean_align_score, power); }
  
  if (to_divide->_nb_bases_in_0_CDS != 0)                { _nb_bases_in_0_CDS                /= pow(to_divide->_nb_bases_in_0_CDS, power); }
  if (to_divide->_nb_bases_in_0_functional_CDS != 0)     { _nb_bases_in_0_functional_CDS     /= pow(to_divide->_nb_bases_in_0_functional_CDS, power); }
  if (to_divide->_nb_bases_in_0_non_functional_CDS != 0) { _nb_bases_in_0_non_functional_CDS /= pow(to_divide->_nb_bases_in_0_non_functional_CDS, power); }
  if (to_divide->_nb_bases_in_0_RNA != 0)                { _nb_bases_in_0_RNA                /= pow(to_divide->_nb_bases_in_0_RNA, power); }
  if (to_divide->_nb_bases_in_0_coding_RNA != 0)         { _nb_bases_in_0_coding_RNA         /= pow(to_divide->_nb_bases_in_0_coding_RNA, power); }
  if (to_divide->_nb_bases_in_0_non_coding_RNA != 0)     { _nb_bases_in_0_non_coding_RNA     /= pow(to_divide->_nb_bases_in_0_non_coding_RNA, power); }
    
  if (to_divide->_nb_bases_non_essential != 0)                    { _nb_bases_non_essential                     /= pow(to_divide->_nb_bases_non_essential, power); }
  if (to_divide->_nb_bases_non_essential_including_nf_genes != 0) { _nb_bases_non_essential_including_nf_genes  /= pow(to_divide->_nb_bases_non_essential_including_nf_genes, power); }
    
  #ifdef __REGUL
    if (to_divide->_nb_influences != 0)                 { _nb_influences                 /= pow(to_divide->_nb_influences, power); }
    if (to_divide->_nb_enhancing_influences != 0)       { _nb_enhancing_influences       /= pow(to_divide->_nb_enhancing_influences, power); }
    if (to_divide->_nb_operating_influences != 0)       { _nb_operating_influences       /= pow(to_divide->_nb_operating_influences, power); }
    if (to_divide->_av_value_influences != 0)           { _av_value_influences           /= pow(to_divide->_av_value_influences, power); }
    if (to_divide->_av_value_enhancing_influences != 0) { _av_value_enhancing_influences /= pow(to_divide->_av_value_enhancing_influences, power); }
    if (to_divide->_av_value_operating_influences != 0) { _av_value_operating_influences /= pow(to_divide->_av_value_operating_influences, power); }
  #endif
}

void StatRecord::add(StatRecord * to_add, int32_t index)
{
  // NB : pop_size is a global values and must not be summed.
  
  _fitness                 += to_add->_fitness;
  
  _metabolic_error         += to_add->_metabolic_error;
  _parent_metabolic_error  += to_add->_parent_metabolic_error;
  _metabolic_fitness       += to_add->_metabolic_fitness;

  _secretion_error         += to_add->_secretion_error;
  _parent_secretion_error  += to_add->_parent_secretion_error;
  
  _secretion_fitness       += to_add->_secretion_fitness;
  _compound_amount         += to_add->_compound_amount;
  
  _amount_of_dna               += to_add->_amount_of_dna;
  _nb_coding_rnas              += to_add->_nb_coding_rnas;
  _nb_non_coding_rnas          += to_add->_nb_non_coding_rnas;
  _av_size_coding_rnas         += to_add->_av_size_coding_rnas;
  _av_size_non_coding_rnas     += to_add->_av_size_non_coding_rnas;
  _nb_functional_genes         += to_add->_nb_functional_genes;
  _nb_non_functional_genes     += to_add->_nb_non_functional_genes;
  _av_size_functional_gene     += to_add->_av_size_functional_gene;
  _av_size_non_functional_gene += to_add->_av_size_non_functional_gene;

  _nb_mut    += to_add->_nb_mut;
  _nb_rear   += to_add->_nb_rear;
  _nb_switch += to_add->_nb_switch;
  _nb_indels += to_add->_nb_indels;
  _nb_dupl   += to_add->_nb_dupl;
  _nb_del    += to_add->_nb_del;
  _nb_trans  += to_add->_nb_trans;
  _nb_inv    += to_add->_nb_inv;
    
  _dupl_rate  += to_add->_dupl_rate;
  _del_rate   += to_add->_del_rate;
  _trans_rate += to_add->_trans_rate;
  _inv_rate   += to_add->_inv_rate;
  //~ printf("%f %f %f %f\n", to_add->_dupl_rate, to_add->_del_rate, to_add->_trans_rate, to_add->_inv_rate);
  _mean_align_score += to_add->_mean_align_score;
  
  _nb_bases_in_0_CDS                += to_add->_nb_bases_in_0_CDS;
  _nb_bases_in_0_functional_CDS     += to_add->_nb_bases_in_0_functional_CDS;
  _nb_bases_in_0_non_functional_CDS += to_add->_nb_bases_in_0_non_functional_CDS;
  _nb_bases_in_0_RNA                += to_add->_nb_bases_in_0_RNA;
  _nb_bases_in_0_coding_RNA         += to_add->_nb_bases_in_0_coding_RNA;
  _nb_bases_in_0_non_coding_RNA     += to_add->_nb_bases_in_0_non_coding_RNA;
    
  _nb_bases_non_essential                     += to_add->_nb_bases_non_essential;
  _nb_bases_non_essential_including_nf_genes  += to_add->_nb_bases_non_essential_including_nf_genes;
    
  #ifdef __REGUL
    _nb_influences                 += to_add->_nb_influences;
    _nb_enhancing_influences       += to_add->_nb_enhancing_influences;
    _nb_operating_influences       += to_add->_nb_operating_influences;
    _av_value_influences           += to_add->_av_value_influences;
    _av_value_enhancing_influences += to_add->_av_value_enhancing_influences;
    _av_value_operating_influences += to_add->_av_value_operating_influences;
  #endif
}

void StatRecord::substract_power(const StatRecord * means,
                                     const StatRecord * to_substract,
                                     double power)
{
  // NB : pop_size is a "global" value and must not be summed.
  _fitness                 += pow(means->_fitness - to_substract->_fitness, power);
  
  _metabolic_error         += pow(means->_metabolic_error - to_substract->_metabolic_error, power);
  _parent_metabolic_error  += pow(means->_parent_metabolic_error - to_substract->_parent_metabolic_error, power);
  _metabolic_fitness         += pow(means->_metabolic_fitness - to_substract->_metabolic_fitness, power);
  
  _secretion_error         += pow(means->_secretion_error - to_substract->_secretion_error, power);
  _parent_secretion_error  += pow(means->_parent_secretion_error - to_substract->_parent_secretion_error, power);
  
  _secretion_fitness       += pow(means->_secretion_fitness - to_substract->_secretion_fitness, power);
  _compound_amount         += pow(means->_compound_amount - to_substract->_compound_amount, power);
  
  _amount_of_dna               += pow(means->_amount_of_dna - to_substract->_amount_of_dna, power);
  _nb_coding_rnas              += pow(means->_nb_coding_rnas - to_substract->_nb_coding_rnas, power);
  _nb_non_coding_rnas          += pow(means->_nb_non_coding_rnas - to_substract->_nb_non_coding_rnas, power);
  _av_size_coding_rnas         += pow(means->_av_size_coding_rnas - to_substract->_av_size_coding_rnas, power);
  _av_size_non_coding_rnas     += pow(means->_av_size_non_coding_rnas - to_substract->_av_size_non_coding_rnas, power);
  _nb_functional_genes         += pow(means->_nb_functional_genes - to_substract->_nb_functional_genes, power);
  _nb_non_functional_genes     += pow(means->_nb_non_functional_genes - to_substract->_nb_non_functional_genes, power);
  _av_size_functional_gene     += pow(means->_av_size_functional_gene - to_substract->_av_size_functional_gene, power);
  _av_size_non_functional_gene += pow(means->_av_size_non_functional_gene - to_substract->_av_size_non_functional_gene, power);

  _nb_mut    += pow(means->_nb_mut - to_substract->_nb_mut, power);
  _nb_rear   += pow(means->_nb_rear - to_substract->_nb_rear, power);
  _nb_switch += pow(means->_nb_switch - to_substract->_nb_switch, power);
  _nb_indels += pow(means->_nb_indels - to_substract->_nb_indels, power);
  _nb_dupl   += pow(means->_nb_dupl - to_substract->_nb_dupl, power);
  _nb_del    += pow(means->_nb_del - to_substract->_nb_del, power);
  _nb_trans  += pow(means->_nb_trans - to_substract->_nb_trans, power);
  _nb_inv    += pow(means->_nb_inv - to_substract->_nb_inv, power);
    
  _dupl_rate  += pow(means->_dupl_rate - to_substract->_dupl_rate, power);
  _del_rate   += pow(means->_del_rate - to_substract->_del_rate, power);
  _trans_rate += pow(means->_trans_rate - to_substract->_trans_rate, power);
  _inv_rate   += pow(means->_inv_rate - to_substract->_inv_rate, power);
  
  _mean_align_score += pow(means->_mean_align_score - to_substract->_mean_align_score, power);
  
  _nb_bases_in_0_CDS                += pow(means->_nb_bases_in_0_CDS - to_substract->_nb_bases_in_0_CDS, power);
  _nb_bases_in_0_functional_CDS     += pow(means->_nb_bases_in_0_functional_CDS - to_substract->_nb_bases_in_0_functional_CDS, power);
  _nb_bases_in_0_non_functional_CDS += pow(means->_nb_bases_in_0_non_functional_CDS - to_substract->_nb_bases_in_0_non_functional_CDS, power);
  _nb_bases_in_0_RNA                += pow(means->_nb_bases_in_0_RNA - to_substract->_nb_bases_in_0_RNA, power);
  _nb_bases_in_0_coding_RNA         += pow(means->_nb_bases_in_0_coding_RNA - to_substract->_nb_bases_in_0_coding_RNA, power);
  _nb_bases_in_0_non_coding_RNA     += pow(means->_nb_bases_in_0_non_coding_RNA - to_substract->_nb_bases_in_0_non_coding_RNA, power);
    
  _nb_bases_non_essential                     += pow(means->_nb_bases_non_essential - to_substract->_nb_bases_non_essential, power);
  _nb_bases_non_essential_including_nf_genes  += pow(means->_nb_bases_non_essential_including_nf_genes - to_substract->_nb_bases_non_essential_including_nf_genes, power);
    
  #ifdef __REGUL
    _nb_influences                 += pow(means->_nb_influences - to_substract->_nb_influences, power);
    _nb_enhancing_influences       += pow(means->_nb_enhancing_influences - to_substract->_nb_enhancing_influences, power);
    _nb_operating_influences       += pow(means->_nb_operating_influences - to_substract->_nb_operating_influences, power);
    _av_value_influences           += pow(means->_av_value_influences - to_substract->_av_value_influences, power);
    _av_value_enhancing_influences += pow(means->_av_value_enhancing_influences - to_substract->_av_value_enhancing_influences, power);
    _av_value_operating_influences += pow(means->_av_value_operating_influences - to_substract->_av_value_operating_influences, power);
  #endif
}


// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
