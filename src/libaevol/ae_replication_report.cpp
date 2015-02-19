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



// =================================================================
//                            Project Files
// =================================================================
#include "ae_replication_report.h"
#include "ae_dna_replic_report.h"
#include "ae_mutation.h"
#include "ae_individual.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                         Class ae_replication_report                         #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_replication_report::ae_replication_report( ae_individual * indiv, ae_individual * parent, ae_individual * donor /* = NULL */ )
{
  _indiv = indiv;
  
  _id   = indiv->get_id();
  _rank = indiv->get_rank();
  
  _parent_id = parent->get_id();
  // _donor_id is set further down
    
  _genome_size        = 0;
  _metabolic_error    = 0.0;
  _nb_genes_activ     = 0;
  _nb_genes_inhib     = 0;
  _nb_non_fun_genes   = 0;
  _nb_coding_RNAs     = 0;
  _nb_non_coding_RNAs = 0;
  
  _parent_metabolic_error = parent->get_dist_to_target_by_feature( METABOLISM );
  _parent_secretion_error = parent->get_dist_to_target_by_feature( SECRETION );
  _parent_genome_size     = parent->get_total_genome_size();
  _mean_align_score       = 0.0;
  
  if ( donor == NULL )
  {
    _donor_id               = -1;
    _donor_metabolic_error  = 0.0;
    _donor_secretion_error	= 0.0;
    _donor_genome_size      = 0;
  }
  else
  {
    _donor_id              = donor->get_id();
    _donor_metabolic_error = donor->get_dist_to_target_by_feature( METABOLISM );
    _donor_secretion_error = donor->get_dist_to_target_by_feature( SECRETION );
    _donor_genome_size     = donor->get_total_genome_size();
  }
    
}


// Creates an independent copy of the original report
ae_replication_report::ae_replication_report( const ae_replication_report &model )
{
  _parent_id  = model._parent_id;
  _donor_id   = model._donor_id;
  
  _id   = model._id;
  _rank = model._rank;
    
  _genome_size        = model._genome_size;
  _metabolic_error    = model._metabolic_error;
  _nb_genes_activ     = model._nb_genes_activ;
  _nb_genes_inhib     = model._nb_genes_inhib;
  _nb_non_fun_genes   = model._nb_non_fun_genes;
  _nb_coding_RNAs     = model._nb_coding_RNAs;
  _nb_non_coding_RNAs = model._nb_non_coding_RNAs;

  _parent_metabolic_error = model._parent_metabolic_error;
  _parent_secretion_error = model._parent_secretion_error;
  _donor_metabolic_error  = model._donor_metabolic_error;
  _donor_secretion_error  = model._donor_secretion_error;
  _parent_genome_size     = model._parent_genome_size;
  _donor_genome_size      = model._donor_genome_size;
  _mean_align_score       = model._mean_align_score;

  for (const auto& rep: model._dna_replic_reports)
    _dna_replic_reports.push_back(new ae_dna_replic_report(*rep));
}


ae_replication_report::ae_replication_report( gzFile tree_file, ae_individual * indiv )
{
  _indiv = indiv;
    
  gzread( tree_file, &_id,        sizeof(_id)         );
  gzread( tree_file, &_rank,      sizeof(_rank)       );
  gzread( tree_file, &_parent_id, sizeof(_parent_id)  );
  gzread( tree_file, &_donor_id,  sizeof(_donor_id)   );
  
  gzread( tree_file, &_genome_size,         sizeof(_genome_size) );
  gzread( tree_file, &_metabolic_error,     sizeof(_metabolic_error) );
  gzread( tree_file, &_nb_genes_activ,      sizeof(_nb_genes_activ) );
  gzread( tree_file, &_nb_genes_inhib,      sizeof(_nb_genes_inhib) );
  gzread( tree_file, &_nb_non_fun_genes,    sizeof(_nb_non_fun_genes) );
  gzread( tree_file, &_nb_coding_RNAs,      sizeof(_nb_coding_RNAs) );
  gzread( tree_file, &_nb_non_coding_RNAs,  sizeof(_nb_non_coding_RNAs) );
  
  int32_t nb_dna_replic_reports;
  gzread( tree_file, &nb_dna_replic_reports, sizeof(nb_dna_replic_reports) );
  
  int32_t mydnareport, myevent;
  int32_t nb_rears, nb_muts, nb_HT;
  ae_dna_replic_report * dnareport = NULL;
  ae_mutation * event = NULL;

  for ( mydnareport = 0 ; mydnareport < nb_dna_replic_reports ; mydnareport++ )
  {
    dnareport = new ae_dna_replic_report();
    
    gzread( tree_file, &nb_HT, sizeof(nb_HT) );
    for ( myevent  = 0 ; myevent < nb_HT ; myevent++ )
    {
      event = new ae_mutation( tree_file );
      dnareport->add_HT( event );
    }
    
    gzread( tree_file, &nb_rears, sizeof(nb_rears) );
    for ( myevent  = 0 ; myevent < nb_rears ; myevent++ )
    {
      event = new ae_mutation( tree_file );
      dnareport->get_rearrangements()->add( event );
    }
    
    gzread( tree_file, &nb_muts, sizeof(nb_muts) );
    for(myevent  = 0 ; myevent < nb_muts ; myevent++ )
    {
      event = new ae_mutation( tree_file );
      dnareport->get_mutations()->add( event );
    }

    dnareport->compute_stats();
    _dna_replic_reports.push_back(dnareport);
  }
  
  _parent_metabolic_error = -1;
  _parent_secretion_error = -1;
  _donor_metabolic_error  = -1;
  _parent_genome_size     = -1;
  _donor_genome_size      = -1;
  _mean_align_score       = 0.0;
}


// =================================================================
//                             Destructors
// =================================================================
ae_replication_report::~ae_replication_report( void )
{
  for (const auto& rep: _dna_replic_reports)
    delete rep;
}

// =================================================================
//                            Public Methods
// =================================================================
/**
 * Method called at the end of the replication. Actions such as finalize the calculation of average values can be done here.
 */
void ae_replication_report::signal_end_of_replication( void )
{
  // Retreive data from the individual
  _genome_size        = _indiv->get_total_genome_size();
  _metabolic_error    = _indiv->get_dist_to_target_by_feature( METABOLISM );
  _nb_genes_activ     = _indiv->get_nb_genes_activ();
  _nb_genes_inhib     = _indiv->get_nb_genes_inhib();
  _nb_non_fun_genes   = _indiv->get_nb_functional_genes();
  _nb_coding_RNAs     = _indiv->get_nb_coding_RNAs();
  _nb_non_coding_RNAs = _indiv->get_nb_non_coding_RNAs();
  
  
  // Compute the mean alignment score
  int32_t nb_align = 0;
  ae_mutation_type        mut_type;
  int32_t*  int32_trash = new int32_t;
  bool*     bool_trash  = new bool;
  int16_t*  align_scores = new int16_t[2];
  for (const auto& dna_rep: _dna_replic_reports) {
    nb_align += dna_rep->get_nb_duplications();
    nb_align += dna_rep->get_nb_deletions();
    nb_align += 2 * dna_rep->get_nb_translocations();
    nb_align += dna_rep->get_nb_inversions();
    
    for (const auto& rear: dna_rep->get_rearrangements_std()) {
      mut_type = rear->get_mut_type();
      
      switch( mut_type )
      {
        case DUPL:
          rear->get_infos_duplication( int32_trash, int32_trash, int32_trash, &align_scores[0] );
          _mean_align_score += align_scores[0];
          break;
        case DEL:
          rear->get_infos_deletion( int32_trash, int32_trash, &align_scores[0] );
          _mean_align_score += align_scores[0];
          break;
        case TRANS:
          rear->get_infos_translocation( int32_trash, int32_trash, int32_trash, int32_trash, bool_trash, &align_scores[0], &align_scores[1] );
          _mean_align_score += align_scores[0] + align_scores[1];
          break;
        case INV:
          rear->get_infos_inversion( int32_trash, int32_trash, &align_scores[0] );
          _mean_align_score += align_scores[0];
          break;
        default:
          fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", mut_type, __FILE__, __LINE__ );
          exit( EXIT_FAILURE );
      }
    }
  }
  
  if ( nb_align != 0 )
  {
    _mean_align_score /= nb_align;
  }
  else
  {
    assert( _mean_align_score == 0.0 );
  }
  
  delete int32_trash;
  delete bool_trash;
  delete [] align_scores;
}

void ae_replication_report::write_to_tree_file( gzFile tree_file ) const
{
  // Store individual identifiers and rank
  gzwrite( tree_file, &_id,         sizeof(_id)         );
  gzwrite( tree_file, &_rank,       sizeof(_rank)       );
  gzwrite( tree_file, &_parent_id,  sizeof(_parent_id)  );
  gzwrite( tree_file, &_donor_id,   sizeof(_donor_id)   );
  
  gzwrite( tree_file, &_genome_size,         sizeof(_genome_size) );
  gzwrite( tree_file, &_metabolic_error,     sizeof(_metabolic_error) );
  gzwrite( tree_file, &_nb_genes_activ,      sizeof(_nb_genes_activ) );
  gzwrite( tree_file, &_nb_genes_inhib,      sizeof(_nb_genes_inhib) );
  gzwrite( tree_file, &_nb_non_fun_genes,    sizeof(_nb_non_fun_genes) );
  gzwrite( tree_file, &_nb_coding_RNAs,      sizeof(_nb_coding_RNAs) );
  gzwrite( tree_file, &_nb_non_coding_RNAs,  sizeof(_nb_non_coding_RNAs) );  
  
  // For each genetic unit, write the mutations and rearrangements undergone during replication
  int32_t nb_dna_replic_reports = _dna_replic_reports.size();
  gzwrite( tree_file, &nb_dna_replic_reports, sizeof(nb_dna_replic_reports) );
  //~ printf( "  nb_dna_replic_reports : %"PRId32"\n", nb_dna_replic_reports );

  for (const auto& report: _dna_replic_reports) {
    // Store HT
    int32_t nb_HT = report->get_nb_HT();
    gzwrite( tree_file, &nb_HT, sizeof(nb_HT) );
    for (const auto& HT: report->get_HT_std())
      HT->save(tree_file);

    // Store rearrangements
    int32_t nb_rears = report->get_nb_rearrangements();
    gzwrite( tree_file, &nb_rears, sizeof(nb_rears) );
    //~ printf( "  nb_rears : %"PRId32"\n", nb_rears );

    for (const auto& rear: report->get_rearrangements_std())
      rear->save(tree_file);

    // Store mutations
    int32_t nb_muts = report->get_nb_small_mutations();
    gzwrite( tree_file, &nb_muts, sizeof(nb_muts) );
    //~ printf( "  nb_muts : %"PRId32"\n", nb_muts );

    for (const auto& mutation: report->get_mutations_std())
      mutation->save( tree_file );
  }
}


// =================================================================
//                           Protected Methods
// =================================================================





// =================================================================
//                          Non inline accessors
// =================================================================

} // namespace aevol
