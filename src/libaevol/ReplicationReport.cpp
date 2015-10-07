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
#include "ReplicationReport.h"
#include "DnaReplicationReport.h"
#include "Mutation.h"
#include "Individual.h"
#include "AeTime.h"
#include "Observable.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                         Class ReplicationReport                         #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ReplicationReport::ReplicationReport(Individual* indiv,
                                     const Individual* parent,
                                     Individual* donor /*= NULL*/)
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
  
  _parent_metabolic_error = parent->get_dist_to_target_by_feature(METABOLISM);
  _parent_secretion_error = parent->get_dist_to_target_by_feature(SECRETION);
  _parent_genome_size     = parent->get_total_genome_size();
  _mean_align_score       = 0.0;
  
  if (donor == NULL)
  {
    _donor_id               = -1;
    _donor_metabolic_error  = 0.0;
    _donor_secretion_error	= 0.0;
    _donor_genome_size      = 0;
  }
  else
  {
    _donor_id              = donor->get_id();
    _donor_metabolic_error = donor->get_dist_to_target_by_feature(METABOLISM);
    _donor_secretion_error = donor->get_dist_to_target_by_feature(SECRETION);
    _donor_genome_size     = donor->get_total_genome_size();
  }
}


// Creates an independent copy of the original report
ReplicationReport::ReplicationReport(const ReplicationReport& other) :
    _dna_replic_report(other._dna_replic_report)
{
  _parent_id  = other._parent_id;
  _donor_id   = other._donor_id;

  _id   = other._id;
  _rank = other._rank;

  _genome_size        = other._genome_size;
  _metabolic_error    = other._metabolic_error;
  _nb_genes_activ     = other._nb_genes_activ;
  _nb_genes_inhib     = other._nb_genes_inhib;
  _nb_non_fun_genes   = other._nb_non_fun_genes;
  _nb_coding_RNAs     = other._nb_coding_RNAs;
  _nb_non_coding_RNAs = other._nb_non_coding_RNAs;

  _parent_metabolic_error = other._parent_metabolic_error;
  _parent_secretion_error = other._parent_secretion_error;
  _donor_metabolic_error  = other._donor_metabolic_error;
  _donor_secretion_error  = other._donor_secretion_error;
  _parent_genome_size     = other._parent_genome_size;
  _donor_genome_size      = other._donor_genome_size;
  _mean_align_score       = other._mean_align_score;
}


ReplicationReport::ReplicationReport(gzFile tree_file, Individual* indiv)
{
  _indiv = indiv;
    
  gzread(tree_file, &_id,        sizeof(_id)        );
  gzread(tree_file, &_rank,      sizeof(_rank)      );
  gzread(tree_file, &_parent_id, sizeof(_parent_id) );
  gzread(tree_file, &_donor_id,  sizeof(_donor_id)  );
  
  gzread(tree_file, &_genome_size,         sizeof(_genome_size));
  gzread(tree_file, &_metabolic_error,     sizeof(_metabolic_error));
  gzread(tree_file, &_nb_genes_activ,      sizeof(_nb_genes_activ));
  gzread(tree_file, &_nb_genes_inhib,      sizeof(_nb_genes_inhib));
  gzread(tree_file, &_nb_non_fun_genes,    sizeof(_nb_non_fun_genes));
  gzread(tree_file, &_nb_coding_RNAs,      sizeof(_nb_coding_RNAs));
  gzread(tree_file, &_nb_non_coding_RNAs,  sizeof(_nb_non_coding_RNAs));

  _dna_replic_report.read_from_tree_file(tree_file);

  _dna_replic_report.compute_stats();
  
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

// =================================================================
//                            Public Methods
// =================================================================
/**
 * Set the individual corresponding to this replication report and
 * the characteristics of its parent
 *
 * This should be called as soon as a replication is started (just after calling
 * the offspring constructor and before doing the mutations)
 */
void ReplicationReport::init(Individual* offspring, Individual* parent)
{
  _indiv = offspring;

  _id = _indiv->get_id();
  _parent_id = parent->get_id();

  _genome_size        = 0;
  _metabolic_error    = 0.0;
  _nb_genes_activ     = 0;
  _nb_genes_inhib     = 0;
  _nb_non_fun_genes   = 0;
  _nb_coding_RNAs     = 0;
  _nb_non_coding_RNAs = 0;

  _parent_metabolic_error = parent->get_dist_to_target_by_feature(METABOLISM);
  _parent_secretion_error = parent->get_dist_to_target_by_feature(SECRETION);
  _parent_genome_size     = parent->get_total_genome_size();
  _mean_align_score       = 0.0;

  // Set ourselves an observer of indiv_'s MUTATION and END_REPLICATION
  _indiv->addObserver(this, MUTATION);
  _indiv->addObserver(this, END_REPLICATION);
}

/**
 * Method called at the end of the replication of an individual.
 * Actions such as finalize the calculation of average values can be done here.
 */
void ReplicationReport::signal_end_of_replication(Individual* indiv) {
  // TODO <david.parsons@inria.fr> tmp patch
  if (_indiv == NULL) _indiv = indiv;

  // Retrieve data from the individual
  _genome_size        = _indiv->get_total_genome_size();
  _metabolic_error    = _indiv->get_dist_to_target_by_feature(METABOLISM);
  _nb_genes_activ     = _indiv->get_nb_genes_activ();
  _nb_genes_inhib     = _indiv->get_nb_genes_inhib();
  _nb_non_fun_genes   = _indiv->get_nb_functional_genes();
  _nb_coding_RNAs     = _indiv->get_nb_coding_RNAs();
  _nb_non_coding_RNAs = _indiv->get_nb_non_coding_RNAs();
}

/**
 * Method called at the end of a generation.
 * Actions such as update the individuals' ranks can be done here.
 */
void ReplicationReport::signal_end_of_generation() {
  _rank = _indiv->get_rank();
}

void ReplicationReport::write_to_tree_file(gzFile tree_file) const
{
  // Store individual identifiers and rank
  gzwrite(tree_file, &_id,         sizeof(_id)        );
  assert(_rank != -1);
  gzwrite(tree_file, &_rank,       sizeof(_rank)      );
  gzwrite(tree_file, &_parent_id,  sizeof(_parent_id) );
  gzwrite(tree_file, &_donor_id,   sizeof(_donor_id)  );
  
  gzwrite(tree_file, &_genome_size,         sizeof(_genome_size));
  gzwrite(tree_file, &_metabolic_error,     sizeof(_metabolic_error));
  gzwrite(tree_file, &_nb_genes_activ,      sizeof(_nb_genes_activ));
  gzwrite(tree_file, &_nb_genes_inhib,      sizeof(_nb_genes_inhib));
  gzwrite(tree_file, &_nb_non_fun_genes,    sizeof(_nb_non_fun_genes));
  gzwrite(tree_file, &_nb_coding_RNAs,      sizeof(_nb_coding_RNAs));
  gzwrite(tree_file, &_nb_non_coding_RNAs,  sizeof(_nb_non_coding_RNAs));

  _dna_replic_report.write_to_tree_file(tree_file);
}


// =================================================================
//                           Protected Methods
// =================================================================





// =================================================================
//                          Non inline accessors
// =================================================================
void ReplicationReport::update(Observable& o, ObservableEvent e, void* arg) {
  switch (e) {
    case END_REPLICATION :
      signal_end_of_replication(dynamic_cast<Individual*>(&o));
      break;
    case MUTATION :
      _dna_replic_report.add_mut(reinterpret_cast<Mutation*>(arg));
      break;
    default :
      Utils::ExitWithDevMsg("Event not handled", __FILE__, __LINE__);
  }
}
} // namespace aevol
