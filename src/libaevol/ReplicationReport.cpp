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
  indiv_ = indiv;

  id_   = indiv->id();
  rank_ = indiv->rank();

  parent_id_ = parent->id();
  // donor_id_ is set further down

  genome_size_        = 0;
  metabolic_error_    = 0.0;
  nb_genes_activ_     = 0;
  nb_genes_inhib_     = 0;
  nb_non_fun_genes_   = 0;
  nb_coding_RNAs_     = 0;
  nb_non_coding_RNAs_ = 0;

  parent_metabolic_error_ = parent->dist_to_target_by_feature(METABOLISM);
  parent_secretion_error_ = parent->dist_to_target_by_feature(SECRETION);
  parent_genome_size_     = parent->total_genome_size();
  mean_align_score_       = 0.0;

  if (donor == NULL)
  {
    donor_id_               = -1;
    donor_metabolic_error_  = 0.0;
    donor_secretion_error_	= 0.0;
    donor_genome_size_      = 0;
  }
  else
  {
    donor_id_              = donor->id();
    donor_metabolic_error_ = donor->dist_to_target_by_feature(METABOLISM);
    donor_secretion_error_ = donor->dist_to_target_by_feature(SECRETION);
    donor_genome_size_     = donor->total_genome_size();
  }
}


// Creates an independent copy of the original report
ReplicationReport::ReplicationReport(const ReplicationReport& other) :
    dna_replic_report_(other.dna_replic_report_)
{
  parent_id_  = other.parent_id_;
  donor_id_   = other.donor_id_;

  id_   = other.id_;
  rank_ = other.rank_;

  genome_size_        = other.genome_size_;
  metabolic_error_    = other.metabolic_error_;
  nb_genes_activ_     = other.nb_genes_activ_;
  nb_genes_inhib_     = other.nb_genes_inhib_;
  nb_non_fun_genes_   = other.nb_non_fun_genes_;
  nb_coding_RNAs_     = other.nb_coding_RNAs_;
  nb_non_coding_RNAs_ = other.nb_non_coding_RNAs_;

  parent_metabolic_error_ = other.parent_metabolic_error_;
  parent_secretion_error_ = other.parent_secretion_error_;
  donor_metabolic_error_  = other.donor_metabolic_error_;
  donor_secretion_error_  = other.donor_secretion_error_;
  parent_genome_size_     = other.parent_genome_size_;
  donor_genome_size_      = other.donor_genome_size_;
  mean_align_score_       = other.mean_align_score_;
}


ReplicationReport::ReplicationReport(gzFile tree_file, Individual* indiv)
{
  indiv_ = indiv;

  gzread(tree_file, &id_,        sizeof(id_));
  gzread(tree_file, &rank_,      sizeof(rank_));
  gzread(tree_file, &parent_id_, sizeof(parent_id_));
  gzread(tree_file, &donor_id_,  sizeof(donor_id_));

  gzread(tree_file, &genome_size_,         sizeof(genome_size_));
  gzread(tree_file, &metabolic_error_,     sizeof(metabolic_error_));
  gzread(tree_file, &nb_genes_activ_,      sizeof(nb_genes_activ_));
  gzread(tree_file, &nb_genes_inhib_,      sizeof(nb_genes_inhib_));
  gzread(tree_file, &nb_non_fun_genes_,    sizeof(nb_non_fun_genes_));
  gzread(tree_file, &nb_coding_RNAs_,      sizeof(nb_coding_RNAs_));
  gzread(tree_file, &nb_non_coding_RNAs_,  sizeof(nb_non_coding_RNAs_));

  dna_replic_report_.read_from_tree_file(tree_file);

  dna_replic_report_.compute_stats();

  parent_metabolic_error_ = -1;
  parent_secretion_error_ = -1;
  donor_metabolic_error_  = -1;
  parent_genome_size_     = -1;
  donor_genome_size_      = -1;
  mean_align_score_       = 0.0;
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
  indiv_ = offspring;

  id_ = indiv_->id();
  parent_id_ = parent->id();

  genome_size_        = 0;
  metabolic_error_    = 0.0;
  nb_genes_activ_     = 0;
  nb_genes_inhib_     = 0;
  nb_non_fun_genes_   = 0;
  nb_coding_RNAs_     = 0;
  nb_non_coding_RNAs_ = 0;

  parent_metabolic_error_ = parent->dist_to_target_by_feature(METABOLISM);
  parent_secretion_error_ = parent->dist_to_target_by_feature(SECRETION);
  parent_genome_size_     = parent->total_genome_size();
  mean_align_score_       = 0.0;

  // Set ourselves an observer of indiv_'s MUTATION and END_REPLICATION
  indiv_->addObserver(this, MUTATION);
  indiv_->addObserver(this, END_REPLICATION);
}

/**
 * Method called at the end of the replication of an individual.
 * Actions such as finalize the calculation of average values can be done here.
 */
void ReplicationReport::signal_end_of_replication(Individual* indiv) {
  // TODO <david.parsons@inria.fr> tmp patch
  if (indiv_ == NULL) indiv_ = indiv;

  // Retrieve data from the individual
  genome_size_        = indiv_->total_genome_size();
  metabolic_error_    = indiv_->dist_to_target_by_feature(METABOLISM);
  nb_genes_activ_     = indiv_->nb_genes_activ();
  nb_genes_inhib_     = indiv_->nb_genes_inhib();
  nb_non_fun_genes_   = indiv_->nb_functional_genes();
  nb_coding_RNAs_     = indiv_->nb_coding_RNAs();
  nb_non_coding_RNAs_ = indiv_->nb_non_coding_RNAs();
}

/**
 * Method called at the end of a generation.
 * Actions such as update the individuals' ranks can be done here.
 */
void ReplicationReport::signal_end_of_generation() {
  rank_ = indiv_->rank();
}

void ReplicationReport::write_to_tree_file(gzFile tree_file) const
{
  // Store individual identifiers and rank
  gzwrite(tree_file, &id_,         sizeof(id_));
  assert(rank_ != -1);
  gzwrite(tree_file, &rank_,       sizeof(rank_));
  gzwrite(tree_file, &parent_id_,  sizeof(parent_id_));
  gzwrite(tree_file, &donor_id_,   sizeof(donor_id_));

  gzwrite(tree_file, &genome_size_,         sizeof(genome_size_));
  gzwrite(tree_file, &metabolic_error_,     sizeof(metabolic_error_));
  gzwrite(tree_file, &nb_genes_activ_,      sizeof(nb_genes_activ_));
  gzwrite(tree_file, &nb_genes_inhib_,      sizeof(nb_genes_inhib_));
  gzwrite(tree_file, &nb_non_fun_genes_,    sizeof(nb_non_fun_genes_));
  gzwrite(tree_file, &nb_coding_RNAs_,      sizeof(nb_coding_RNAs_));
  gzwrite(tree_file, &nb_non_coding_RNAs_,  sizeof(nb_non_coding_RNAs_));

  dna_replic_report_.write_to_tree_file(tree_file);
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
      dna_replic_report_.add_mut(reinterpret_cast<Mutation*>(arg));
      break;
    default :
      Utils::ExitWithDevMsg("Event not handled", __FILE__, __LINE__);
  }
}
} // namespace aevol
