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
#include <algorithm>
#include "Tree.h"

#include "macros.h"
#include "ExpManager.h"
#include "ExpSetup.h"
#include "Individual.h"
#include "Utils.h"
#include "SIMD_Individual.h"

namespace aevol {


//##############################################################################
//                                                                             #
//                                Class Tree                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
const int32_t Tree::NO_PARENT = -1;

// =================================================================
//                             Constructors
// =================================================================
Tree::Tree(ExpManager* exp_m, int64_t tree_step) {
  exp_m_ = exp_m;
  tree_step_ = tree_step;

  replics_ = new ReplicationReport** [tree_step_];

  for (int32_t time = 0 ; time < tree_step ; time++) {
    replics_[time] = new ReplicationReport* [exp_m->nb_indivs()];
    for (int32_t num_indiv = 0 ;
         num_indiv < exp_m->nb_indivs() ;
         num_indiv++) {
      replics_[time][num_indiv] = nullptr;
    }
  }
}


/**
 *
 */
Tree::Tree(ExpManager* exp_m, char* tree_file_name) {
  exp_m_ = exp_m;
  tree_step_ = exp_m_->tree_step();

  gzFile tree_file = gzopen(tree_file_name, "r");
  if (tree_file == Z_NULL) {
    printf("ERROR : Could not read tree file %s\n", tree_file_name);
    exit(EXIT_FAILURE);
  }

        printf("---- ] Allocate tree %ld\n",tree_step_);
  replics_ = new ReplicationReport** [tree_step_];
  //for (int64_t t = AeTime::time()-tree_step_+1 ; t <= AeTime::time() ; t++) {
  for (int64_t t = 0 ; t < tree_step_ ; t++) {
    replics_[t] = new ReplicationReport* [exp_m_->nb_indivs()];

    for (int32_t indiv_i = 0 ;
         indiv_i < exp_m_->nb_indivs() ;
         indiv_i++) {
      // Retrieve a replication report
      ReplicationReport* replic_report = new ReplicationReport(tree_file,
                                                               nullptr);

      // Put it at its rightful position
        replics_[t][replic_report->id()] = replic_report;
    }
  }
  gzclose(tree_file);
}




// =================================================================
//                             Destructors
// =================================================================
Tree::~Tree() {
  if (replics_ != NULL)  {
    for (int32_t i = 0 ; i < tree_step_ ; i++)
      if (replics_[i] != NULL) {
        for (int32_t j = 0 ; j < exp_m_->nb_indivs() ; j++)
          delete replics_[i][j];
        delete [] replics_[i];
      }
    delete [] replics_;
  }
}

// =================================================================
//                            Public Methods
// =================================================================
    ReplicationReport** Tree::reports(int64_t t) const {
        return replics_[Utils::mod(t - 1, tree_step_)];
    }

    ReplicationReport* Tree::report_by_index(int64_t t, int32_t index) const {
        return replics_[Utils::mod(t - 1, tree_step_)][index];
    }

    ReplicationReport* Tree::report_by_rank(int64_t t, int32_t rank) const {
        int32_t nb_indivs = exp_m_->nb_indivs();
        assert(rank <= nb_indivs);

        for (int32_t i = 0 ; i < nb_indivs ; i++) {
            if (replics_[Utils::mod(t - 1, tree_step_)][i]->rank() == rank) {
                return replics_[Utils::mod(t - 1, tree_step_)][i];
            }
        }

        fprintf(stderr,
                "ERROR: Couldn't find indiv with rank %" PRId32 " in file %s:%d\n",
                rank, __FILE__, __LINE__);
        return NULL;
    }

    void Tree::signal_end_of_generation() {
        auto cur_reports = reports(AeTime::time());
        for (int32_t i = 0; i < exp_m_->nb_indivs(); i++) {
            if (cur_reports[i] == nullptr){
                printf("error!!!\n");exit(-1);
            }

            cur_reports[i]->signal_end_of_generation();
        }
    }

    void Tree::write_to_tree_file(char* tree_file_name) {
        gzFile tree_file = gzopen( tree_file_name, "w" );
        printf("Tree Write %d : \n",time());//,t,indiv_i);
        // Write the tree in the backup
        for (int64_t t = 0 ; t < tree_step_ ; t++)
            for (int32_t indiv_i = 0 ; indiv_i < exp_m_->nb_indivs() ; indiv_i++) {

                assert(replics_[t][indiv_i] != NULL);
                replics_[t][indiv_i]->write_to_tree_file(tree_file);
            }

        gzclose(tree_file);


        // Reinitialize the tree
        for (int32_t time = 0 ; time < tree_step_ ; time++) {
            for (int32_t num_indiv = 0 ;
                 num_indiv < exp_m_->nb_indivs() ;
                 num_indiv++) {
                delete replics_[time][num_indiv];
                replics_[time][num_indiv] = nullptr;
            }
        }

    }

void Tree::update_new_indiv(NewIndivEvent* evt) {

    replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
                                    evt->simd_child->exp_m_->grid_height()
                                    + evt->y] = new ReplicationReport();
    replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
                                                         evt->simd_child->exp_m_->grid_height()
                                                         + evt->y]->
            init(this, evt->simd_child, evt->simd_parent, evt->indiv_id_, evt->parent_id_);
}

void Tree::update_end_replication(EndReplicationEvent* evt) {
    replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
                                                         evt->simd_child->exp_m_->grid_height()
                                                         + evt->y]->signal_end_of_replication(evt->simd_child);
}


void Tree::update_end_generation() {
    signal_end_of_generation();
}

void Tree::update(Observable& o, ObservableEvent e, void* arg) {

  switch (e) {
    case NEW_INDIV : {
      // Initialize the replication report corresponding to the new individual

      auto ievent = reinterpret_cast<NewIndivEvent*>(arg);
       if (SIMD_Individual::standalone_simd) {
           replics_[Utils::mod(AeTime::time() - 1, tree_step_)][ievent->x *
                                              ievent->simd_child->exp_m_->grid_height()
                                              + ievent->y]->
                      init(this, ievent->simd_child, ievent->simd_parent, ievent->indiv_id_, ievent->parent_id_);
          } else {
           replics_[Utils::mod(AeTime::time() - 1, tree_step_)][ievent->x *
                                              ievent->child->exp_m()->grid_height()
                                              + ievent->y]->
                      init(this, ievent->child, ievent->parent, ievent->indiv_id_, ievent->parent_id_);
          }

      break;
    }
    case END_GENERATION : {
      signal_end_of_generation();
      break;
    }
    case END_REPLICATION : {
      auto ievent = reinterpret_cast<EndReplicationEvent*>(arg);

            if (SIMD_Individual::standalone_simd) {
                //printf("EoR %d : %p -- %p\n",ievent->simd_child->indiv_id,ievent->simd_child, replics_[AeTime::time()][ievent->simd_child->indiv_id]);
                replics_[Utils::mod(AeTime::time() - 1, tree_step_)][ievent->x *
                                                ievent->simd_child->exp_m_->grid_height()
                                                + ievent->y]->signal_end_of_replication(
                        ievent->simd_child);
            } else {
                replics_[Utils::mod(AeTime::time() - 1, tree_step_)][ievent->x *
                                                ievent->child->exp_m()->grid_height()
                                                + ievent->y]->signal_end_of_replication(
                        ievent->child);
            }

      break;
    }
    default : {
      Utils::ExitWithDevMsg("Event not handled", __FILE__, __LINE__);
    }
  }
}



// =================================================================
//                  Non-inline accessors' definition
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
