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
//                              Includes
// =================================================================
#include "OutputManager.h"

#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <fstream>
#include <string>

#include "ExpManager.h"
#include "AeTime.h"

using std::string;
using std::endl;

namespace aevol {



//##############################################################################
//                                                                             #
//                           Class OutputManager                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
OutputManager::OutputManager(ExpManager * exp_m) {
  exp_m_  = exp_m;
  stats_  = nullptr;
  tree_   = nullptr;
  dump_   = nullptr;
  compute_phen_contrib_by_GU_ = false;
  record_tree_ = false;
  make_dumps_ = false;
  dump_step_ = 0;
  logs_ = new Logging();
}

// =================================================================
//                             Destructors
// =================================================================
OutputManager::~OutputManager() {
  delete stats_;
  delete tree_;
  delete dump_;
  delete logs_;
}

// =================================================================
//                            Public Methods
// =================================================================
void OutputManager::InitStats() {
  stats_ = new Stats(exp_m_);
}

void OutputManager::WriteSetupFile(gzFile setup_file) const {
  // Write the backup steps
  gzwrite(setup_file, &backup_step_,      sizeof(backup_step_));
  gzwrite(setup_file, &big_backup_step_,  sizeof(big_backup_step_));

  // Stats
  gzwrite(setup_file, &compute_phen_contrib_by_GU_,  sizeof(compute_phen_contrib_by_GU_));

  // Tree
  int8_t record_tree = record_tree_;
  gzwrite(setup_file, &record_tree, sizeof(record_tree));
  if (record_tree_)
  {
    auto tmp_tree_step = tree_->tree_step();
    gzwrite(setup_file, &tmp_tree_step, sizeof(tmp_tree_step));
  }

  // Dumps
  int8_t make_dumps = make_dumps_;
  gzwrite(setup_file, &make_dumps,  sizeof(make_dumps));
  gzwrite(setup_file, &dump_step_,  sizeof(dump_step_));

  // Logs
  int8_t logs = logs_->logs();
  gzwrite(setup_file, &logs,  sizeof(logs));
}

/**
 * This is a temporary patch for experiment propagation, it shall become
 * obsolete or need to be adapted when in/out dirs are managed properly
 */
void OutputManager::PropagateStats(const std::string& outdir,
                                   int64_t propagated_timestep) const {
  Stats stats("stat");
  stats.Propagate(outdir, propagated_timestep);
}

void OutputManager::load(gzFile setup_file, bool verbose, bool to_be_run) {
  // Write the backup steps
  gzread(setup_file, &backup_step_,      sizeof(backup_step_));
  gzread(setup_file, &big_backup_step_,  sizeof(big_backup_step_));

  // Stats
  if (to_be_run) {
    delete stats_;
    stats_ = new Stats(exp_m_, AeTime::time());
  }
  gzread(setup_file, &compute_phen_contrib_by_GU_,  sizeof(compute_phen_contrib_by_GU_));

  // Tree
  int8_t record_tree;
  gzread(setup_file, &record_tree, sizeof(record_tree));
  record_tree_ = record_tree;
  if (record_tree_) {
    int32_t tmp_tree_step;
    gzread(setup_file, &tmp_tree_step, sizeof(tmp_tree_step));

    tree_ = new Tree(exp_m_, tmp_tree_step);
  }

  // Dumps
  int8_t make_dumps;
  gzread(setup_file, &make_dumps, sizeof(make_dumps));
  make_dumps_ = make_dumps;
  gzread(setup_file, &dump_step_,  sizeof(dump_step_));
  if(make_dumps_) {
    dump_ = new Dump(exp_m_);
  }

  // Logs
  int8_t logs;
  gzread(setup_file, &logs, sizeof(logs));
  if (to_be_run) {
    logs_->load(logs, AeTime::time());
  }
}

void OutputManager::write_current_generation_outputs() const {
  // Write stats
  stats_->write_current_generation_statistics();

  // Manage tree
  if (record_tree_ &&
      AeTime::time() > 0 &&
      (AeTime::time() % tree_->tree_step() == 0)) {
    write_tree();
  }

  // Write backup
  if (AeTime::time() % backup_step_ == 0) {
    stats_->flush();
    exp_m_->WriteDynamicFiles();

    WriteLastGenerFile();
  }

  // Write dumps
  if (make_dumps_) {
    if(AeTime::time() % dump_step_ == 0) {
      dump_->write_current_generation_dump();
    }
  }
}

// TODO <david.parsons@inria.fr> we need an output_dir attribute in this class !
void OutputManager::WriteLastGenerFile(const string& output_dir) const {
  std::ofstream last_gener_file(output_dir + "/" + LAST_GENER_FNAME,
                                std::ofstream::out);
  if (last_gener_file.fail()) {
    Utils::ExitWithUsrMsg(string("could not open file ") + LAST_GENER_FNAME);
  }
  else {
    last_gener_file << AeTime::time() << endl;
    last_gener_file.close();
  }
}

// TODO <david.parsons@inria.fr> we need an input_dir attribute in this class !
int64_t OutputManager::last_gener() {
  int64_t time;
  FILE* lg_file = fopen(LAST_GENER_FNAME, "r");
  if (lg_file != NULL) {
    if (fscanf(lg_file, "%" PRId64 "\n", &time) == EOF) {
      Utils::ExitWithDevMsg("failed to read last generation", __FILE__,
                            __LINE__);
    }
    fclose(lg_file);
  }
  return time;
}

// =================================================================
//                           Protected Methods
// =================================================================
void OutputManager::write_tree() const
{
  // Create the tree directory if it doesn't exist
  int status;
  status = mkdir(TREE_DIR, 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, "Impossible to create the directory %s", TREE_DIR);
  }

  char tree_file_name[50];

  sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae", AeTime::time());

  
  gzFile tree_file = gzopen( tree_file_name, "w" );
  // Write phylogenetic data (tree)
  tree_->write_to_tree_file(tree_file);

  gzclose(tree_file);
}

// =================================================================
//                          Non inline accessors
// =================================================================
} // namespace aevol
