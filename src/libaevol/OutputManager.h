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


#ifndef AEVOL_OUPUT_MANAGER_H_
#define AEVOL_OUPUT_MANAGER_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_enums.h"
#include "Stats.h"
#include "Tree.h"
#include "Dump.h"
#include "Logging.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;






class OutputManager {
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  OutputManager() = delete;
  OutputManager(const OutputManager&) = delete;
  OutputManager(ExpManager* exp_m);

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~OutputManager();

  // =================================================================
  //                        Accessors: getters
  // =================================================================

  // Backup
  inline int64_t	backup_step() const;
  inline int64_t	big_backup_step() const;

  // Tree
  inline bool record_tree() const;
  inline int64_t tree_step() const;
  inline Tree* tree() const;

  // Logs
  inline FILE* log(LogType log_type) const;
  inline bool is_logged(LogType log_type) const;

  // Stats
  inline bool compute_phen_contrib_by_GU() const;

  // =================================================================
  //                        Accessors: setters
  // =================================================================
  inline void set_backup_step(int64_t backup_step);
  inline void set_big_backup_step(int64_t big_backup_step);
  inline void init_tree(ExpManager* exp_m, int64_t tree_step_);
  inline void set_dump_step(int64_t dump_step);
  inline void set_compute_phen_contrib_by_GU(bool compute_phen_contrib_by_GU);
  inline void set_logs (int8_t logs);

  // =================================================================
  //                              Operators
  // =================================================================

  // =================================================================
  //                            Public Methods
  // =================================================================
  void InitStats();
  void WriteSetupFile(gzFile setup_file) const;
  void WriteLastGenerFile(const std::string& input_dir = ".") const;
  void PropagateStats(const std::string& outdir,
                      int64_t propagated_timestep) const;
  void load(gzFile file, bool verbose, bool to_be_run);
  void write_current_generation_outputs() const;
  inline void flush();

  static int64_t last_gener();

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  void write_tree() const;

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  ExpManager* exp_m_;

  // Backups
  int64_t backup_step_;
  int64_t big_backup_step_;

  // Stats
  Stats* stats_;
  bool compute_phen_contrib_by_GU_;

  // Tree
  bool record_tree_;
  Tree* tree_;

  // Dumps
  bool make_dumps_;
  int64_t dump_step_;
  Dump* dump_;

  // Logs
  Logging* logs_;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// Backup
inline int64_t OutputManager::backup_step() const {
  return backup_step_;
}

inline int64_t OutputManager::big_backup_step() const {
  return big_backup_step_;
}

// Tree
inline bool OutputManager::record_tree() const {
  return record_tree_;
}

inline int64_t OutputManager::tree_step() const {
  return tree_->tree_step();
}

inline Tree *OutputManager::tree() const {
  return tree_;
}

// Logs
inline FILE* OutputManager::log(LogType log_type) const {
  return logs_->log(log_type);
}

inline bool  OutputManager::is_logged(LogType log_type) const {
  return logs_->is_logged(log_type);
}

// Stats
inline bool OutputManager::compute_phen_contrib_by_GU() const {
  return compute_phen_contrib_by_GU_;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
void OutputManager::set_backup_step(int64_t backup_step) {
  backup_step_ = backup_step;
}

void OutputManager::set_big_backup_step(int64_t big_backup_step) {
  big_backup_step_ = big_backup_step;
}

void OutputManager::init_tree(ExpManager* exp_m, int64_t tree_step_) {
  record_tree_ = true;
  tree_ = new Tree(exp_m, tree_step_);
}

void OutputManager::set_dump_step(int64_t dump_step) {
  make_dumps_ = true;
  dump_step_  = dump_step;
}

void OutputManager::set_compute_phen_contrib_by_GU(
    bool compute_phen_contrib_by_GU) {
  compute_phen_contrib_by_GU_ = compute_phen_contrib_by_GU;
}

void OutputManager::set_logs(int8_t logs) {
  logs_->set_logs(logs);
}

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void OutputManager::flush() {
  stats_->flush();
}

} // namespace aevol

#endif // AEVOL_OUPUT_MANAGER_H_
