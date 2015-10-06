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


#ifndef AEVOL_OUPUT_MANAGER_H__
#define AEVOL_OUPUT_MANAGER_H__


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
  inline int64_t	get_backup_step() const;
  inline int64_t	get_big_backup_step() const;
  
  // Tree
  inline bool get_record_tree() const;
  inline int64_t get_tree_step() const;
  inline Tree* get_tree() const;

  // Logs
  inline FILE* get_log(LogType log_type) const;
  inline bool is_logged(LogType log_type) const;

  // Stats
  inline bool get_compute_phen_contrib_by_GU() const;

  // =================================================================
  //                        Accessors: setters
  // =================================================================
  inline void set_backup_step(int64_t backup_step);
  inline void set_big_backup_step(int64_t big_backup_step);
  inline void init_tree(ExpManager* exp_m, int64_t _tree_step);
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
  void CopyStats(const std::string& outdir, int64_t time) const;
  void load(gzFile file, bool verbose, bool to_be_run);
  void write_current_generation_outputs(void) const;
  inline void flush(void);

  static int64_t get_last_gener();

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  void write_tree(void) const;

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  ExpManager* _exp_m;
  
  // Backups
  int64_t _backup_step;
  int64_t _big_backup_step;
  
  // Stats
  Stats* _stats;
  bool _compute_phen_contrib_by_GU;

  // Tree
  bool _record_tree;
  Tree* _tree;
  
  // Dumps
  bool _make_dumps;
  int64_t _dump_step;
  Dump* _dump;
  
  // Logs
  Logging* _logs;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// Backup
inline int64_t OutputManager::get_backup_step() const {
  return _backup_step;
}

inline int64_t OutputManager::get_big_backup_step() const {
  return _big_backup_step;
}

// Tree
inline bool OutputManager::get_record_tree() const {
  return _record_tree;
}

inline int64_t OutputManager::get_tree_step() const {
  return _tree->get_tree_step();
}

inline Tree *OutputManager::get_tree() const {
  return _tree;
}

// Logs
inline FILE* OutputManager::get_log(LogType log_type) const {
  return _logs->get_log(log_type);
}

inline bool  OutputManager::is_logged(LogType log_type) const {
  return _logs->is_logged(log_type);
}

// Stats
inline bool OutputManager::get_compute_phen_contrib_by_GU() const {
  return _compute_phen_contrib_by_GU;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void OutputManager::set_backup_step(int64_t backup_step) {
  _backup_step = backup_step;
}

inline void OutputManager::set_big_backup_step(int64_t big_backup_step) {
  _big_backup_step = big_backup_step;
}

inline void OutputManager::init_tree(ExpManager * exp_m, int64_t _tree_step) {
  _record_tree = true;
  _tree = new Tree(exp_m, _tree_step);
}

inline void OutputManager::set_dump_step(int64_t dump_step) {
  _make_dumps = true;
  _dump_step  = dump_step;
}

inline void OutputManager::set_compute_phen_contrib_by_GU(
    bool compute_phen_contrib_by_GU) {
  _compute_phen_contrib_by_GU = compute_phen_contrib_by_GU;
}

inline void OutputManager::set_logs (int8_t logs) {
  _logs->set_logs(logs);
}

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void OutputManager::flush(void) {
  _stats->flush();
}

} // namespace aevol

#endif // AEVOL_OUPUT_MANAGER_H__
