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
#ifndef AEVOL_EXP_MANAGER_H_
#define AEVOL_EXP_MANAGER_H_

#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <list>

#include "AeTime.h"
#include "JumpingMT.h"
#include "ExpSetup.h"
#include "OutputManager.h"
#include "World.h"
#include "Observer.h"
#include "ObservableEvent.h"

namespace aevol {
// =================================================================
//                          Class declarations
// =================================================================

/// Allows for high-level experiment management. (This is Aevol's top-level class.)
///
/// An experiment manager allows one to... manage an experiment.
/// It owns a population and an experimental_setup that can be loaded from a
/// pair of aevol binary files (pop and exp_setup)
class ExpManager : public Observer {
 public:
  // =======================================================================
  //                                Constructors
  // =======================================================================
  ExpManager();

  // =======================================================================
  //                                Destructors
  // =======================================================================
  virtual ~ExpManager();

  // =======================================================================
  //                                 Algorithms
  // =======================================================================

  // =======================================================================
  //                           Accessors: getters
  // =======================================================================
  ExpSetup* get_exp_s(void) const { return _exp_s; }
  Selection* get_sel(void) const { return get_exp_s()->get_sel(); }
  OutputManager* get_output_m(void) const { return _output_m; }
  bool quit_signal_received(void) const { return _quit_signal_received; }
  SelectionScheme get_selection_scheme(void) const { return get_sel()->get_selection_scheme(); }
  double get_selection_pressure(void) const { return get_sel()->get_selection_pressure(); }

  // Spatial structure
  World* world(void) const { return world_; }
  GridCell* get_grid_cell(int16_t x, int16_t y) const { return world()->grid(x, y); }
  int16_t get_grid_width(void) const { return world()->width(); }
  int16_t get_grid_height(void) const { return world()->height(); }
  GridCell*** get_pop_grid(void) const { return world()->grid(); }

  // Global settings
  double get_repl_HT_detach_rate(void) const { return get_exp_s()->get_repl_HT_detach_rate();}

  // The ability to own a plasmid is a property of the individuals (allow_plasmids_) because it is used during mutations.
  // However, the experimental setup's member variable _with_plasmids indicates whether plasmids are used
  // because the replication and loading/writting processes need this information.
  // For now when plasmids are used each individual has one and only one plasmid (so these variables should always be
  // equal), this may change in the future, though.
  // Member variable _with_plasmids_HT has been removed because the ability to transfer is evolvable and may thus depend
  // on the plasmid itself

  bool   get_with_plasmids(void) const { return get_exp_s()->get_with_plasmids(); }
  double get_prob_plasmid_HT(void) const { return get_exp_s()->get_prob_plasmid_HT(); }
  double get_tune_donor_ability(void) const { return get_exp_s()->get_tune_donor_ability(); }
  double get_tune_recipient_ability(void) const { return get_exp_s()->get_tune_recipient_ability(); }
  bool   get_swap_GUs(void) const { return get_exp_s()->get_swap_GUs(); }
  bool   get_with_secretion(void) const { return get_exp_s()->get_with_secretion(); }
  double get_secretion_contrib_to_fitness(void) const { return get_exp_s()->get_secretion_contrib_to_fitness(); }
  double get_secretion_cost(void) const { return get_exp_s()->get_secretion_cost(); }

  // Accessors to population stuff
  std::list<Individual*> get_indivs() const { return world()->get_indivs(); }
  int32_t get_nb_indivs(void) const { return world()->get_nb_indivs(); }
  Individual* get_best_indiv(void) const { return world()->get_best_indiv(); }
  Individual* get_indiv_by_id(int32_t id) const;

  // Accessors to output manager stuff
  int64_t	get_backup_step(void) const { return get_output_m()->get_backup_step(); }
  bool get_record_tree(void) const { return get_output_m()->get_record_tree(); }
  int32_t get_tree_step(void) const { return get_output_m()->get_tree_step(); }
  Tree* get_tree(void) const { return get_output_m()->get_tree(); }

  // =======================================================================
  //                          Accessors: setters
  // =======================================================================
  void set_t_end(int64_t _t_end) { t_end = _t_end; }
  void set_HT_ins_rate(double HT_ins_rate) { _exp_s->set_HT_ins_rate(HT_ins_rate); }
  void set_HT_repl_rate(double HT_repl_rate) { _exp_s->set_HT_repl_rate(HT_repl_rate); }

  // =======================================================================
  //                                 Operators
  // =======================================================================

  // =======================================================================
  //                               Public Methods
  // =======================================================================
  void InitializeWorld(int16_t grid_width,
                       int16_t grid_height,
                       std::shared_ptr<JumpingMT> prng,
                       const Habitat& habitat,
                       bool share_phenotypic_target);
  void Save() const;
  void WriteSetupFiles() const;
  void WriteDynamicFiles() const;
  void save_copy(char* dir, int64_t time = 0) const;
  void load(int64_t first_gener, bool verbose = false, bool to_be_run = true) { load(".", first_gener, verbose, to_be_run); }
  void load(const char* dir, int64_t t0,
            bool verbose = false, bool to_be_run = true);
  /// Load an experiment with default files from the current directory
  void load(int64_t t0, char* exp_s_file, char* exp_backup_file, char* sp_struct_file, char* out_p_file,
            bool verbose = false, bool to_be_run = true);
  void run_evolution(void);
  virtual void display(void) {};
  void update_best(void);
  void FillGridWithClones(Individual& dolly) { world_->FillGridWithClones(dolly); }
  void update(Observable& o, ObservableEvent e, void* arg) override {}

 protected:
  // =======================================================================
  //                              Protected Methods
  // =======================================================================
  void step_to_next_generation();

  void load(gzFile& exp_s_file,
            gzFile& exp_backup_file,
            gzFile& sp_struct_file,
            gzFile& out_p_file,
            bool verbose = false,
            bool to_be_run = true);

  void create_missing_directories(const char* dir = ".") const;
  void open_backup_files(gzFile& sel_file,
                         gzFile& sp_struct_file,
                         int64_t t,
                         const char mode[3],
                         const char* dir = ".") const;
  void close_backup_files(gzFile& sel_file,
                          gzFile& sp_struct_file) const;
  void open_setup_files(gzFile& exp_s_gzfile,
                        gzFile& out_p_gzfile,
                        int64_t t,
                        const char mode[3],
                        const char* dir = ".") const;
  void close_setup_files(gzFile& exp_s_gzfile, gzFile& out_p_gzfile) const;

  // =======================================================================
  //                             Protected Attributes
  // =======================================================================
  /// Experimental setup
  ExpSetup* _exp_s;

  /// Spatial structure
  World* world_;

  /// Output manager
  OutputManager* _output_m;

  /// Time step up to which we want to simulate
  int64_t t_end;

  // Set to true when ctrl-Q is received. Will cause the simulation
  // to be ended after the current time step is completed
  bool _quit_signal_received;
};

// ===========================================================================
//                             Getters' definitions
// ===========================================================================

// ===========================================================================
//                             Setters' definitions
// ===========================================================================

// ===========================================================================
//                            Operators' definitions
// ===========================================================================

// ===========================================================================
//                         Inline methods' definition
// ===========================================================================

} // namespace aevol
#endif // AEVOL_EXP_MANAGER_H_
