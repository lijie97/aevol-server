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


/*!
  \class ExpManager

  \brief This is aevol's top-level class. It allows for
  high-level experiment management

  An experiment manager allows one to... manage an experiment.
  It owns a population and an experimental_setup that can be loaded from a
  pair of aevol binary files (pop and exp_setup)
*/


// =================================================================
//                              Includes
// =================================================================
#ifndef AEVOL_EXP_MANAGER_H__
#define AEVOL_EXP_MANAGER_H__

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


class ExpManager : public Observer
{
  public :
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
    // TODO: tmp disabled
    // void foreach_indiv(void (*processor)(Individual& indiv)) const;

    // =======================================================================
    //                           Accessors: getters
    // =======================================================================
    inline ExpSetup * get_exp_s(void) const;
    inline Selection * get_sel(void) const;
    inline OutputManager * get_output_m(void) const;

    inline bool quit_signal_received(void) const;

    inline SelectionScheme get_selection_scheme(void) const;
    inline double get_selection_pressure(void) const;

    // ------------------------------------------------------ Spatial structure
    inline World*             world(void) const;
    inline GridCell *      get_grid_cell(int16_t x, int16_t y) const;
    inline int16_t            get_grid_width(void) const;
    inline int16_t            get_grid_height(void) const;
    inline GridCell ***    get_pop_grid(void) const;

    // -------------------------------------------------------- Global settings
    inline bool   get_with_HT(void) const;
    inline bool   get_repl_HT_with_close_points(void) const;
    inline double get_HT_ins_rate(void) const;
    inline double get_HT_repl_rate(void) const;
    inline double get_repl_HT_detach_rate(void) const;

    // The ability to own a plasmid is a property of the individuals (_allow_plasmids) because it is used during mutations
    // However there is also a property of the experimental setup (_with_plasmids) that indicates whether plasmids are used because we need this during replication and during loading/writting
    // For now when plasmids are used each individual has one and only one plasmid (so these variables should always be equals), however this may change in the future
    // There is no longer property _with_plasmids_HT because the ability to transfer is evolvable and thus may depend on the plasmid itself

    inline bool   get_with_plasmids(void) const;
    inline double get_prob_plasmid_HT(void) const;
    inline double get_tune_donor_ability(void) const;
    inline double get_tune_recipient_ability(void) const;
    inline double get_donor_cost(void) const;
    inline double get_recipient_cost(void) const;
    inline bool   get_swap_GUs(void) const;

    inline bool   get_with_secretion(void) const;
    inline double get_secretion_contrib_to_fitness(void) const;
    inline double get_secretion_cost(void) const;

    //~ inline bool   get_with_alignments(void) const;

    // Accessors to population stuff
    inline std::list<Individual *> get_indivs() const;
    inline int32_t get_nb_indivs(void) const;

    inline Individual * get_best_indiv(void) const;
    Individual* get_indiv_by_id(int32_t id) const;
    // inline Individual* get_indiv_by_rank(int32_t rank) const;

    // Accessors to output manager stuff
    inline int64_t	get_backup_step(void) const;
    inline int64_t	get_big_backup_step(void) const;
    inline bool         get_record_tree(void) const;
    inline int32_t      get_tree_step(void) const;
    inline Tree*     get_tree(void) const;

    // =======================================================================
    //                          Accessors: setters
    // =======================================================================
    inline void set_t_end(int64_t _t_end) { t_end = _t_end; };
    //~ inline void set_min_genome_length(int32_t min_genome_length);
    //~ inline void set_max_genome_length(int32_t max_genome_length);
    inline void set_with_HT(bool with_HT);
    inline void set_repl_HT_with_close_points (bool repl_HT_with_close_points) ;
    inline void set_HT_ins_rate(double HT_ins_rate) ;
    inline void set_HT_repl_rate(double HT_repl_rate) ;
    inline void set_repl_HT_detach_rate(double repl_HT_detach_rate) ;

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
    inline void load(int64_t first_gener,
                     bool verbose = false, bool to_be_run = true);
    void load(const char* dir, int64_t t0,
        bool verbose = false, bool to_be_run = true);
    void load(int64_t t0,
              char* exp_s_file,
              char* exp_backup_file,
              char* sp_struct_file,
              char* out_p_file,
              bool verbose = false,
              bool to_be_run = true);
    void run_evolution(void);
    virtual void display(void) {};
    void update_best(void);

    void FillGridWithClones(Individual & dolly) {
      world_->FillGridWithClones(dolly);
    }

    void update(Observable& o, ObservableEvent e, void* arg) override {

    };




  protected :
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
    // ---------------------------------------------------- Experimental setup
    ExpSetup* _exp_s;

    // ----------------------------------------------------- Spatial structure
    World* world_;

    // -------------------------------------------------------- Output manager
    OutputManager* _output_m;

    // ----------------------------- Time step up to which we want to simulate
    int64_t t_end;

    // Set to true when ctrl-Q is received. Will cause the simulation
    // to be ended after the current time step is completed
    bool _quit_signal_received;
};


// ===========================================================================
//                             Getters' definitions
// ===========================================================================
inline ExpSetup *ExpManager::get_exp_s(void) const
{
  return _exp_s;
}

inline OutputManager *ExpManager::get_output_m(void) const
{
  return _output_m;
}

inline bool ExpManager::quit_signal_received(void) const
{
  return _quit_signal_received;
}

inline Selection *ExpManager::get_sel(void) const
{
  return get_exp_s()->get_sel();
}

inline SelectionScheme ExpManager::get_selection_scheme(void) const
{
  return get_sel()->get_selection_scheme();
}

inline double ExpManager::get_selection_pressure(void) const
{
  return get_sel()->get_selection_pressure();
}

// Global settings
inline World* ExpManager::world(void) const
{
  return world_;
}

inline GridCell *ExpManager::get_grid_cell(int16_t x, int16_t y) const
{
  return world()->grid(x, y);
}

inline int16_t ExpManager::get_grid_width(void) const
{
  return world()->width();
}

inline int16_t ExpManager::get_grid_height(void) const
{
  return world()->height();
}

inline GridCell ***ExpManager::get_pop_grid(void) const
{
  return world()->grid();
}

inline bool ExpManager::get_with_HT(void) const
{
  return get_exp_s()->get_with_HT();
}

inline bool ExpManager::get_repl_HT_with_close_points (void) const
{
  return get_exp_s()->get_repl_HT_with_close_points();
}

inline double ExpManager::get_HT_ins_rate(void) const
{
  return get_exp_s()->get_HT_ins_rate();
}

inline double ExpManager::get_HT_repl_rate(void) const
{
  return get_exp_s()->get_HT_repl_rate();
}

inline double ExpManager::get_repl_HT_detach_rate(void) const
{
  return get_exp_s()->get_repl_HT_detach_rate();
}

inline bool ExpManager::get_with_plasmids(void) const
{
  return get_exp_s()->get_with_plasmids();
}

inline double ExpManager::get_prob_plasmid_HT(void) const
{
  return get_exp_s()->get_prob_plasmid_HT();
}

inline double ExpManager::get_tune_donor_ability(void) const
{
  return get_exp_s()->get_tune_donor_ability();
}

inline double ExpManager::get_tune_recipient_ability(void) const
{
  return get_exp_s()->get_tune_recipient_ability();
}

inline double ExpManager::get_donor_cost(void) const
{
  return get_exp_s()->get_donor_cost();
}

inline double ExpManager::get_recipient_cost(void) const
{
  return get_exp_s()->get_recipient_cost();
}

inline bool ExpManager::get_swap_GUs(void) const
{
  return get_exp_s()->get_swap_GUs();
}

inline bool ExpManager::get_with_secretion(void) const
{
  return get_exp_s()->get_with_secretion();
}

inline double ExpManager::get_secretion_contrib_to_fitness(void) const
{
  return get_exp_s()->get_secretion_contrib_to_fitness();
}

inline double ExpManager::get_secretion_cost(void) const
{
  return get_exp_s()->get_secretion_cost();
}

//~ inline bool ExpManager::get_with_alignments(void) const
//~ {
  //~ return _exp_s->get_with_alignments();
//~ }


// Accessors to population stuff
inline int32_t ExpManager::get_nb_indivs(void) const
{
  return world()->get_nb_indivs();
}

inline Individual *ExpManager::get_best_indiv(void) const
{
  return world()->get_best_indiv();
}

inline std::list<Individual *> ExpManager::get_indivs() const
{
  return world()->get_indivs();
}


// Accessors to output manager stuff
inline int64_t ExpManager::get_backup_step(void) const
{
	return get_output_m()->get_backup_step();
}

inline int64_t ExpManager::get_big_backup_step(void) const
{
	return get_output_m()->get_big_backup_step();
}

inline bool ExpManager::get_record_tree(void) const
{
	return get_output_m()->get_record_tree();
}

inline int32_t ExpManager::get_tree_step(void) const
{
	return get_output_m()->get_tree_step();
}

inline Tree *ExpManager::get_tree(void) const
{
	return get_output_m()->get_tree();
}

// ===========================================================================
//                             Setters' definitions
// ===========================================================================
// Global constraints
//~ inline void ExpManager::set_min_genome_length(int32_t min_genome_length)
//~ {
  //~ _exp_s->set_min_genome_length(min_genome_length);
//~ }

//~ inline void ExpManager::set_max_genome_length(int32_t max_genome_length)
//~ {
  //~ _exp_s->set_max_genome_length(max_genome_length);
//~ }


inline void ExpManager::set_with_HT(bool with_HT)
{
  _exp_s->set_with_HT(with_HT);
}

inline void ExpManager::set_repl_HT_with_close_points (bool repl_HT_with_close_points)
{
  _exp_s->set_repl_HT_with_close_points(repl_HT_with_close_points);
}

inline void ExpManager::set_HT_ins_rate(double HT_ins_rate)
{
  _exp_s->set_HT_ins_rate(HT_ins_rate);
}

inline void ExpManager::set_HT_repl_rate(double HT_repl_rate)
{
  _exp_s->set_HT_repl_rate(HT_repl_rate);
}

inline void ExpManager::set_repl_HT_detach_rate(double repl_HT_detach_rate)
{
  _exp_s->set_repl_HT_detach_rate(repl_HT_detach_rate);
}

// ===========================================================================
//                            Operators' definitions
// ===========================================================================

// ===========================================================================
//                         Inline methods' definition
// ===========================================================================
/*!
  \brief Load an experiment with default files from the current directory
 */
inline void ExpManager::load(int64_t first_gener,
                             bool verbose /*= false*/,
                             bool to_be_run /*= true*/)
{
  load(".", first_gener, verbose, to_be_run);
}

} // namespace aevol
#endif // AEVOL_EXP_MANAGER_H__
