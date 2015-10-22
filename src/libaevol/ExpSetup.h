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

#ifndef AEVOL_EXP_SETUP_H__
#define AEVOL_EXP_SETUP_H__

// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdlib>

#include "Selection.h"
#include "Stats.h"
#include "Logging.h"
#include "Tree.h"
#include "Dump.h"
#include "JumpingMT.h"

namespace aevol {

// ===========================================================================
//                          Class declarations
// ===========================================================================
class ae_param_loader;

class ExpSetup {
  friend class ExpManager;
  
  public:
    // =======================================================================
    //                             Constructors
    // =======================================================================
    ExpSetup() = delete;
    ExpSetup(const ExpSetup&) = delete;
    ExpSetup(ExpManager * exp_m);

    // =======================================================================
    //                             Destructors
    // =======================================================================
    virtual ~ExpSetup();

    // =======================================================================
    //                         Accessors: getters
    // =======================================================================
    // ----------------------------------------------------- Selection context
    Selection * get_sel() const { return _sel; }
  
    // --------------------------------------------------------------- Transfer
    double get_repl_HT_detach_rate() const { return repl_HT_detach_rate_; }
    
    // --------------------------------------------------------------- Plasmids
    // See comments in ExpManager.h on how plasmids are handled
    bool   get_with_plasmids() const { return with_plasmids_; }
    double get_prob_plasmid_HT() const { return prob_plasmid_HT_; }
    double get_tune_donor_ability() const { return tune_donor_ability_; }
    double get_tune_recipient_ability() const { return tune_recipient_ability_; }
    bool   get_swap_GUs() const { return swap_GUs_; }
    
    // -------------------------------------------------------------- Secretion
    bool   get_with_secretion() const { return with_secretion_; }
    double get_secretion_contrib_to_fitness() const { return secretion_contrib_to_fitness_; }
    double get_secretion_cost() const { return secretion_cost_; }
  
    // =======================================================================
    //                         Accessors: setters
    // =======================================================================
    // --------------------------------------------------------------- Transfer
    void set_with_HT(bool with_HT) { with_HT_ = with_HT; }
    void set_repl_HT_with_close_points(bool repl_HT_with_close_points) { repl_HT_with_close_points_ = repl_HT_with_close_points; }
    void set_HT_ins_rate(double HT_ins_rate) { HT_ins_rate_ = HT_ins_rate; }
    void set_HT_repl_rate(double HT_repl_rate) { HT_repl_rate_ = HT_repl_rate; }
    void set_repl_HT_detach_rate(double repl_HT_detach_rate) { repl_HT_detach_rate_ = repl_HT_detach_rate; }
  
    // --------------------------------------------------------------- Plasmids
    void set_with_plasmids(bool with_p) { with_plasmids_ = with_p; }
    void set_prob_plasmid_HT(double prob_p_HT) { prob_plasmid_HT_ = prob_p_HT; }
    void set_tune_donor_ability(double tune_donor_ability) { tune_donor_ability_ = tune_donor_ability; }
    void set_tune_recipient_ability(double tune_recipient_ability) { tune_recipient_ability_ = tune_recipient_ability; }
    void set_donor_cost(double donor_cost) { donor_cost_ = donor_cost; }
    void set_recipient_cost(double recipient_cost) { recipient_cost_ = recipient_cost; }
    void set_swap_GUs(bool swap_GUs) { swap_GUs_ = swap_GUs; }
    
    // -------------------------------------------------------------- Secretion
    void set_with_secretion(bool with_secretion) { with_secretion_ = with_secretion; }
    void set_secretion_contrib_to_fitness(double secretion_contrib) { secretion_contrib_to_fitness_ = secretion_contrib; }
    void set_secretion_cost(double secretion_cost) { secretion_cost_ = secretion_cost; }

    // =======================================================================
    //                            Public Methods
    // =======================================================================
    void write_setup_file(gzFile exp_setup_file) const;
    void save(gzFile backup_file) const;
    void load(gzFile setup_file, gzFile backup_file, bool verbose);
    /// Make the individuals reproduce
    void step_to_next_generation() { _sel->step_to_next_generation(); }

    // =======================================================================
    //                           Public Attributes
    // =======================================================================
  
  protected :
    // =======================================================================
    //                           Protected Methods
    // =======================================================================
    virtual void display() {};
  
    // =======================================================================
    //                          Protected Attributes
    // =======================================================================
    ExpManager* _exp_m;
      
    // ----------------------------------------------------- Selection context
    Selection* _sel;
    
    // --------------------------------------------------- Transfer parameters
    bool   with_HT_;
    bool   repl_HT_with_close_points_;
    double HT_ins_rate_;
    double HT_repl_rate_;
    double repl_HT_detach_rate_;
  
    // --------------------------------------------------- Plasmids parameters
    bool   with_plasmids_;
    double prob_plasmid_HT_; // Base transfer ability independent of evolvable donor and recipient ability
    double tune_donor_ability_; // How much the individuals can tune their ability to send plasmids
    double tune_recipient_ability_; // How much the individuals can tune their ability to receive plasmids
    double donor_cost_;
    double recipient_cost_;
    bool   swap_GUs_; // Whether plasmid HT is uni- or bidirectional
    
    // -------------------------------------------------- Secretion parameters
    bool   with_secretion_;
    double secretion_contrib_to_fitness_;
    double secretion_cost_;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// =====================================================================
//                           Setters' definitions
// =====================================================================

// =====================================================================
//                       functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_EXP_SETUP_H__
