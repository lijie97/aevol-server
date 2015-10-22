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
    double get_repl_HT_detach_rate() const { return _repl_HT_detach_rate; }
    
    // --------------------------------------------------------------- Plasmids
    // See comments in ExpManager.h on how plasmids are handled
    bool   get_with_plasmids() const { return _with_plasmids; }
    double get_prob_plasmid_HT() const { return _prob_plasmid_HT; }
    double get_tune_donor_ability() const { return _tune_donor_ability; }
    double get_tune_recipient_ability() const { return _tune_recipient_ability; }
    bool   get_swap_GUs() const { return _swap_GUs; }
    
    // -------------------------------------------------------------- Secretion
    bool   get_with_secretion() const { return _with_secretion; }
    double get_secretion_contrib_to_fitness() const { return _secretion_contrib_to_fitness; }
    double get_secretion_cost() const { return _secretion_cost; }
  
    // =======================================================================
    //                         Accessors: setters
    // =======================================================================
    // --------------------------------------------------------------- Transfer
    void set_with_HT(bool with_HT) { _with_HT = with_HT; }
    void set_repl_HT_with_close_points(bool repl_HT_with_close_points) { _repl_HT_with_close_points = repl_HT_with_close_points; }
    void set_HT_ins_rate(double HT_ins_rate) { _HT_ins_rate = HT_ins_rate; }
    void set_HT_repl_rate(double HT_repl_rate) { _HT_repl_rate = HT_repl_rate; }
    void set_repl_HT_detach_rate(double repl_HT_detach_rate) { _repl_HT_detach_rate = repl_HT_detach_rate; }
  
    // --------------------------------------------------------------- Plasmids
    void set_with_plasmids(bool with_p) { _with_plasmids = with_p; }
    void set_prob_plasmid_HT(double prob_p_HT) { _prob_plasmid_HT = prob_p_HT; }
    void set_tune_donor_ability(double tune_donor_ability) { _tune_donor_ability = tune_donor_ability; }
    void set_tune_recipient_ability(double tune_recipient_ability) { _tune_recipient_ability = tune_recipient_ability; }
    void set_donor_cost(double donor_cost) { _donor_cost = donor_cost; }
    void set_recipient_cost(double recipient_cost) { _recipient_cost = recipient_cost; }
    void set_swap_GUs(bool swap_GUs) { _swap_GUs = swap_GUs; }
    
    // -------------------------------------------------------------- Secretion
    void set_with_secretion(bool with_secretion) { _with_secretion = with_secretion; }
    void set_secretion_contrib_to_fitness(double secretion_contrib) { _secretion_contrib_to_fitness = secretion_contrib; }
    void set_secretion_cost(double secretion_cost) { _secretion_cost = secretion_cost; }

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
    bool    _with_HT;
    bool    _repl_HT_with_close_points;
    double  _HT_ins_rate;
    double  _HT_repl_rate;
    double  _repl_HT_detach_rate;
  
    // --------------------------------------------------- Plasmids parameters
    bool    _with_plasmids;
    double  _prob_plasmid_HT; // Base transfer ability independent of evolvable donor and recipient ability
    double  _tune_donor_ability; // How much the individuals can tune their ability to send plasmids
    double  _tune_recipient_ability; // How much the individuals can tune their ability to receive plasmids
    double  _donor_cost;
    double  _recipient_cost;
    bool    _swap_GUs; // Whether plasmid HT is uni- or bidirectional
    
    // -------------------------------------------------- Secretion parameters
    bool    _with_secretion;
    double  _secretion_contrib_to_fitness;
    double  _secretion_cost;
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
