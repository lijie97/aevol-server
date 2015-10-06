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




 
class ExpSetup
{
  friend class ExpManager;
  
  public :
  
    // =======================================================================
    //                             Constructors
    // =======================================================================
    ExpSetup() = delete;
    ExpSetup(const ExpSetup&) = delete;
    ExpSetup(ExpManager * exp_m);
  
  
    // =======================================================================
    //                             Destructors
    // =======================================================================
    virtual ~ExpSetup(void);
  
  
    // =======================================================================
    //                         Accessors: getters
    // =======================================================================
    // ----------------------------------------------------- Selection context
    inline Selection * get_sel(void) const;
  
    // --------------------------------------------------------------- Transfer
    inline bool   get_with_HT(void) const;
    inline bool   get_repl_HT_with_close_points (void) const;
    inline double get_HT_ins_rate(void) const;
    inline double get_HT_repl_rate(void) const;
    inline double get_repl_HT_detach_rate(void) const;
    
    // --------------------------------------------------------------- Plasmids
    // See comments in ExpManager.h on how plasmids are handled
    inline bool   get_with_plasmids(void) const;
    inline double get_prob_plasmid_HT(void) const;
    inline double get_tune_donor_ability(void) const;
    inline double get_tune_recipient_ability(void) const;
    inline double get_donor_cost(void) const;
    inline double get_recipient_cost(void) const;
    inline bool   get_swap_GUs(void) const;
    
    // -------------------------------------------------------------- Secretion
    inline bool   get_with_secretion(void) const;
    inline double get_secretion_contrib_to_fitness(void) const;
    inline double get_secretion_cost(void) const;
  
  
    // =======================================================================
    //                         Accessors: setters
    // =======================================================================
    // --------------------------------------------------------------- Transfer
    inline void set_with_HT(bool with_HT);
    inline void set_repl_HT_with_close_points(bool repl_HT_with_close_points);
    inline void set_HT_ins_rate(double HT_ins_rate);
    inline void set_HT_repl_rate(double HT_repl_rate);
    inline void set_repl_HT_detach_rate(double repl_HT_detach_rate);
  
    // --------------------------------------------------------------- Plasmids
    inline void set_with_plasmids(bool with_p);
    //~ inline void set_nb_plasmid_HT(int16_t nb_p_HT);
    inline void set_prob_plasmid_HT(double prob_p_HT);
    inline void set_tune_donor_ability(double tune_donor_ability);
    inline void set_tune_recipient_ability(double tune_recipient_ability);
    inline void set_donor_cost(double donor_cost);
    inline void set_recipient_cost(double recipient_cost);
    inline void set_swap_GUs(bool swap_GUs);
    
    // -------------------------------------------------------------- Secretion
    inline void set_with_secretion(bool with_secretion);
    inline void set_secretion_contrib_to_fitness(double secretion_contrib);
    inline void set_secretion_cost(double secretion_cost);
  
  
    // =======================================================================
    //                            Public Methods
    // =======================================================================
    void write_setup_file(gzFile exp_setup_file) const;
    void save(gzFile backup_file) const;
    void load(gzFile setup_file, gzFile backup_file, bool verbose);
    
    inline void step_to_next_generation(void);
    
    
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
inline Selection *ExpSetup::get_sel(void) const
{
  return _sel;
}

inline bool ExpSetup::get_with_HT(void) const
{
  return _with_HT;
}

inline bool ExpSetup::get_repl_HT_with_close_points(void) const
{
  return _repl_HT_with_close_points;
}

inline double ExpSetup::get_HT_ins_rate(void) const
{
  return _HT_ins_rate;
}

inline double ExpSetup::get_HT_repl_rate(void) const
{
  return _HT_repl_rate;
}

inline double ExpSetup::get_repl_HT_detach_rate(void) const
{
  return _repl_HT_detach_rate;
}

inline bool ExpSetup::get_with_plasmids(void) const
{
  return _with_plasmids;
}

inline double ExpSetup::get_prob_plasmid_HT(void) const
{
  return _prob_plasmid_HT;
}

inline double ExpSetup::get_tune_donor_ability(void) const
{
  return _tune_donor_ability;
}

inline double ExpSetup::get_tune_recipient_ability(void) const
{
  return _tune_recipient_ability;
}

inline double ExpSetup::get_donor_cost(void) const
{
  return _donor_cost;
}

inline double ExpSetup::get_recipient_cost(void) const
{
  return _recipient_cost;
}

inline bool   ExpSetup::get_swap_GUs(void) const
{
  return _swap_GUs;
}

inline bool ExpSetup::get_with_secretion(void) const
{
  return _with_secretion;
}

inline double ExpSetup::get_secretion_contrib_to_fitness(void) const
{
  return _secretion_contrib_to_fitness;
}

inline double ExpSetup::get_secretion_cost(void) const
{
  return _secretion_cost;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
// --------------------------------------------------------------- Transfer
inline void ExpSetup::set_with_HT(bool with_HT)
{
  _with_HT = with_HT;
}

inline void ExpSetup::set_repl_HT_with_close_points(bool repl_HT_with_close_points)
{
  _repl_HT_with_close_points = repl_HT_with_close_points;
}

inline void ExpSetup::set_HT_ins_rate(double HT_ins_rate)
{
  _HT_ins_rate = HT_ins_rate;
}

inline void ExpSetup::set_HT_repl_rate(double HT_repl_rate)
{
  _HT_repl_rate = HT_repl_rate;
}

inline void ExpSetup::set_repl_HT_detach_rate(double repl_HT_detach_rate)
{
  _repl_HT_detach_rate = repl_HT_detach_rate;
}

inline void ExpSetup::set_with_plasmids(bool with_p)
{
  _with_plasmids = with_p;
}

inline void ExpSetup::set_prob_plasmid_HT(double prob_p_HT)
{
  _prob_plasmid_HT = prob_p_HT;
}

inline void ExpSetup::set_tune_donor_ability(double tune_donor_ability)
{
  _tune_donor_ability = tune_donor_ability;
}

inline void ExpSetup::set_tune_recipient_ability(double tune_recipient_ability)
{
  _tune_recipient_ability = tune_recipient_ability;
}

inline void ExpSetup::set_donor_cost(double donor_cost)
{
  _donor_cost = donor_cost;
}

inline void ExpSetup::set_recipient_cost(double recipient_cost)
{
  _recipient_cost = recipient_cost;
}

inline void ExpSetup::set_swap_GUs(bool swap_GUs)
{
  _swap_GUs = swap_GUs;
}

// -------------------------------------------------------------- Secretion
inline void ExpSetup::set_with_secretion(bool with_secretion)
{
  _with_secretion = with_secretion;
}

inline void ExpSetup::set_secretion_contrib_to_fitness(double secretion_contrib)
{
  _secretion_contrib_to_fitness = secretion_contrib;
}

inline void ExpSetup::set_secretion_cost(double secretion_cost)
{
  _secretion_cost = secretion_cost;
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ExpSetup::step_to_next_generation(void)
{ 
  // Make the individuals reproduce
  _sel->step_to_next_generation();
}

} // namespace aevol

#endif // AEVOL_EXP_SETUP_H__
