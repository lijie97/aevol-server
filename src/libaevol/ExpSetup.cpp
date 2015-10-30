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
#include <cstdio>
#include <zlib.h>

// =================================================================
//                            Project Files
// =================================================================
#include "ExpSetup.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                              Class ExpSetup                                 #
//                                                                             #
//##############################################################################

// ===========================================================================
//                         Definition of static attributes
// ===========================================================================

// ===========================================================================
//                                  Constructors
// ===========================================================================

// ===========================================================================
//                                 Destructor
// ===========================================================================

// ===========================================================================
//                                 Public Methods
// ===========================================================================

template<typename Field>
void gzwrite(gzFile file, Field field) {
  gzwrite(file, &field, static_cast<unsigned>(sizeof(field)));
}

/*!
*/
void ExpSetup::write_setup_file(gzFile exp_setup_file) const
{
  // --------------------------------------------------------------- Transfer
  int8_t tmp_with_HT = with_HT_;
  gzwrite(exp_setup_file, &tmp_with_HT, sizeof(tmp_with_HT));
  int8_t tmp_repl_HT_with_close_points = repl_HT_with_close_points_;
  gzwrite(exp_setup_file, &tmp_repl_HT_with_close_points, sizeof(tmp_repl_HT_with_close_points));
  if (with_HT_)
  {
    gzwrite(exp_setup_file, &HT_ins_rate_,  sizeof(HT_ins_rate_));
    gzwrite(exp_setup_file, &HT_repl_rate_, sizeof(HT_repl_rate_));
  }
  if(repl_HT_with_close_points_)
  {
    gzwrite(exp_setup_file, &repl_HT_detach_rate_,  sizeof(repl_HT_detach_rate_));
  }

  // --------------------------------------------------------------- Plasmids
  int8_t tmp_with_plasmids = with_plasmids();
  gzwrite(exp_setup_file, &tmp_with_plasmids, sizeof(tmp_with_plasmids));
  if (tmp_with_plasmids)
  {
    gzwrite(exp_setup_file, &prob_plasmid_HT_,  sizeof(prob_plasmid_HT_));
    gzwrite(exp_setup_file, &tune_donor_ability_,  sizeof(tune_donor_ability_));
    gzwrite(exp_setup_file, &tune_recipient_ability_,  sizeof(tune_recipient_ability_));
    gzwrite(exp_setup_file, &donor_cost_,  sizeof(donor_cost_));
    gzwrite(exp_setup_file, &recipient_cost_,  sizeof(recipient_cost_));
    int8_t tmp_swap_GUs = swap_GUs_;
    gzwrite(exp_setup_file, &tmp_swap_GUs, sizeof(tmp_swap_GUs));
  }

  // -------------------------------------------------------------- Secretion
  int8_t tmp_with_secretion = with_secretion_;
  gzwrite(exp_setup_file, &tmp_with_secretion, sizeof(tmp_with_secretion));
  gzwrite(exp_setup_file, &secretion_contrib_to_fitness_, sizeof(secretion_contrib_to_fitness_));
  gzwrite(exp_setup_file, &secretion_cost_, sizeof(secretion_cost_));

  sel()->write_setup_file(exp_setup_file);
}

void ExpSetup::load(gzFile setup_file, gzFile backup_file, bool verbose)
{
  // -------------------------------------------- Retrieve transfer parameters
  int8_t tmp_with_HT;
  gzread(setup_file, &tmp_with_HT, sizeof(tmp_with_HT));
  with_HT_ = static_cast<bool>(tmp_with_HT);
  int8_t tmp_repl_HT_with_close_points;
  gzread(setup_file, &tmp_repl_HT_with_close_points, sizeof(tmp_repl_HT_with_close_points));
  repl_HT_with_close_points_ = static_cast<bool>(tmp_repl_HT_with_close_points);
  if (with_HT_)
  {
    gzread(setup_file, &HT_ins_rate_,  sizeof(HT_ins_rate_));
    gzread(setup_file, &HT_repl_rate_, sizeof(HT_repl_rate_));
  }
   if(repl_HT_with_close_points_)
  {
    gzread(setup_file, &repl_HT_detach_rate_,  sizeof(repl_HT_detach_rate_));
  }


  // -------------------------------------------- Retrieve plasmid parameters
  int8_t tmp_with_plasmids;
  gzread(setup_file, &tmp_with_plasmids, sizeof(tmp_with_plasmids));
  with_plasmids_ = static_cast<bool>(tmp_with_plasmids);
  if (with_plasmids_)
  {
    gzread(setup_file, &prob_plasmid_HT_,  sizeof(prob_plasmid_HT_));
    gzread(setup_file, &tune_donor_ability_,  sizeof(tune_donor_ability_));
    gzread(setup_file, &tune_recipient_ability_,  sizeof(tune_recipient_ability_));
    gzread(setup_file, &donor_cost_,  sizeof(donor_cost_));
    gzread(setup_file, &recipient_cost_,  sizeof(recipient_cost_));
    int8_t tmp_swap_GUs;
    gzread(setup_file, &tmp_swap_GUs, sizeof(tmp_swap_GUs));
    swap_GUs_ = static_cast<bool>(tmp_swap_GUs);
  }

  // ------------------------------------------ Retrieve secretion parameters
  int8_t tmp_with_secretion;
  gzread(setup_file, &tmp_with_secretion, sizeof(tmp_with_secretion));
  with_secretion_ = static_cast<bool>(tmp_with_secretion);
  gzread(setup_file, &secretion_contrib_to_fitness_, sizeof(secretion_contrib_to_fitness_));
  gzread(setup_file, &secretion_cost_, sizeof(secretion_cost_));

  // ---------------------------------------------- Retrieve selection context
  sel()->load(setup_file, backup_file, verbose);
}




// ===========================================================================
//                                Protected Methods
// ===========================================================================
} // namespace aevol
