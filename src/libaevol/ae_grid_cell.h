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


#ifndef __AE_GRID_CELL_H__
#define __AE_GRID_CELL_H__


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include "ae_individual.h"
#include "habitat.h"

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;





class ae_grid_cell
{
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  ae_grid_cell(void) = delete;
  ae_grid_cell(const ae_grid_cell&) = delete;
  ae_grid_cell(int16_t x, int16_t y,
               std::unique_ptr<Habitat>&& habitat,
               ae_individual* indiv);
  ae_grid_cell(gzFile backup_file,
               ae_exp_manager* exp_m,
               std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler_);

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~ae_grid_cell(void);


  // =================================================================
  //                        Accessors: getters
  // =================================================================
  inline int16_t x() const {return x_;};
  inline int16_t y() const {return y_;};
  inline double compound_amount(void) const;
  inline ae_individual* get_individual(void) const;

  inline double get_secreted_amount(void) const;
  inline double get_metabolic_fitness(void) const;
  inline double get_total_fitness(void) const;

  const Habitat& habitat() const {
    return *habitat_;
  }
  const PhenotypicTarget& phenotypic_target() const {
    return habitat_->phenotypic_target();
  }

  // =================================================================
  //                        Accessors: setters
  // =================================================================
  inline void set_compound_amount(double compound_amount);
  inline void set_individual(ae_individual* indiv);

  // =================================================================
  //                            Public Methods
  // =================================================================
  void save(gzFile backup_file,
            bool skip_phenotypic_target = false) const;
  void load(gzFile backup_file,
            ae_exp_manager* exp_m,
            std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler);

  // =================================================================
  //                           Public Attributes
  // =================================================================





 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  // Position on the grid
  int16_t x_;
  int16_t y_;

  // Pointer to the individual in this cell
  ae_individual* individual_ = NULL;

  std::unique_ptr<Habitat> habitat_;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline double ae_grid_cell::compound_amount(void) const
{
  return habitat_->compound_amount();
}

inline ae_individual* ae_grid_cell::get_individual(void) const
{
  return individual_;
}

inline double ae_grid_cell::get_secreted_amount(void) const
{
  return individual_->get_fitness_by_feature(SECRETION);
}

inline double ae_grid_cell::get_metabolic_fitness(void) const
{
  return individual_->get_fitness_by_feature(METABOLISM);
}

inline double ae_grid_cell::get_total_fitness(void) const
{
  return individual_->get_fitness();
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void ae_grid_cell::set_compound_amount(double compound_amount)
{
  habitat_->set_compound_amount(compound_amount);
}

inline void ae_grid_cell::set_individual(ae_individual* indiv)
{
  individual_ = indiv;
  if (individual_->get_grid_cell() != this)
  {
    individual_->set_grid_cell(this);
  }
}



// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // __AE_GRID_CELL_H__
