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


#ifndef AEVOL_GRID_CELL_H__
#define AEVOL_GRID_CELL_H__


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
#include "Individual.h"
#include "Habitat.h"

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;





class GridCell
{
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  GridCell(void) = delete;
  GridCell(const GridCell &) = delete;
  GridCell(int16_t x, int16_t y,
               std::unique_ptr<Habitat>&& habitat,
               Individual * indiv);
  GridCell(gzFile backup_file,
               ExpManager * exp_m,
               std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler_);

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~GridCell(void);


  // =================================================================
  //                        Accessors: getters
  // =================================================================
  inline int16_t x() const {return x_;};
  inline int16_t y() const {return y_;};
  inline double compound_amount(void) const;
  inline Individual * get_individual(void) const;

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
  inline void set_individual(Individual * indiv);

  // =================================================================
  //                            Public Methods
  // =================================================================
  void ApplyHabitatVariation();
  void save(gzFile backup_file,
            bool skip_phenotypic_target = false) const;
  void load(gzFile backup_file,
            ExpManager * exp_m,
            std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler);

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
  Individual* individual_ = NULL;

  std::unique_ptr<Habitat> habitat_ = nullptr;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline double GridCell::compound_amount(void) const
{
  return habitat_->compound_amount();
}

inline Individual *GridCell::get_individual(void) const
{
  return individual_;
}

inline double GridCell::get_secreted_amount(void) const
{
  return individual_->get_fitness_by_feature(SECRETION);
}

inline double GridCell::get_metabolic_fitness(void) const
{
  return individual_->get_fitness_by_feature(METABOLISM);
}

inline double GridCell::get_total_fitness(void) const
{
  return individual_->get_fitness();
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void GridCell::set_compound_amount(double compound_amount)
{
  habitat_->set_compound_amount(compound_amount);
}

inline void GridCell::set_individual(Individual * indiv)
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
#endif // AEVOL_GRID_CELL_H__
