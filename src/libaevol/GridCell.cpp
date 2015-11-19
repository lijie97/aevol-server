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
//*****************************************************************************




// =================================================================
//                              Includes
// =================================================================
#include "GridCell.h"

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include <iostream>

#include "Individual.h"

#ifdef __X11
#include "Individual_X11.h"
#endif

#ifdef __REGUL
#include "raevol/Individual_R.h"
#include "raevol/Habitat_R.h"
#endif

using std::cout;
using std::endl;


namespace aevol {



//##############################################################################
//                                                                             #
//                              Class GridCell                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
GridCell::GridCell(int16_t x, int16_t y,
                   std::unique_ptr<Habitat>&& habitat,
                   Individual* indiv)
{
  x_ = x;
  y_ = y;
  
  individual_ = indiv;
  habitat_ = std::move(habitat);
}

GridCell::GridCell(gzFile backup_file,
                   ExpManager* exp_m,
                   PhenotypicTargetHandler* phenotypic_target_handler)
{
  load(backup_file, exp_m, phenotypic_target_handler);
}

// =================================================================
//                             Destructors
// =================================================================
GridCell::~GridCell(void)
{
  delete individual_;
}

// =================================================================
//                            Public Methods
// =================================================================
void GridCell::ApplyHabitatVariation() {
  habitat_->ApplyVariation();
}

void GridCell::save(gzFile backup_file,
                        bool skip_phenotypic_target /*=false*/) const
{
  gzwrite(backup_file, &x_, sizeof(x_));
  gzwrite(backup_file, &y_, sizeof(y_));

  #ifndef __REGUL
  habitat_->save(backup_file, skip_phenotypic_target);
  individual_->save(backup_file);
  #else
  (dynamic_cast<Habitat_R*> (habitat_.get()))->save(backup_file, skip_phenotypic_target);
  (dynamic_cast<Individual_R*> (individual_))->save(backup_file);
  #endif
}

void GridCell::load(gzFile backup_file,
                        ExpManager * exp_m,
                        PhenotypicTargetHandler* phenotypic_target_handler)
{
  gzread(backup_file, &x_, sizeof(x_));
  gzread(backup_file, &y_, sizeof(y_));

#if __cplusplus == 201103L
  habitat_ = make_unique<Habitat>(backup_file, phenotypic_target_handler);
#else
  habitat_ = std::make_unique<Habitat>(backup_file, phenotypic_target_handler);
#endif

  individual_ = Individual::CreateIndividual(exp_m, backup_file);

  individual_->set_grid_cell(this);
}

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
