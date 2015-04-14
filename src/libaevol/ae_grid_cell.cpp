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
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include "ae_grid_cell.h"
#include "ae_individual.h"
#include "ae_individual_X11.h"

namespace aevol {



//##############################################################################
//                                                                             #
//                              Class ae_grid_cell                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_grid_cell::ae_grid_cell(int16_t x, int16_t y, ae_individual* indiv)
{
  x_ = x;
  y_ = y;
  
  individual_ = indiv;
  habitat_ = new Habitat();
}

ae_grid_cell::ae_grid_cell(gzFile backup_file, ae_exp_manager* exp_m)
{
  gzread(backup_file, &x_, sizeof(x_));
  gzread(backup_file, &y_, sizeof(y_));

  habitat_ = new Habitat(backup_file);
  
  individual_ = ae_individual::CreateIndividual(exp_m, backup_file);

  individual_->set_grid_cell(this);
  
  habitat_ = new Habitat();
}

// =================================================================
//                             Destructors
// =================================================================
ae_grid_cell::~ae_grid_cell(void)
{
  delete individual_;
  delete habitat_;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_grid_cell::save(gzFile backup_file) const
{
  gzwrite(backup_file, &x_, sizeof(x_));
  gzwrite(backup_file, &y_, sizeof(y_));

  habitat_->save(backup_file);

  individual_->save(backup_file);
}

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
