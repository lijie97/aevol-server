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
#include "SaveWorld.h"

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#ifdef __REGUL
#include "raevol/Individual_R.h"
#endif


namespace aevol {


//##############################################################################
//                                                                             #
//                                Class SaveWorld                              #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

// =================================================================
//                             Destructor
// =================================================================
SaveWorld::~SaveWorld()
{
  for (int16_t x = 0 ; x < width_ ; x++) {
    for (int16_t y = 0 ; y < height_ ; y++)
      delete grid_[x][y];
    delete[] grid_[x];
  }

  // grid_ is 2D accessible but 1D allocated, there were only 2 new
  // statements and these are the corresponding deletes
  delete [] grid_;
}

// =================================================================
//                            Public Methods
// =================================================================

void SaveWorld::MallocGrid()
{
  grid_ = new SaveGridCell ** [width_];
  for(int16_t x = 0 ; x < width_ ; x++)
    grid_[x] = new SaveGridCell * [height_];
}

void SaveWorld::save(gzFile backup_file)
{
  if (prng_ == nullptr)
  {
    printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  prng_->save(backup_file);

  mut_prng_->save(backup_file);

  int8_t tmp_with_stoch = static_cast<int8_t>(stoch_prng_ == nullptr ? 0 : 1);
  gzwrite(backup_file, &tmp_with_stoch, sizeof(tmp_with_stoch));
  if (tmp_with_stoch)
  {
    stoch_prng_->save(backup_file);
  }
  if (grid_ == nullptr)
  {
    printf("%s:%d: error: grid not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  // Manage shared or private phenotypic targets
  int8_t tmp_phenotypic_target_shared =
      static_cast<int8_t>(phenotypic_target_shared_ ? 1 : 0);
  gzwrite(backup_file,
          &tmp_phenotypic_target_shared,
          sizeof(tmp_phenotypic_target_shared));
  if (phenotypic_target_shared_) {
    phenotypic_target_handler_->save(backup_file);
  }


  gzwrite(backup_file, &width_,   sizeof(width_));
  gzwrite(backup_file, &height_,  sizeof(height_));

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      grid_[x][y]->save(backup_file, phenotypic_target_shared_);

  gzwrite(backup_file, &x_best, sizeof(x_best));
  gzwrite(backup_file, &y_best, sizeof(y_best));

  gzwrite(backup_file, &is_well_mixed_,
          sizeof(is_well_mixed_));
  gzwrite(backup_file, &partial_mix_nb_permutations_,
          sizeof(partial_mix_nb_permutations_));
  gzwrite(backup_file, &secretion_diffusion_prop_,
          sizeof(secretion_diffusion_prop_));
  gzwrite(backup_file, &secretion_degradation_prop_,
          sizeof(secretion_degradation_prop_));
}
} // namespace aevol
