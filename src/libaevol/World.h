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


#ifndef AEVOL_WORLD_H__
#define AEVOL_WORLD_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>

#include <zlib.h>

#include "GridCell.h"


namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;






class World
{
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  World(void) = default;
  World(const World &) = delete;

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~World(void);

  // =================================================================
  //                        Accessors: getters
  // =================================================================
  // PRNGs
  std::shared_ptr<JumpingMT> get_prng(void) const;
  std::shared_ptr<JumpingMT> get_mut_prng(void) const;
  std::shared_ptr<JumpingMT> get_stoch_prng(void) const;

  std::list<Individual *> get_indivs(void) const;
  inline int32_t          get_nb_indivs(void) const;
  inline Individual *   get_best_indiv(void) const;
  int16_t          width()  const {return width_;};
  int16_t          height() const {return height_;};
  inline int32_t          partial_mix_nb_permutations(void) const;
  GridCell ***  grid(void) const {return grid_;};
  inline GridCell*    grid(int16_t x, int16_t y) const;
  inline Individual*   get_indiv_at(int16_t x, int16_t y) const;
  Individual* get_indiv_by_id(int32_t id) const;

  inline double** get_secretion_present_grid(void) const;
  inline double** get_secreted_amount_grid(void) const;
  inline double** get_metabolic_fitness_grid(void) const;
  inline double** get_total_fitness_grid(void) const;

  bool phenotypic_target_shared() const {
    return phenotypic_target_shared_;
  }
  const std::shared_ptr<PhenotypicTargetHandler>&
  phenotypic_target_handler() const {
    return phenotypic_target_handler_;
  }

  // =================================================================
  //                        Accessors: setters
  // =================================================================
  // PRNGs
  inline void set_prng(std::shared_ptr<JumpingMT> prng);
  void set_mut_prng(std::shared_ptr<JumpingMT> prng);
  void set_stoch_prng(std::shared_ptr<JumpingMT> prng);

  inline void set_is_well_mixed(bool is_well_mixed);
  inline void set_partial_mix_nb_permutations(int32_t nb_permutations);
  inline void set_secretion_degradation_prop(double degradation_prop);
  inline void set_secretion_diffusion_prop(double diffusion_prop);
  inline void set_best(int16_t x, int16_t y);

  // =================================================================
  //                              Operators
  // =================================================================

  // =================================================================
  //                            Public Methods
  // =================================================================
  void InitGrid(int16_t width, int16_t height,
                const Habitat& habitat,
                bool share_phenotypic_target);
  void PlaceIndiv(Individual * indiv, int16_t x, int16_t y);
  void FillGridWithClones(Individual & dolly);
  void evaluate_individuals();
  void update_secretion_grid(void);
  void MixIndivs(void);
  void update_best(void);
  void ApplyHabitatVariation();

  void save(gzFile backup_file) const;
  void load(gzFile backup_file, ExpManager * exp_man);

  void set_phen_target_prngs(std::shared_ptr<JumpingMT> var_prng,
                             std::shared_ptr<JumpingMT> noise_prng);

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  void MallocGrid(void);
  void WellMixIndivs(void);
  void PartiallyMixIndivs(void);
  void backup_stoch_prng(void);

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  std::shared_ptr<JumpingMT> _prng = nullptr;

  std::shared_ptr<JumpingMT> _mut_prng = nullptr;
  std::shared_ptr<JumpingMT> _stoch_prng = nullptr;
  std::unique_ptr<JumpingMT> _stoch_prng_bak = nullptr;

  int16_t width_  = -1;
  int16_t height_ = -1;

  int16_t x_best = -1;
  int16_t y_best = -1;

  GridCell*** grid_ = nullptr;
  GridCell** grid_1d_ = nullptr;

  bool is_well_mixed_ = false;
  int32_t partial_mix_nb_permutations_ = 0;

  bool phenotypic_target_shared_ = true;
  std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler_ = nullptr;

  double  _secretion_diffusion_prop = -1;
  double  _secretion_degradation_prop = -1;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline int32_t World::get_nb_indivs(void) const
{
  return width_ * height_;
}

inline Individual *World::get_best_indiv(void) const
{
  return grid_[x_best][y_best]->get_individual();
}

inline int32_t World::partial_mix_nb_permutations(void) const
{
  return partial_mix_nb_permutations_;
}

inline GridCell *World::grid(int16_t x, int16_t y) const
{
  return grid_[x][y];
}

inline Individual *World::get_indiv_at(int16_t x, int16_t y) const
{
  return grid_[x][y]->get_individual();
}

inline double**World::get_secretion_present_grid(void) const
{
  double** ret = new double*[width_];
  
  for (int16_t x = 0; x < width_ ; x++)
  {
    ret[x] = new double[height_];
    for (int16_t y = 0; y < height_ ; y++)
    { 
      ret[x][y] = grid_[x][y]->compound_amount(); 
    }
  }
  
  return ret;
}

inline double**World::get_secreted_amount_grid(void) const
{
  double** ret = new double*[width_];
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    ret[x] = new double[height_];
    for (int16_t y = 0; y < height_ ; y++)
    {
      ret[x][y] = grid_[x][y]->get_secreted_amount();
    }
  }
  
  return ret;
}

inline double**World::get_metabolic_fitness_grid(void) const
{
  double** ret = new double*[width_];
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    ret[x] = new double[height_];
    for (int16_t y = 0; y < height_ ; y++)
    {
      ret[x][y] = grid_[x][y]->get_metabolic_fitness();
    }
  }
  
  return ret;
}

inline double**World::get_total_fitness_grid(void) const
{
  double** ret = new double*[width_];
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    ret[x] = new double[height_];
    for (int16_t y = 0; y < height_ ; y++)
    {
      ret[x][y] = grid_[x][y]->get_total_fitness();
    }
  }
  
  return ret;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void World::set_prng(std::shared_ptr<JumpingMT> prng) {
  _prng = prng;
}

inline void World::set_is_well_mixed(bool is_well_mixed) {
  is_well_mixed_ = is_well_mixed;
}

inline void World::set_partial_mix_nb_permutations(int32_t nb_permutations) {
  partial_mix_nb_permutations_ = nb_permutations;
}

inline void World::set_secretion_degradation_prop(double degradation_prop) {
  _secretion_degradation_prop=degradation_prop;
}

inline void World::set_secretion_diffusion_prop(double diffusion_prop) {
  _secretion_diffusion_prop=diffusion_prop;
}

inline void World::set_best(int16_t x, int16_t y) {
  x_best = x;
  y_best = y;
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_WORLD_H__
