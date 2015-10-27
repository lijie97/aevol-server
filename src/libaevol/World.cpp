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
#include "World.h"

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include <list>
#include <iostream>

using std::cout;
using std::endl;
using std::list;


namespace aevol {


//##############################################################################
//                                                                             #
//                                Class World                                  #
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
World::~World(void)
{
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      delete grid_[x][y];
  
  // grid_ is 2D accessible but 1D allocated, there were only 2 new
  // statements and these are the corresponding deletes
  delete [] grid_1d_;
  delete [] grid_;
}

// =================================================================
//                            Public Methods
// =================================================================
void World::InitGrid(int16_t width, int16_t height,
                     const Habitat& habitat,
                     bool share_phenotypic_target)
{
  assert(share_phenotypic_target);

  if (share_phenotypic_target)
    phenotypic_target_handler_ = std::make_shared<PhenotypicTargetHandler>
        (habitat.phenotypic_target_handler());

  width_  = width;
  height_ = height;

  MallocGrid();

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
    {
      if (share_phenotypic_target)
        grid_[x][y] =
            new GridCell(x, y,
#if __cplusplus == 201103L
                             make_unique<Habitat>
                                 (habitat, share_phenotypic_target),
#else
                             std::make_unique<Habitat>
                                 (habitat, share_phenotypic_target),
#endif
                             NULL);
    }
}

void World::MallocGrid(void)
{
  // Although grid_ is a 2D array, we want all its cells to be contiguous
  // in memory. However, we also want it to be 2D-accessible i.e. we want to
  // be able to access a cell with grid_[x][y].
  // The following code does just this
  grid_1d_ = new GridCell * [width_ * height_];
  grid_ = new GridCell ** [width_];
  for (int16_t x = 0 ; x < width_ ; x++)
    grid_[x] = &(grid_1d_[x * height_]);
}

void World::PlaceIndiv(Individual * indiv, int16_t x, int16_t y) {
  grid_[x][y]->set_individual(indiv);
}

void World::FillGridWithClones(Individual & dolly)
{
  int32_t id_new_indiv = 0;
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      PlaceIndiv(Individual::CreateClone(&dolly, id_new_indiv++), x, y);
}

void World::evaluate_individuals()
{
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      get_indiv_at(x, y)->Evaluate();
      get_indiv_at(x, y)->compute_statistical_data();
    }
}

void World::update_secretion_grid(void)
{
  int16_t cur_x, cur_y;

  double ** new_secretion = new double*[width_];
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    new_secretion[x] = new double[height_];
    for (int16_t y = 0 ; y < height_ ; y++)
    {
      new_secretion[x][y] = grid_[x][y]->compound_amount();
    }
  }

  for (int16_t x = 0 ; x < width_ ; x++)
  {
    for (int16_t y = 0 ; y < height_ ; y++)
    {
      // look at the entire neighborhood
      for (int8_t i = -1 ; i < 2 ; i++)
      {
        for (int8_t j = -1 ; j < 2 ; j ++)
        {
          cur_x = (x + i + width_)  % width_;
          cur_y = (y + j + height_) % height_;

          // add the diffusion from the neighboring cells
          new_secretion[x][y] +=
              grid_[cur_x][cur_y]->compound_amount() *
                  secretion_diffusion_prop_;
        }
      }
    }
  }

  // Substract what has diffused from each cell, and calculate the
  // compound degradation
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    for (int16_t y = 0 ; y < height_ ; y++)
    {
      grid_[x][y]->set_compound_amount(new_secretion[x][y] -
          9 * grid_[x][y]->compound_amount() * secretion_diffusion_prop_);
      grid_[x][y]->set_compound_amount(grid_[x][y]->compound_amount() *
          (1 - secretion_degradation_prop_));
    }
  }
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    delete [] new_secretion[x];
  }
  delete [] new_secretion;
}

/*
 * Perform mixing of individuals
 *
 * Depending on parameters this can either well-mix or partially-mix
 * the population
 */
void World::MixIndivs(void)
{
  if (is_well_mixed_)
    WellMixIndivs();
  else if (partial_mix_nb_permutations_ > 0)
    PartiallyMixIndivs();
}

/*
 * Suffle individuals randomly using Fisher-Yates shuffle
 */
void World::WellMixIndivs(void)
{
  for (int16_t i = width_ * height_ - 1 ; i > 0 ; i--) {
    int16_t j = prng_->random(i + 1); // random in [0, 1]

    // Swap individuals btw cells i and j
    Individual * tmp = grid_1d_[i]->get_individual();
    grid_1d_[i]->set_individual(grid_1d_[j]->get_individual());
    grid_1d_[j]->set_individual(tmp);
  }
}

/*
 * Perform permutations between individuals randomly
 *
 * The number of permutations is given by partial_mix_nb_permutations_
 */
void World::PartiallyMixIndivs(void)
{
  for (int32_t i = 0 ; i < partial_mix_nb_permutations_ ; i++)
  {
    int16_t old_x = prng_->random(width_);
    int16_t old_y = prng_->random(height_);
    int16_t new_x = prng_->random(width_);
    int16_t new_y = prng_->random(height_);

    // Swap the individuals in these grid cells...
    Individual * tmp_swap = grid_[old_x][old_y]->get_individual();
    grid_[old_x][old_y]->set_individual(grid_[new_x][new_y]->get_individual());
    grid_[new_x][new_y]->set_individual(tmp_swap);
  }
}

void World::update_best(void)
{
  x_best = y_best = 0;
  double fit_best = get_indiv_at(0, 0)->get_fitness();
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      if (get_indiv_at(x, y)->get_fitness() > fit_best) {
        x_best = x;
        y_best = y;
        fit_best = get_indiv_at(x, y)->get_fitness();
      }
}

void World::ApplyHabitatVariation() {
  if (phenotypic_target_shared_)
    phenotypic_target_handler_->ApplyVariation();
  else
    for (int16_t x = 0 ; x < width_ ; x++)
      for (int16_t y = 0 ; y < height_ ; y++)
        grid_[x][y]->ApplyHabitatVariation();
}

void World::save(gzFile backup_file) const
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
  if (phenotypic_target_shared_)
    phenotypic_target_handler_->save(backup_file);

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

void World::load(gzFile backup_file, ExpManager * exp_man)
{
  // Retrieve PRNGs
  prng_ = std::make_shared<JumpingMT>(backup_file);
  mut_prng_ = std::make_shared<JumpingMT>(backup_file);

  int8_t tmp_with_stoch;
  gzread(backup_file, &tmp_with_stoch, sizeof(tmp_with_stoch));
  if (tmp_with_stoch)
  {
    stoch_prng_ = std::make_shared<JumpingMT>(backup_file);
  }

  // Manage shared or private phenotypic targets
  int8_t tmp_phenotypic_target_shared;
  gzread(backup_file,
          &tmp_phenotypic_target_shared,
          sizeof(tmp_phenotypic_target_shared));
  phenotypic_target_shared_ = tmp_phenotypic_target_shared;
  if (phenotypic_target_shared_)
    phenotypic_target_handler_ =
        std::make_shared<PhenotypicTargetHandler>(backup_file);
  phenotypic_target_handler_->BuildPhenotypicTarget();

  gzread(backup_file, &width_,  sizeof(width_));
  gzread(backup_file, &height_, sizeof(height_));

  MallocGrid();

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      grid_[x][y] = new GridCell(backup_file,
                                 exp_man,
                                 phenotypic_target_handler_);

  gzread(backup_file, &x_best, sizeof(x_best));
  gzread(backup_file, &y_best, sizeof(y_best));

  gzread(backup_file, &is_well_mixed_,
         sizeof(is_well_mixed_));
  gzread(backup_file, &partial_mix_nb_permutations_,
         sizeof(partial_mix_nb_permutations_));
  gzread(backup_file, &secretion_diffusion_prop_,
         sizeof(secretion_diffusion_prop_));
  gzread(backup_file, &secretion_degradation_prop_,
         sizeof(secretion_degradation_prop_));
}

// =================================================================
//                           Protected Methods
// =================================================================
void World::backup_stoch_prng(void)
{
  // Store a copy of stoch_prng_ in stoch_prng_bak_
#if __cplusplus == 201103L
  stoch_prng_bak_ = make_unique<JumpingMT>(*stoch_prng_);
#else
  stoch_prng_bak_ = std::make_unique<JumpingMT>(*stoch_prng_);
#endif
}

// =================================================================
//                          Non inline accessors
// =================================================================
std::shared_ptr<JumpingMT> World::get_prng(void) const
{
  return prng_;
}

std::shared_ptr<JumpingMT> World::get_mut_prng(void) const
{
  return mut_prng_;
}

std::shared_ptr<JumpingMT> World::get_stoch_prng(void) const
{
  return stoch_prng_;
}

list<Individual *> World::get_indivs(void) const
{
  list<Individual *> r;

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      r.push_back(get_indiv_at(x, y));

  return r;
}

void World::set_mut_prng(std::shared_ptr<JumpingMT> prng)
{
  mut_prng_ = prng;

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      Individual * indiv;
      if ((indiv = get_indiv_at(x, y)))
        indiv->set_mut_prng(mut_prng_);
    }
}

void World::set_stoch_prng(std::shared_ptr<JumpingMT> prng)
{
  stoch_prng_ = prng;

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      Individual * indiv;
      if ((indiv = get_indiv_at(x, y)))
        indiv->set_stoch_prng(stoch_prng_);
    }
}

void World::set_phen_target_prngs(std::shared_ptr<JumpingMT> var_prng,
                                  std::shared_ptr<JumpingMT> noise_prng) {
  assert(phenotypic_target_shared_);
  phenotypic_target_handler_->set_var_prng(var_prng);
  phenotypic_target_handler_->set_noise_prng(noise_prng);
}

Individual* World::get_indiv_by_id(int32_t id) const {
  Individual* indiv = grid_1d_[id]->get_individual();
  // When the population isn't mixed at all, the individual with id n is in
  // grid_1d_[n]. Try this first...
  if (indiv->get_id() == id)
    return indiv;
  // ... If it isn't, do a basic search
  int32_t nb_indivs = width_ * height_;
  for (int32_t i = 0 ; i < nb_indivs ; i++) {
    if (grid_1d_[i]->get_individual()->get_id() == id)
      return grid_1d_[i]->get_individual();
  }
  return nullptr;
}
} // namespace aevol
