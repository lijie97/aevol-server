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
#include <list>


// =================================================================
//                            Project Files
// =================================================================
#include "world.h"


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
World::World(void)
{
  _prng = NULL;

  #ifndef DISTRIBUTED_PRNG
    _mut_prng       = NULL;
    _stoch_prng     = NULL;
    _stoch_prng_bak = NULL;
  #endif

  width_ = -1;
  height_ = -1;

  grid_ = NULL;

  _migration_number = -1;
  _secretion_diffusion_prop = -1;
  _secretion_degradation_prop = -1;
}

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

  delete _prng;
  #ifndef DISTRIBUTED_PRNG
    delete _mut_prng;
    delete _stoch_prng;
    delete _stoch_prng_bak;
  #endif
}

// =================================================================
//                            Public Methods
// =================================================================
void World::InitGrid(int16_t width, int16_t height)
{
  width_  = width;
  height_ = height;

  MallocGrid();

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      grid_[x][y] = new ae_grid_cell(x, y, NULL);
}

void World::MallocGrid(void)
{
  // Although grid_ is a 2D array, we want all its cells to be contiguous
  // in memory. However, we also want it to be 2D-accessible i.e. we want to
  // be able to access a cell with grid_[x][y].
  // The following code does just this
  grid_1d_ = new ae_grid_cell* [width_ * height_];
  grid_ = new ae_grid_cell** [width_];
  for (int16_t x = 0 ; x < width_ ; x++)
    grid_[x] = &(grid_1d_[x * height_]);
}

void World::PlaceIndiv(ae_individual* indiv,
                                       int16_t x, int16_t y)
{
  grid_[x][y]->set_individual(indiv);
}

void World::FillGridWithClones(ae_individual& dolly)
{
  int32_t id_new_indiv = 0;
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      PlaceIndiv(ae_individual::create_clone(&dolly, id_new_indiv++), x, y);
}

void World::evaluate_individuals(Environment* envir)
{
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      get_indiv_at(x, y)->evaluate(envir);
      get_indiv_at(x, y)->compute_statistical_data();
    }
}

void World::update_secretion_grid(void)
{
  int16_t cur_x, cur_y;

  double ** new_secretion = new double*[width_];
  for ( int16_t x = 0 ; x < width_ ; x++ )
  {
    new_secretion[x] = new double[height_];
    for ( int16_t y = 0 ; y < height_ ; y++ )
    {
      new_secretion[x][y] = grid_[x][y]->get_compound_amount();
    }
  }

  for ( int16_t x = 0 ; x < width_ ; x++ )
  {
    for ( int16_t y = 0 ; y < height_ ; y++ )
    {
      // look at the entire neighborhood
      for ( int8_t i = -1 ; i < 2 ; i++ )
      {
        for ( int8_t j = -1 ; j < 2 ; j ++ )
        {
          cur_x = (x + i + width_)  % width_;
          cur_y = (y + j + height_) % height_;

          // add the diffusion from the neighboring cells
          new_secretion[x][y] += grid_[cur_x][cur_y]->get_compound_amount() * _secretion_diffusion_prop;
        }
      }
    }
  }

  // substract what has diffused from each cell, and calculate the compound degradation
  for ( int16_t x = 0 ; x < width_ ; x++ )
  {
    for ( int16_t y = 0 ; y < height_ ; y++ )
    {
      grid_[x][y]->set_compound_amount(new_secretion[x][y] -
          9 * grid_[x][y]->get_compound_amount() * _secretion_diffusion_prop);
      grid_[x][y]->set_compound_amount(grid_[x][y]->get_compound_amount() *
          (1 - _secretion_degradation_prop));
    }
  }
  for ( int16_t x = 0 ; x < width_ ; x++ )
  {
    delete [] new_secretion[x];
  }
  delete [] new_secretion;
}

/*
 * Suffle individuals randomly using Fisher-Yates shuffle
 */
void World::shuffle_indivs(void)
{
  // for (i from n − 1 downto 1) {
  //   j ← random integer with 0 ≤ j ≤ i
  //   exchange a[j] and a[i]
  // }


  ae_individual * tmp_swap;

  int16_t old_x; int16_t old_y; int16_t new_x; int16_t new_y;
  for (int16_t i = 0 ; i < _migration_number ; i++)
  {
    old_x = (int16_t) (_prng->random() * width_);
    old_y = (int16_t) (_prng->random() * height_);
    new_x = (int16_t) (_prng->random() * width_);
    new_y = (int16_t) (_prng->random() * height_);


    // swap the individuals in these grid cells...
    tmp_swap = grid_[old_x][old_y]->get_individual();
    grid_[old_x][old_y]->set_individual( grid_[new_x][new_y]->get_individual() );
    grid_[new_x][new_y]->set_individual( tmp_swap );
  }
}

void World::update_best(void)
{
  x_best = y_best = 0;
  double fit_best = get_indiv_at(0, 0)->get_fitness();
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
  {
    if (get_indiv_at(x, y)->get_fitness() > fit_best)
    {
      x_best = x;
      y_best = y;
      fit_best = get_indiv_at(x, y)->get_fitness();
    }
  }
}

void World::save(gzFile backup_file) const
{
  if (_prng == NULL)
  {
    printf( "%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }

  _prng->save(backup_file);

  #ifndef DISTRIBUTED_PRNG
    _mut_prng->save(backup_file);

    int8_t tmp_with_stoch = _stoch_prng == NULL ? 0 : 1;
    gzwrite( backup_file, &tmp_with_stoch, sizeof(tmp_with_stoch) );
    if (tmp_with_stoch)
    {
      _stoch_prng->save(backup_file);
    }
  #endif
  if ( grid_ == NULL )
  {
    printf( "%s:%d: error: grid not initialized.\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }

  gzwrite(backup_file, &width_,   sizeof(width_));
  gzwrite(backup_file, &height_,  sizeof(height_));

  for ( int16_t x = 0 ; x < width_ ; x++ )
  {
    for ( int16_t y = 0 ; y < height_ ; y++ )
    {
      grid_[x][y]->save( backup_file );
    }
  }

  gzwrite(backup_file, &x_best, sizeof(x_best));
  gzwrite(backup_file, &y_best, sizeof(y_best));

  gzwrite(backup_file, &_migration_number,           sizeof(_migration_number));
  gzwrite(backup_file, &_secretion_diffusion_prop,   sizeof(_secretion_diffusion_prop));
  gzwrite(backup_file, &_secretion_degradation_prop, sizeof(_secretion_degradation_prop));
}

void World::load(gzFile backup_file, ae_exp_manager* exp_man)
{
  _prng = new ae_jumping_mt(backup_file);
  #ifndef DISTRIBUTED_PRNG
    _mut_prng   = new ae_jumping_mt(backup_file);
    int8_t tmp_with_stoch;
    gzread( backup_file, &tmp_with_stoch, sizeof(tmp_with_stoch) );
    if (tmp_with_stoch)
    {
      _stoch_prng = new ae_jumping_mt(backup_file);
    }
  #endif

  gzread(backup_file, &width_,  sizeof(width_));
  gzread(backup_file, &height_, sizeof(height_));

  MallocGrid();

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      grid_[x][y] = new ae_grid_cell(backup_file, exp_man);

  gzread(backup_file, &x_best, sizeof(x_best));
  gzread(backup_file, &y_best, sizeof(y_best));

  gzread(backup_file, &_migration_number,           sizeof(_migration_number));
  gzread(backup_file, &_secretion_diffusion_prop,   sizeof(_secretion_diffusion_prop));
  gzread(backup_file, &_secretion_degradation_prop, sizeof(_secretion_degradation_prop));
}

// =================================================================
//                           Protected Methods
// =================================================================
#ifndef DISTRIBUTED_PRNG
  void World::backup_stoch_prng( void )
  {
    delete _stoch_prng_bak;
    _stoch_prng_bak = new ae_jumping_mt(*_stoch_prng);
  }
#endif

// =================================================================
//                          Non inline accessors
// =================================================================
ae_jumping_mt* World::get_prng(void) const
{
  return _prng;
}

ae_jumping_mt* World::get_mut_prng( void ) const
{
  return _mut_prng;
}

ae_jumping_mt* World::get_stoch_prng( void ) const
{
  return _stoch_prng;
}

list<ae_individual*>&& World::get_indivs_std(void) const
{
  list<ae_individual*> r;

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      r.push_back(get_indiv_at(x, y));

  return std::move(r);
}

void World::set_mut_prng(ae_jumping_mt* prng)
{
  if (_mut_prng != NULL) delete _mut_prng;
  _mut_prng = prng;

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      ae_individual* indiv;
      if ((indiv = get_indiv_at(x, y)))
        indiv->set_mut_prng(_mut_prng);
    }
}

void World::set_stoch_prng(ae_jumping_mt* prng)
{
  if (_stoch_prng != NULL)
    delete _stoch_prng;
  _stoch_prng = prng;

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      ae_individual* indiv;
      if ((indiv = get_indiv_at(x, y)))
        indiv->set_stoch_prng(_stoch_prng);
    }
}
} // namespace aevol
