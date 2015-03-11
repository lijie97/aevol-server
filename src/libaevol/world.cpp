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

  _grid_width = -1;
  _grid_height = -1;

  _pop_grid = NULL;

  _migration_number = -1;
  _secretion_diffusion_prop = -1;
  _secretion_degradation_prop = -1;
}

// =================================================================
//                             Destructor
// =================================================================
World::~World(void)
{
  for (int16_t x = 0 ; x < _grid_width ; x++)
  {
    for (int16_t y = 0 ; y < _grid_height ; y++)
    {
      delete _pop_grid[x][y];
    }

    delete [] _pop_grid[x];
  }
  delete [] _pop_grid;

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
void World::place_indiv(ae_individual* indiv,
                                       int16_t x, int16_t y)
{
  _pop_grid[x][y]->set_individual(indiv);
}

void World::FillGridWithClones(ae_individual& dolly)
{
  int32_t id_new_indiv = 0;
  for (int16_t x = 0 ; x < _grid_width ; x++)
    for (int16_t y = 0 ; y < _grid_height ; y++)
  {
    place_indiv(ae_individual::create_clone(&dolly, id_new_indiv++), x, y);
  }
}

void World::evaluate_individuals(Environment* envir)
{
  for (int16_t x = 0 ; x < _grid_width ; x++)
    for (int16_t y = 0 ; y < _grid_height ; y++)
    {
      get_indiv_at(x, y)->evaluate(envir);
      get_indiv_at(x, y)->compute_statistical_data();
    }
}

void World::update_secretion_grid(void)
{
  int16_t cur_x, cur_y;

  double ** new_secretion = new double*[_grid_width];
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    new_secretion[x] = new double[_grid_height];
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      new_secretion[x][y] = _pop_grid[x][y]->get_compound_amount();
    }
  }

  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      // look at the entire neighborhood
      for ( int8_t i = -1 ; i < 2 ; i++ )
      {
        for ( int8_t j = -1 ; j < 2 ; j ++ )
        {
          cur_x = (x + i + _grid_width)  % _grid_width;
          cur_y = (y + j + _grid_height) % _grid_height;

          // add the diffusion from the neighboring cells
          new_secretion[x][y] += _pop_grid[cur_x][cur_y]->get_compound_amount() * _secretion_diffusion_prop;
        }
      }
    }
  }

  // substract what has diffused from each cell, and calculate the compound degradation
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      _pop_grid[x][y]->set_compound_amount( new_secretion[x][y] - 9 * _pop_grid[x][y]->get_compound_amount()
                                                                    * _secretion_diffusion_prop );
      _pop_grid[x][y]->set_compound_amount(   _pop_grid[x][y]->get_compound_amount()
                                            * (1 - _secretion_degradation_prop) );
    }
  }
  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    delete [] new_secretion[x];
  }
  delete [] new_secretion;
}

void World::do_random_migrations ( void )
{
  ae_individual * tmp_swap;

  int16_t old_x; int16_t old_y; int16_t new_x; int16_t new_y;
  for ( int16_t i = 0 ; i < _migration_number ; i++ )
  {
    old_x = (int16_t) (_prng->random() * _grid_width);
    old_y = (int16_t) (_prng->random() * _grid_height);
    new_x = (int16_t) (_prng->random() * _grid_width);
    new_y = (int16_t) (_prng->random() * _grid_height);


    // swap the individuals in these grid cells...
    tmp_swap = _pop_grid[old_x][old_y]->get_individual();
    _pop_grid[old_x][old_y]->set_individual( _pop_grid[new_x][new_y]->get_individual() );
    _pop_grid[new_x][new_y]->set_individual( tmp_swap );
  }
}

void World::update_best(void)
{
  x_best = y_best = 0;
  double fit_best = get_indiv_at(0, 0)->get_fitness();
  for (int16_t x = 0 ; x < _grid_width ; x++)
    for (int16_t y = 0 ; y < _grid_height ; y++)
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
  if ( _pop_grid == NULL )
  {
    printf( "%s:%d: error: grid not initialized.\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }

  gzwrite(backup_file, &_grid_width,   sizeof(_grid_width));
  gzwrite(backup_file, &_grid_height,  sizeof(_grid_height));

  for ( int16_t x = 0 ; x < _grid_width ; x++ )
  {
    for ( int16_t y = 0 ; y < _grid_height ; y++ )
    {
      _pop_grid[x][y]->save( backup_file );
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

  gzread(backup_file, &_grid_width,  sizeof(_grid_width));
  gzread(backup_file, &_grid_height, sizeof(_grid_height));

  _pop_grid = new ae_grid_cell** [_grid_width];
  for (int16_t x = 0 ; x < _grid_width ; x++)
  {
    _pop_grid[x] = new ae_grid_cell* [_grid_height];
    for (int16_t y = 0 ; y < _grid_height ; y++)
    {
      _pop_grid[x][y] = new ae_grid_cell(backup_file, exp_man);
    }
  }

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
    _stoch_prng_bak = new ae_jumping_mt( *_stoch_prng );
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

  for (int16_t x = 0 ; x < _grid_width ; x++)
    for (int16_t y = 0 ; y < _grid_height ; y++)
  {
    r.push_back(get_indiv_at(x, y));
  }

  return std::move(r);
}

void World::set_mut_prng(ae_jumping_mt* prng)
{
  if (_mut_prng != NULL) delete _mut_prng;
  _mut_prng = prng;

  for (int16_t x = 0 ; x < _grid_width ; x++)
    for (int16_t y = 0 ; y < _grid_height ; y++)
  {
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

  for (int16_t x = 0 ; x < _grid_width ; x++)
    for (int16_t y = 0 ; y < _grid_height ; y++)
  {
    ae_individual* indiv;
    if ((indiv = get_indiv_at(x, y)))
      indiv->set_stoch_prng(_stoch_prng);
  }
}
} // namespace aevol
