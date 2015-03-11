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


#ifndef __AE_WORLD_H__
#define __AE_WORLD_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_grid_cell.h"


namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
  class ae_exp_manager;






class World
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    World(void);
    World(const World&) = delete;

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~World(void);

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    // PRNGs
    ae_jumping_mt* get_prng(void) const;
    ae_jumping_mt* get_mut_prng(void) const;
    ae_jumping_mt* get_stoch_prng(void) const;

    std::list<ae_individual*>&& get_indivs_std(void) const;
    inline int32_t          get_nb_indivs(void) const;
    inline ae_individual*   get_best_indiv(void) const;
    inline int16_t          width()  const {return width_;};
    inline int16_t          height() const {return height_;};
    inline int32_t          get_migration_number(void) const;
    inline ae_grid_cell***  grid(void) const {return grid_;};
    inline ae_grid_cell*    get_grid_cell(int16_t x, int16_t y) const;
    inline ae_individual*   get_indiv_at(int16_t x, int16_t y) const;
  
    inline double** get_secretion_present_grid(void) const;
    inline double** get_secreted_amount_grid(void) const;
    inline double** get_metabolic_fitness_grid(void) const;
    inline double** get_total_fitness_grid(void) const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    // PRNGs
    inline void set_prng(ae_jumping_mt* prng);
    void set_mut_prng(ae_jumping_mt* prng);
    void set_stoch_prng(ae_jumping_mt* prng);

    inline void set_migration_number(int32_t migration_number);
    inline void set_secretion_degradation_prop(double degradation_prop);
    inline void set_secretion_diffusion_prop(double diffusion_prop);

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void InitGrid(int16_t width, int16_t height);
    void PlaceIndiv(ae_individual* indiv, int16_t x, int16_t y);
    void FillGridWithClones(ae_individual& dolly);
    void evaluate_individuals(Environment* envir);
    void update_secretion_grid(void); 
    void ShuffleIndivs(void);
    void update_best(void);
    void save(gzFile backup_file) const;
    void load(gzFile backup_file, ae_exp_manager* exp_man);

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================
    void MallocGrid(void);
    #ifndef DISTRIBUTED_PRNG
      void backup_stoch_prng(void);
    #endif

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_jumping_mt* _prng;

    #ifndef DISTRIBUTED_PRNG
      ae_jumping_mt* _mut_prng;
      ae_jumping_mt* _stoch_prng;
      ae_jumping_mt* _stoch_prng_bak;
    #endif
    
    int16_t width_; 
    int16_t height_;

    int16_t x_best, y_best;
    
    ae_grid_cell*** grid_;
    ae_grid_cell** grid_1d_;
    
    int32_t _migration_number;
    double  _secretion_diffusion_prop;
    double  _secretion_degradation_prop;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline int32_t World::get_nb_indivs(void) const
{
  return width_ * height_;
}

inline ae_individual* World::get_best_indiv(void) const
{
  return grid_[x_best][y_best]->get_individual();
}

inline int32_t World::get_migration_number(void) const
{
  return _migration_number;
}

inline ae_grid_cell* World::get_grid_cell(int16_t x,
                                                         int16_t y) const
{
  return grid_[x][y];
}

inline ae_individual* World::get_indiv_at(int16_t x, int16_t y) const
{
  return grid_[x][y]->get_individual();
}

inline double** World::get_secretion_present_grid(void) const
{
  double** ret = new double*[width_];
  
  for (int16_t x = 0; x < width_ ; x++)
  {
    ret[x] = new double[height_];
    for (int16_t y = 0; y < height_ ; y++)
    { 
      ret[x][y] = grid_[x][y]->get_compound_amount(); 
    }
  }
  
  return ret;
}

inline double** World::get_secreted_amount_grid(void) const
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

inline double** World::get_metabolic_fitness_grid(void) const
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

inline double** World::get_total_fitness_grid(void) const
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
inline void World::set_prng(ae_jumping_mt* prng)
{
  if (_prng != NULL) delete _prng;
  _prng = prng;
}

inline void World::set_migration_number(int32_t migration_number)
{
  _migration_number = migration_number;
}

inline void World::set_secretion_degradation_prop(double degradation_prop)
{
  _secretion_degradation_prop=degradation_prop;
}
inline void World::set_secretion_diffusion_prop(double diffusion_prop)
{
  _secretion_diffusion_prop=diffusion_prop;
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // __AE_WORLD_H__
