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


#ifndef AEVOL_SELECTION_H_
#define AEVOL_SELECTION_H_


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>



// =================================================================
//                            Project Files
// =================================================================
#include "World.h"
#include "Observable.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;






class Selection : public Observable
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    Selection() = delete;
    Selection(const Selection&) = delete;
    Selection(ExpManager* exp_m);

    // =================================================================
    //                             Destructors
    // =================================================================

    virtual ~Selection();

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    inline SelectionScheme selection_scheme() const;
    inline double               selection_pressure() const;
    inline double*              prob_reprod() const;
    // inline std::unique_ptr<JumpingMT> prng() const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    // ----------------------------------------- Pseudo-random number generator
    inline void set_prng(std::unique_ptr<JumpingMT>&& prng);

    // -------------------------------------------------------------- Selection
    inline void set_selection_scheme(SelectionScheme sel_scheme);
    inline void set_selection_pressure(double sel_pressure);

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void step_to_next_generation();
    void PerformPlasmidTransfers();
    void write_setup_file(gzFile setup_file) const;
    void save(gzFile& backup_file) const;
    void load(gzFile& exp_setup_file, gzFile& backup_file, bool verbose);

    Individual* do_replication(Individual* parent, int32_t index = -1);
    Individual* do_replication(Individual* parent,
                               int32_t index,
                               int16_t x,
                               int16_t y);

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================
    void compute_prob_reprod();
    void compute_local_prob_reprod();
    Individual* do_local_competition(int16_t x, int16_t y);

    // =======================================================================
    //                             Protected Attributes
    // =======================================================================
    ExpManager* exp_m_;

    // ----------------------------------------- Pseudo-random number generator
    std::unique_ptr<JumpingMT> prng_;

    // -------------------------------------------------------------- Selection
    SelectionScheme selection_scheme_;
    double selection_pressure_;

    // --------------------------- Probability of reproduction of each organism
    double* prob_reprod_;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
// inline std::unique_ptr<JumpingMT> Selection::prng() const
// {
//   return prng_;
// }

inline SelectionScheme Selection::selection_scheme() const
{
  return selection_scheme_;
}

inline double Selection::selection_pressure() const
{
  return selection_pressure_;
}

inline double*Selection::prob_reprod() const
{
  if (prob_reprod_ == NULL)
  {
    printf("ERROR, prob_reprod_ has not been computed %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  return prob_reprod_;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
// ----------------------------------------- Pseudo-random number generator
inline void Selection::set_prng(std::unique_ptr<JumpingMT>&& prng)
{
  prng_ = std::move(prng);
}

// -------------------------------------------------------------- Selection
inline void Selection::set_selection_scheme(SelectionScheme sel_scheme)
{
  selection_scheme_ = sel_scheme;
}

inline void Selection::set_selection_pressure(double sel_pressure)
{
  selection_pressure_ = sel_pressure;
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_SELECTION_H_
