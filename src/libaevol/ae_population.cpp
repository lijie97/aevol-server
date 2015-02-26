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
#include <stdio.h>
#include <math.h>

#include <list>
#include <algorithm>

// =================================================================
//                            Project Files
// =================================================================
#include "ae_population.h"

#include "ae_exp_manager.h"
#include "ae_exp_setup.h"
#include "ae_individual.h"

#ifdef __NO_X
  #ifndef __REGUL

  #else
    #include "ae_individual_R.h"
  #endif
#elif defined __X11
  #ifndef __REGUL
    #include "ae_individual_X11.h"
  #else
    #include "ae_individual_R_X11.h"
  #endif
#endif

#include "ae_vis_a_vis.h"
#include "ae_align.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                             Class ae_population                             #
//                                                                             #
//##############################################################################



// =====================================================================
//                           Getters' definitions
// =====================================================================
std::list<ae_individual*> ae_population::get_indivs() const {
  return _indivs;
}

int32_t ae_population::get_nb_indivs( void ) const
{
  return _indivs.size();
}

ae_individual* ae_population::get_best( void ) const
{
  return _indivs.back();
}

/// Get the indiv corresponding to the given rank (1 for the worst indiv, POP_SIZE for the best)
///
/// Warning, be sure you call sort_individuals() before using get_indiv_by_rank
///
/// TODO vld: If this function is used often, suggests that _indivs
/// should be a vector instead of a list.
ae_individual* ae_population::get_indiv_by_rank( int32_t rank ) const
{
  auto indiv = _indivs.begin();
  std::advance(indiv, rank);
  assert((*indiv)->get_rank() == rank);
  return *indiv;
}

ae_jumping_mt* ae_population::get_mut_prng( void ) const
{
  return _mut_prng;
}

ae_jumping_mt* ae_population::get_stoch_prng( void ) const
{
  return _stoch_prng;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
void ae_population::add_indiv( ae_individual* indiv )
{
  _indivs.push_back(indiv);
}


void ae_population::set_mut_prng( ae_jumping_mt* prng )
{
  if (_mut_prng != NULL) delete _mut_prng;
  _mut_prng = prng;
  for (auto& indiv: _indivs)
    indiv->set_mut_prng(_mut_prng);
}

void ae_population::set_stoch_prng( ae_jumping_mt* prng )
{
  if (_stoch_prng != NULL)
    delete _stoch_prng;
  _stoch_prng = prng;

  for (auto& indiv: _indivs)
    indiv->set_stoch_prng(_stoch_prng);
}


// Mutation rates etc...
void ae_population::set_overall_point_mutation_rate( double point_mutation_rate )
{
  for (auto& indiv: _indivs)
    indiv->set_point_mutation_rate(point_mutation_rate);
}

void ae_population::set_overall_small_insertion_rate( double small_insertion_rate )
{
  for (auto& indiv: _indivs)
    indiv->set_small_insertion_rate(small_insertion_rate);
}

void ae_population::set_overall_small_deletion_rate( double small_deletion_rate )
{
  for (auto& indiv: _indivs)
    indiv->set_small_deletion_rate(small_deletion_rate);
}

void ae_population::set_overall_max_indel_size( int16_t max_indel_size )
{
  for (auto& indiv: _indivs)
    indiv->set_max_indel_size(max_indel_size);
}

void ae_population::set_overall_duplication_rate( double duplication_rate )
{
  for (auto& indiv: _indivs)
    indiv->set_duplication_rate(duplication_rate);
}

void ae_population::set_overall_deletion_rate( double deletion_rate)
{
  for (auto& indiv: _indivs)
    indiv->set_deletion_rate(deletion_rate);
}

void ae_population::set_overall_translocation_rate( double translocation_rate)
{
  for (auto& indiv: _indivs)
    indiv->set_translocation_rate(translocation_rate);
}

void ae_population::set_overall_inversion_rate( double inversion_rate)
{
  for (auto& indiv: _indivs)
    indiv->set_inversion_rate( inversion_rate );
}

void ae_population::set_overall_transfer_ins_rate (double transfer_ins_rate)
{
  for (auto& indiv: _indivs)
    indiv->set_HT_ins_rate( transfer_ins_rate );
}

void ae_population::set_overall_transfer_repl_rate (double transfer_repl_rate)
{
  for (auto& indiv: _indivs)
    indiv->set_HT_repl_rate( transfer_repl_rate );
}

void ae_population::set_overall_neighbourhood_rate( double neighbourhood_rate)
{
  for (auto& indiv: _indivs)
    indiv->set_neighbourhood_rate( neighbourhood_rate );
}

void ae_population::set_overall_duplication_proportion( double duplication_proportion)
{
  for (auto& indiv: _indivs)
    indiv->set_duplication_proportion( duplication_proportion );
}

void ae_population::set_overall_deletion_proportion( double deletion_proportion)
{
  for (auto& indiv: _indivs)
    indiv->set_deletion_proportion( deletion_proportion );
}

void ae_population::set_overall_translocation_proportion( double translocation_proportion)
{
  for (auto& indiv: _indivs)
    indiv->set_translocation_proportion( translocation_proportion );
}

void ae_population::set_overall_inversion_proportion( double inversion_proportion)
{
  for (auto& indiv: _indivs)
    indiv->set_inversion_proportion( inversion_proportion );
}

void ae_population::set_replication_reports( ae_tree* tree, int32_t num_gener)
{
  for (auto& indiv: _indivs)
    indiv->set_replication_report( tree->get_report_by_index( num_gener, indiv->get_id()));
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
void ae_population::evaluate_individuals( Environment* envir )
{
  for (auto& indiv: _indivs) {
    indiv->evaluate( envir );
    indiv->compute_statistical_data();
  }
}


// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_population::ae_population( ae_exp_manager* exp_m )
{
  _exp_m = exp_m;

  #ifndef DISTRIBUTED_PRNG
    _mut_prng       = NULL;
    _stoch_prng     = NULL;
    _stoch_prng_bak = NULL;
  #endif
}



// =================================================================
//                             Destructors
// =================================================================
ae_population::~ae_population( void )
{
  #ifndef DISTRIBUTED_PRNG
    delete _mut_prng;
    delete _stoch_prng;
    delete _stoch_prng_bak;
  #endif

    for (auto& indiv: _indivs)
      delete indiv;
    _indivs.clear();
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_population::set_nb_indivs(size_t nb_indivs)
{
  int32_t index_to_duplicate;
  ae_individual* indiv = NULL;
  if(nb_indivs > _indivs.size())
  {
    int32_t initial_pop_size = _indivs.size();
    for(size_t i = initial_pop_size; i < nb_indivs; i++)
    {
      index_to_duplicate = _exp_m->get_sel()->get_prng()->random( initial_pop_size );
      indiv = new ae_individual(*get_indiv_by_id(index_to_duplicate), true);
      indiv->set_id(i);
      add_indiv(indiv);
    }
  }
  else if(nb_indivs < _indivs.size())
  {
    std::list<ae_individual*> new_population;
    for(size_t i = 0; i < nb_indivs; i++)
    {
      index_to_duplicate = _exp_m->get_sel()->get_prng()->random(static_cast<int32_t>(_indivs.size()));
      indiv = new ae_individual(*get_indiv_by_id(index_to_duplicate), true);
      indiv->set_id(i);
      new_population.emplace_back(indiv);
    }
    update_population(std::move(new_population));
  }
  sort_individuals();
}

void ae_population::replace_population(std::list<ae_individual*>&& new_indivs)
{
  // First replace the former indivs by the new ones
  update_population(std::move(new_indivs));

  // Then reconciliate any possible inconsistency...

  // Update pointer to exp_manager in each individual
  // Replace indivs id by a new one
  int32_t id = 0;
  for (auto& indiv: _indivs) {
    indiv->set_exp_m(_exp_m);
    indiv->set_id(id++);
  }
}

void ae_population::save( gzFile backup_file ) const
{
  // Write population intrinsic data
  #ifndef DISTRIBUTED_PRNG
    _mut_prng->save( backup_file );
    int8_t tmp_with_stoch = _stoch_prng == NULL ? 0 : 1;
    gzwrite( backup_file, &tmp_with_stoch, sizeof(tmp_with_stoch) );
    if ( tmp_with_stoch )
    {
      _stoch_prng->save( backup_file );
    }
  #endif
  const int32_t nb_indivs = _indivs.size();
  gzwrite( backup_file, &nb_indivs, sizeof(nb_indivs) );

  // Write individuals
  for (const auto& indiv: _indivs)
    indiv->save( backup_file );
}

void ae_population::load( gzFile backup_file, bool verbose )
{
  // --------------------------------------- Retreive population intrinsic data
  #ifndef DISTRIBUTED_PRNG
    _mut_prng   = new ae_jumping_mt( backup_file );
    int8_t tmp_with_stoch;
    gzread( backup_file, &tmp_with_stoch, sizeof(tmp_with_stoch) );
    if ( tmp_with_stoch )
    {
      _stoch_prng = new ae_jumping_mt( backup_file );
    }
  #endif
  int32_t nb_indivs = _indivs.size();
  gzread( backup_file, &nb_indivs, sizeof(nb_indivs) );

  // ----------------------------------------------------- Retreive individuals
  if ( verbose ) printf( "  Loading individuals " );
  ae_individual* indiv = NULL;
  int32_t nb_ind_div_10 = nb_indivs / 10;
  for ( size_t i = 0 ; i < static_cast<size_t>(nb_indivs) ; i++ )
  {
    if ( verbose && i && i % nb_ind_div_10 == 0 )
    {
      putchar( '*' );
      fflush( stdout );
    }

    #ifdef __NO_X
      #ifndef __REGUL
        indiv = new ae_individual( _exp_m, backup_file );
      #else
        indiv = new ae_individual_R( _exp_m, backup_file );
      #endif
    #elif defined __X11
      #ifndef __REGUL
        indiv = new ae_individual_X11( _exp_m, backup_file );
      #else
        indiv = new ae_individual_R_X11( _exp_m, backup_file );
      #endif
    #endif

    _indivs.emplace_back( indiv );
  }
}

void ae_population::load(const char* backup_file_name, bool verbose)
{
  gzFile backup_file = gzopen(backup_file_name, "r");
  if ( backup_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open pop file %s\n",
            __FILE__, __LINE__, backup_file_name );
    exit( EXIT_FAILURE );
  }
  this->load(backup_file, verbose);
}

#ifndef DISTRIBUTED_PRNG
  void ae_population::backup_stoch_prng( void )
  {
    delete _stoch_prng_bak;
    _stoch_prng_bak = new ae_jumping_mt( *_stoch_prng );
  }
#endif


// =================================================================
//                           Protected Methods
// =================================================================
void ae_population::sort_individuals( void )
{
  _indivs.sort([](const ae_individual* i1,
                  const ae_individual* i2) { return i1->get_fitness() < i2->get_fitness(); });

  // Update the rank of the individuals
  int32_t rank = 1;
  for (auto& indiv: _indivs) {
    indiv->set_rank(rank);
    ++rank;
  }
}

/// Find the best individual and put it at the end of the list: this
/// is quicker than sorting the whole list in case we only need the
/// best individual, for example when we have spatial structure.
void ae_population::update_best( void )
{
  auto best = std::max_element(_indivs.begin(),
                               _indivs.end(),
                               [](const ae_individual* i1,
                                  const ae_individual* i2) -> bool {
                                 return i1->get_fitness() < i2->get_fitness();
                               });
  _indivs.splice(_indivs.end(), _indivs, best);
  (*best)->set_rank(_indivs.size());
}

// The new pop must be consistent and belong to the same experiment as the one
// it replaces, otherwise use replace_population
void ae_population::update_population(std::list<ae_individual*>&& new_indivs)
{
  // Delete the former indivs
  for (auto& indiv: _indivs)
    delete indiv;
  _indivs.clear();

  // Replace with new indivs
  _indivs = new_indivs;
}

ae_individual* ae_population::create_clone( ae_individual* dolly, int32_t id )
{
  ae_individual* indiv;

  #ifdef __NO_X
    #ifndef __REGUL
      indiv = new ae_individual( *dolly, false );
    #else
      indiv = new ae_individual_R( *(dynamic_cast<ae_individual_R*>(dolly)) , false);
    #endif
  #elif defined __X11
    #ifndef __REGUL
      indiv = new ae_individual_X11( *(dynamic_cast<ae_individual_X11*>(dolly)), false );
    #else
      indiv = new ae_individual_R_X11( *(dynamic_cast<ae_individual_R_X11*>(dolly)), false );
    #endif
  #endif

  indiv->set_id( id );

  return indiv;
}





// =================================================================
//                          Non inline accessors
// =================================================================
ae_individual* ae_population::get_indiv_by_id( int32_t id ) const
{
  for (const auto& indiv: _indivs)
    if (indiv->get_id() == id)
      return indiv;
  return nullptr;
}

} // namespace aevol
