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


// =================================================================
//                            Project Files
// =================================================================
#include <ae_population.h>

#include <ae_exp_manager.h>
#include <ae_exp_setup.h>
#include <ae_individual.h>

#ifdef __NO_X
  #ifndef __REGUL

  #else
    #include <ae_individual_R.h>
  #endif
#elif defined __X11
  #ifndef __REGUL
    #include <ae_individual_X11.h>
  #else
    #include <ae_individual_R_X11.h>
  #endif
#endif

#include <ae_vis_a_vis.h>
#include <ae_align.h>

//##############################################################################
//                                                                             #
//                             Class ae_population                             #
//                                                                             #
//##############################################################################

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

  // Individuals
  _nb_indivs  = 0;
  _indivs     = new ae_list<ae_individual*>();
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

  _indivs->erase( true );
  delete _indivs;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_population::set_nb_indivs(int32_t nb_indivs)
{
	int32_t index_to_duplicate;
	ae_individual* indiv = NULL;
	if(nb_indivs > _nb_indivs)
	{
		int32_t initial_pop_size = _nb_indivs;
		for(int32_t i = initial_pop_size; i < nb_indivs; i++)
		{
			index_to_duplicate = _exp_m->get_sel()->get_prng()->random( initial_pop_size );
			indiv = new ae_individual(*get_indiv_by_id(index_to_duplicate), true);
			indiv->set_id(i);
			add_indiv(indiv);
		}
	}
	else if(nb_indivs < _nb_indivs)
	{
		ae_list<ae_individual*>* new_population = new ae_list<ae_individual*>();
		for(int32_t i = 0; i < nb_indivs; i++)
		{
			index_to_duplicate = _exp_m->get_sel()->get_prng()->random( _nb_indivs );
			indiv = new ae_individual(*get_indiv_by_id(index_to_duplicate), true);
			indiv->set_id(i);
			new_population->add(indiv);
		}
		replace_population(new_population);
	}
	sort_individuals();
}

void ae_population::replace_population( ae_list<ae_individual*>* new_indivs )
{
  _indivs->erase( true );
  delete _indivs;

  _indivs = new_indivs;
  _nb_indivs = _indivs->get_nb_elts();
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
  gzwrite( backup_file, &_nb_indivs, sizeof(_nb_indivs) );

  // Write individuals
  ae_list_node<ae_individual*>*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->save( backup_file );
    indiv_node = indiv_node->get_next();
  }
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
  gzread( backup_file, &_nb_indivs, sizeof(_nb_indivs) );

  // ----------------------------------------------------- Retreive individuals
  if ( verbose ) printf( "  Loading individuals " );
  ae_individual* indiv = NULL;
  int32_t nb_ind_div_10 = _nb_indivs / 10;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
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

    _indivs->add( indiv );
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
  // Insertion sort
  ae_list_node<ae_individual*>* last_sorted   = _indivs->get_first();
  ae_list_node<ae_individual*>* comp          = NULL;
  ae_list_node<ae_individual*>* item_to_sort  = NULL;
  double  fit_comp;
  double  fitness_to_sort;
  bool    is_sorted;

  // only "pop_size - 1" iterations since the first item is already "sorted"
  for ( int32_t nb_sorted = 1 ; nb_sorted < _nb_indivs ; nb_sorted++ )
  {
    item_to_sort    = last_sorted->get_next();
    fitness_to_sort = item_to_sort->get_obj()->get_fitness();
    is_sorted       = false;

    comp = last_sorted;

    // looking for its place among the sorted items
    while ( !is_sorted )
    {
      fit_comp = comp->get_obj()->get_fitness();

      if ( fitness_to_sort >= fit_comp ) // The right place of the item is after comp
      {
        if ( item_to_sort->get_prev() == comp )
        {
          // item_to_sort is already at the right place
          is_sorted = true;
          last_sorted = item_to_sort;
        }

        // item_to_sort has to be inserted just after comp
        _indivs->remove( item_to_sort, false, false );
        _indivs->insert_after( item_to_sort, comp );
        is_sorted = true;
      }
      else
      {
        // move on to compare with the next item
        comp = comp->get_prev();
      }

      if ( comp == NULL )
      {
        // item_to_sort has to be inserted at the beginning of the list
        _indivs->remove( item_to_sort, false, false );
        _indivs->add_front( item_to_sort );
        is_sorted = true;
      }
    }
  }

  // Update the rank of the individuals
  ae_list_node<ae_individual*>* indiv_node  = _indivs->get_first();
  ae_individual* indiv       = NULL;

  for ( int32_t rank = 1 ; rank <= _nb_indivs ; rank++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_rank( rank );
    indiv_node = indiv_node->get_next();
  }

}

// Find the best individual and put it at the end of the list: this is quicker than sorting the whole list in case we only need the best individual, for example when we have spatial structure.
void ae_population::update_best( void )
{
  ae_list_node<ae_individual*>* current_best  = _indivs->get_first();
  ae_list_node<ae_individual*>* candidate     = _indivs->get_first();

  while ( candidate != NULL )
  {
    if (  candidate->get_obj()->get_fitness() >= current_best->get_obj()->get_fitness() )
    {
      current_best = candidate;
    }
    candidate = candidate->get_next();
  }

  _indivs->remove( current_best, false, false );
  _indivs->add( current_best );

  current_best->get_obj()->set_rank( _nb_indivs );
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
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  while ( indiv_node != NULL )
  {
    indiv = indiv_node->get_obj();

    if ( indiv->get_id() == id )
    {
      return indiv;
    }

    indiv_node = indiv_node->get_next();
  }

  return NULL;
}
