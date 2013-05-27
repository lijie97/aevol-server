//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, Virginie Lefort
//                    David Parsons
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************


/*! \class ae_selection
    \brief
*/
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <math.h>


// =================================================================
//                            Project Files
// =================================================================
#include <ae_selection.h>
#include <ae_exp_manager.h>
#include <ae_vis_a_vis.h>

#ifdef __NO_X
  #ifndef __REGUL
    #include <ae_individual.h>
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




//##############################################################################
//                                                                             #
//                              Class ae_selection                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_selection::ae_selection( ae_exp_manager* exp_m )
{
  _exp_m = exp_m;

  // ----------------------------------------- Pseudo-random number generator
  _prng = NULL;

  // -------------------------------------------------------------- Selection
  _selection_scheme   = RANK_EXPONENTIAL;
  _selection_pressure = 0.998;

  // --------------------------- Probability of reproduction of each organism
  _prob_reprod = NULL;
}

// =================================================================
//                             Destructors
// =================================================================
ae_selection::~ae_selection( void )
{
  delete _prng;
  delete [] _prob_reprod;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_selection::step_to_next_generation( void )
{
  // To create the new generation, we must create nb_indivs new individuals
  // (offsprings) and "kill" the existing ones.
  // The number of offspring on a given individual will be given by a stochastic
  // process biased on it's fitness value (the selection process).
  // There are 3 possible selection schemes :
  //    * Linear Ranking
  //    * Exponential Ranking
  //    * Fitness proportionate
  //
  // Whichever method is chosen, we will
  // 1) Compute the probability of reproduction of each individual in the population
  // 2) Simulate the stochastic process by a multinomial drawing (based upon the probabilities computed in 1)
  // 3) Make the selected individuals reproduce, thus creating the new generation
  // 4) Replace the current generation by the newly created one.
  // 5) Sort the newly created population*
  
  if ( _prng == NULL )
  {
    printf( "%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  if ( _exp_m->is_spatially_structured() )
  {
    step_to_next_generation_grid();
    return;
  }
  
  
  if ( _exp_m->get_pop()->get_indiv_by_id( 0 )->get_with_stochasticity() )
  {
    _exp_m->get_pop()->backup_stoch_prng();
  }


  // -------------------------------------------------------------------------------
  // 1) Compute the probability of reproduction of each individual in the population
  // -------------------------------------------------------------------------------
  #ifndef FIXED_POPULATION_SIZE
    #error this method is not ready for variable population size
    compute_prob_reprod();
  #else
    if ( _selection_scheme == FITNESS_PROPORTIONATE || _prob_reprod == NULL )
    {
      compute_prob_reprod();
    }
  #endif

  // --------------------------------------------------------------------------------------------------------
  // 2) Simulate the stochastic process by a multinomial drawing (based upon the probabilities computed in 1)
  // --------------------------------------------------------------------------------------------------------
  int32_t  nb_indivs = _exp_m->get_pop()->get_nb_indivs();
  int32_t* nb_offsprings = new int32_t[nb_indivs];
  _prng->multinomial_drawing( nb_offsprings, _prob_reprod, nb_indivs, nb_indivs );

  // ------------------------------------------------------------------------------
  // 3) Make the selected individuals "reproduce", thus creating the new generation
  // ------------------------------------------------------------------------------
  ae_list<ae_individual*>*      new_generation  = new ae_list<ae_individual*>();
  ae_list<ae_individual*>*      old_generation  = _exp_m->get_indivs();
  ae_list_node<ae_individual*>* indiv_node      = old_generation->get_first();
  ae_list_node<ae_individual*>* next_indiv_node = NULL;
  ae_individual*  indiv           = NULL;
  int32_t         index_new_indiv = 0;

  for ( int32_t i = 0 ; i < nb_indivs ; i++ )
  {
    // Make indiv i reproduce (nb_offsprings[i] offsprings)
    indiv = indiv_node->get_obj();
    
    next_indiv_node = indiv_node->get_next();
    
    for ( int32_t j = 0 ; j < nb_offsprings[i] ; j++ )
    {
      #ifdef DISTRIBUTED_PRNG
        #error Not implemented yet !
        // For all but the first time, Take a jump in the PRNG
        if ( j > 0 ) indiv->do_prng_jump();
      #endif
      
      // Create a new individual (evaluated at the end of do_replication)
      new_generation->add( do_replication( indiv, index_new_indiv++ ) );
    }
    
    // All the offsprings of this individual have been generated, if there is no transfer,
    // the indiv will not be used any more and can hence be deleted
    if ( (not _exp_m->get_with_HT()) and (not _exp_m->get_with_plasmids()) )
    {
      old_generation->remove( indiv_node, true, true );
    }
        
    indiv_node = next_indiv_node;
  }
  
  if ( _exp_m->get_with_HT() or _exp_m->get_with_plasmids() )
  {
    // The individuals have not yet been deleted, do it now.
    old_generation->erase( true );
  }
  

  delete [] nb_offsprings;



  // -------------------------------------------------------------
  //  4) Replace the current generation by the newly created one.
  // -------------------------------------------------------------
  assert( _exp_m->get_pop()->get_indivs()->is_empty() );
  _exp_m->get_pop()->replace_population( new_generation );


  // --------------------------------------
  //  5) Sort the newly created population
  // --------------------------------------
  _exp_m->get_pop()->sort_individuals();
}

void ae_selection::step_to_next_generation_grid( void )
{
  // -------------------------------------------------------------------------------
  // 1) Compute the probability of reproduction of each individual in the population
  // -------------------------------------------------------------------------------
  #ifndef FIXED_POPULATION_SIZE
    #error this method is not ready for variable population size
    compute_local_prob_reprod();
  #else
    if ( _selection_scheme == FITNESS_PROPORTIONATE || _prob_reprod == NULL )
    {
      compute_local_prob_reprod();
    }
  #endif
  
  if ( _prng == NULL )
  {
    printf( "%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  ae_spatial_structure* sp_struct = _exp_m->get_spatial_structure();
  
  int16_t grid_width  = _exp_m->get_grid_width();
  int16_t grid_height = _exp_m->get_grid_height();
  
  ae_grid_cell*** pop_grid = _exp_m->get_pop_grid();
  
  // create a temporary grid to store new individuals
  ae_individual*** new_indiv_grid = new ae_individual** [grid_width];
  for ( int16_t i = 0 ; i < grid_width ; i++ )
  {
    new_indiv_grid[i] = new ae_individual* [grid_height];
  }

        
  // Do local competitions
  for ( int16_t x = 0 ; x < grid_width ; x++ )
  {
    for ( int16_t y = 0 ; y < grid_height ; y++ )
    {   
      new_indiv_grid[x][y] = calculate_local_competition( x, y );
      new_indiv_grid[x][y]->set_grid_cell( pop_grid[x][y] );
    }
  }  
  
  
  // Add the compound secreted by the individuals
  if ( _exp_m->get_with_secretion() )
  {
    double tmp_secretion; 
    for ( int16_t x = 0 ; x < grid_width ; x++ )
    {
      for ( int16_t y = 0 ; y < grid_height ; y++ )
      {
        tmp_secretion = pop_grid[x][y]->get_compound_amount() + pop_grid[x][y]->get_individual()->get_fitness_by_feature(SECRETION);        
        pop_grid[x][y]->set_compound_amount( tmp_secretion );
      }
    }
    
    // Diffusion and degradation of compound in the environment
    sp_struct->update_secretion_grid();
  }
  
  
  // Create the new generation
  ae_list<ae_individual*>* new_generation = new ae_list<ae_individual*>();
  int32_t index_new_indiv = 0;
  for ( int16_t x = 0 ; x < grid_width ; x++ )
  {
    for ( int16_t y = 0 ; y < grid_height ; y++ )
    {
      pop_grid[x][y]->set_individual( do_replication( new_indiv_grid[x][y], index_new_indiv++, x, y ) );
      #ifdef DISTRIBUTED_PRNG
        #error Not implemented yet !
        new_indiv_grid[x][y]->do_prng_jump();
      #endif
      new_generation->add( pop_grid[x][y]->get_individual() );
    }
  }
  
  // Replace the old population by the newly created one
  _exp_m->get_pop()->replace_population( new_generation );

  // delete the temporary grid
  for ( int16_t x = 0 ; x < grid_width ; x++ )
  {
    delete [] new_indiv_grid[x];    
  }
  delete [] new_indiv_grid;
  
  // randomly migrate some organisms, if necessary 
  if ( sp_struct->get_migration_number() > 0 )
  {
    sp_struct->do_random_migrations();
  }
  
  // Perform plasmid transfer
  if ( _exp_m->get_with_plasmids() && ( (_exp_m->get_prob_plasmid_HT() != 0.0) || (_exp_m->get_tune_donor_ability() != 0.0) || (_exp_m->get_tune_recipient_ability() != 0.0) ) )
  {
    int16_t x_offset, y_offset, new_x, new_y;

    
    // Shuffle the grid:
    int16_t total_size = ((grid_width)*(grid_height));
    int16_t** shuffled_table = new int16_t* [total_size];
    for ( int16_t z = 0 ; z < total_size ; z++ )
    {
      shuffled_table[z] = new int16_t[2];
      int16_t quotient = z / grid_width;
      int16_t remainder = z % grid_width;
      shuffled_table[z][0] = (int16_t) remainder;
      shuffled_table[z][1] = (int16_t) quotient;
    }
    
    for ( int16_t z = 0 ;z < total_size - 1 ; z++ )
    {
      int16_t rand_nb = _prng->random((int16_t) (total_size-z));
      int16_t* tmp=shuffled_table[z+rand_nb];
      shuffled_table[z+rand_nb]=shuffled_table[z];
      shuffled_table[z]=tmp;		
    }
    
    
    // First transfer all the plasmids, but just add them at the end of the list of the GUs
    for ( int16_t z = 0 ; z < total_size ; z++ ) // for each individual x
    {      
      int16_t x=shuffled_table[z][0];
      int16_t y=shuffled_table[z][1];
      
      for ( int16_t n = 0 ; n < 9 ; n++ ) // for each neighbour n of x
      {
        x_offset = ( n / 3 ) - 1;
        y_offset = ( n % 3 ) - 1;
        
        new_x = (x+x_offset+grid_width) % grid_width;
        new_y = (y+y_offset+grid_height) % grid_height;
        
        if ((new_x != x)||(new_y != y))
        {
          double ptransfer = _exp_m->get_prob_plasmid_HT() + _exp_m->get_tune_donor_ability()
                            * sp_struct->get_indiv_at(x, y)->get_fitness_by_feature(DONOR)
                            + _exp_m->get_tune_recipient_ability() * sp_struct->get_indiv_at(new_x, new_y)->get_fitness_by_feature(RECIPIENT) ;
          if ( _prng->random() < ptransfer ) // will x give a plasmid to n ?
          {
            if ( _exp_m->get_swap_GUs() )
            {
              sp_struct->get_indiv_at(new_x, new_y)->inject_2GUs( sp_struct->get_indiv_at(x, y) );
            }
            else
            {
              sp_struct->get_indiv_at(new_x, new_y)->inject_GU( sp_struct->get_indiv_at(x, y) );
            }
          }
        }
      }
    }
    
    for(int16_t z=0;z <total_size;z++)
    {
      delete [] shuffled_table[z];
    }
    delete [] shuffled_table;
    
	
	
    // If an individual has more than 2 GUs, we keep only the first (main chromosome) and the last one
    // and re-evaluate the individual
    for ( int16_t x = 0 ; x < grid_width ; x++ )
    {
      for ( int16_t y = 0 ; y < grid_height ; y++ )
      { 
        bool reevaluate = false;
        while ( sp_struct->get_indiv_at(x, y)->get_genetic_unit_list()->get_nb_elts() > 2 ) 
        {
          reevaluate = true;
          sp_struct->get_indiv_at(x, y)->get_genetic_unit_list()->remove(
                  sp_struct->get_indiv_at(x, y)->get_genetic_unit_list()->get_first()->get_next(), true, true);
        }
        if (reevaluate)
        {
          sp_struct->get_indiv_at(x, y)->reevaluate();
          //previous code avoided full re-evaluation, just updated the "full" proteins list:
          //_pop_grid[x][y]->get_individual()->reevaluate_after_GU_transfer(ae_common::sim->get_env());
        }
      }
    }
  }
  
  // Update the best individual
  _exp_m->get_pop()->update_best();
}

/*!
*/
void ae_selection::write_setup_file( gzFile exp_setup_file ) const
{
  // ---------------------------------------------------- Selection Parameters
  int8_t tmp_sel_scheme = _selection_scheme;
  gzwrite( exp_setup_file, &tmp_sel_scheme,      sizeof(tmp_sel_scheme) );
  gzwrite( exp_setup_file, &_selection_pressure, sizeof(_selection_pressure) );
}

/*!
*/
void ae_selection::write_setup_file( FILE* exp_setup_file ) const
{
  // TODO
}


/*!
*/
void ae_selection::save( gzFile& backup_file ) const
{
  if ( _prng == NULL )
  {
    printf( "%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  // ----------------------------------------- Pseudo-random number generator
  _prng->save( backup_file );
}

void ae_selection::load( gzFile& exp_setup_file,
                         gzFile& backup_file,
                         bool verbose )
{
  // ---------------------------------------------------- Selection parameters
  int8_t tmp_sel_scheme;
  gzread( exp_setup_file, &tmp_sel_scheme, sizeof(tmp_sel_scheme) );
  _selection_scheme = (ae_selection_scheme) tmp_sel_scheme;
  gzread( exp_setup_file, &_selection_pressure, sizeof(_selection_pressure) );

  // ----------------------------------------- Pseudo-random number generator
  _prng = new ae_jumping_mt( backup_file );
}

void ae_selection::load( FILE*& exp_setup_file,
                         gzFile& backup_file,
                         bool verbose )
{
  // TODO
}


// =================================================================
//                           Protected Methods
// =================================================================
void ae_selection::compute_prob_reprod( void )
{
  if ( _prob_reprod != NULL )
  {
    delete [] _prob_reprod;
  }
  
  int32_t nb_indivs = _exp_m->get_pop()->get_nb_indivs();
  _prob_reprod = new double[nb_indivs];

  if ( _selection_scheme == RANK_LINEAR )
  {
    // The probability of reproduction for an individual is given by
    // ( 2-SP + 2 * (SP-1) * (R-1)/(N-1) ) / N
    // With :
    //      SP : selective pressure. Linear ranking allows values of SP in [1.0, 2.0].
    //      R  : the rank of the individual in the population (1 for the worst individual)
    //      N  : the number of individuals in the population
    //
    // We can transform this expression into (2-SP)/N + ((2*(SP-1)) / (N*(N-1))) * (R-1)
    // Furthermore, (R-1) is given directly by <i> (the index of our probability table)
    //
    // probs[0] will hence be given by (2-SP)/N
    // probs[i+1] can then be expressed by probs[i] + (2*(SP-1)) / (N*(N-1))

    double increment  = (2 * (_selection_pressure-1)) / (nb_indivs * (nb_indivs-1));
    _prob_reprod[0]   = (2 - _selection_pressure) / nb_indivs;

    for ( int32_t i = 1 ; i < nb_indivs ; i++ )
    {
      _prob_reprod[i] = _prob_reprod[i-1] + increment;
    }

    // No need to normalize: The sum is always 1 for linear ranking
  }
  else if ( _selection_scheme == RANK_EXPONENTIAL )
  {
    // The probability of reproduction for an individual is given by
    // ( (SP-1) * SP^(N-R) ) / ( SP^N - 1 )
    // Which is equivalent to
    // ( (SP-1) * SP^N ) / ( (SP^N - 1) * SP^R )
    // With :
    //      SP : selective pressure. Exponential ranking allows values of SP in ]0.0, 1.0[
    //      R  : the rank of the individual in the population (1 for the worst individual)
    //      N  : the number of individuals in the population
    //
    // NB : The only rank-dependent term is SP^R
    //
    // Because we don't allow ex-aequo,
    // probs[i+1] can hence be expressed as (probs[i] / SP)
    // We will hence compute probs[0] with the original formula and infer the remaining values

    double SP_N = pow(_selection_pressure, nb_indivs); // SP^N
    _prob_reprod[0] = ( (_selection_pressure - 1) * SP_N ) /
                      ( (SP_N - 1) * _selection_pressure );

    for ( int32_t i = 1 ; i < nb_indivs ; i++ )
    {
      _prob_reprod[i] = _prob_reprod[i-1] / _selection_pressure;
    }

    // No need to normalize: We don't allow ex-aequo
  }
  else if ( _selection_scheme == FITNESS_PROPORTIONATE ) // Fitness Proportionate
  {
    // The probability of reproduction for an individual is given by
    // exp( -SP * gap ) / sum of this measure on all individuals
    //    SP : selective pressure. Fitness proportionate allows values of SP in ]0, +inf[
    //                             The closer SP to 0, the closer the selection to being linear.

    double* fitnesses = new double[nb_indivs];
    double  sum       = 0;
    
    ae_list_node<ae_individual*>* indiv_node = _exp_m->get_pop()->get_indivs()->get_first();
    ae_individual*  indiv       = NULL;
   
    for ( int32_t i = 0 ; i < nb_indivs ; i++ )
    {
      indiv = indiv_node->get_obj();
      fitnesses[i] = indiv->get_fitness();
      sum += fitnesses[i];
      indiv_node = indiv_node->get_next();
    }

    for ( int32_t i = 0 ; i < nb_indivs ; i++ )
    {
      _prob_reprod[i] = fitnesses[i] / sum;
    }

    delete [] fitnesses;
  }
  else if ( _selection_scheme == FITTEST) //  Fittest individual
  {
    printf( "ERROR, fittest selection scheme is meant to be used for spatially structured populations %s:%d\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  else
  {
    printf( "ERROR, invalid selection scheme in file %s:%d\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
}

void ae_selection::compute_local_prob_reprod( void )
{
  int16_t neighborhood_size = 9;
  
  if ( _prob_reprod != NULL )
  {
    printf ("Warning, already defined %s:%d\n", __FILE__, __LINE__);
    delete [] _prob_reprod;
  }
  
  _prob_reprod = new double[neighborhood_size];
  
  if ( _selection_scheme == RANK_LINEAR )
  {
    double increment = (2 * (_selection_pressure-1)) / (neighborhood_size * (neighborhood_size-1));
    double init_prob = (2 - _selection_pressure) / neighborhood_size;
    
    for ( int16_t i = 0 ; i < neighborhood_size ; i++ )
    {
      _prob_reprod[i] = init_prob + increment * i;
    }   
  }
  else if ( _selection_scheme == RANK_EXPONENTIAL )
  {
    double SP_N = pow(_selection_pressure, neighborhood_size); 
    _prob_reprod[0] = ( (_selection_pressure - 1) * SP_N ) /
    ( (SP_N - 1) * _selection_pressure );
    
    for ( int16_t i = 1 ; i < neighborhood_size ; i++ )
    {
      _prob_reprod[i] =  _prob_reprod[i-1] /  _selection_pressure; 
    }
  }
  else if ( _selection_scheme == FITTEST) //  Fittest individual
  {
    for ( int16_t i = 0 ; i < neighborhood_size-1 ; i++ )
    {
      _prob_reprod[i] = 0.;
    }
    _prob_reprod[neighborhood_size-1] = 1.;
  }
  else if ( _selection_scheme == FITNESS_PROPORTIONATE ) // Fitness Proportionate
  {
    printf( "ERROR, this function is not intented to be use with this selection scheme %s:%d\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  else
  {
    printf( "ERROR, invalid selection scheme in file %s:%d\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
}

ae_individual* ae_selection::do_replication( ae_individual* parent, int32_t index, int16_t x /*= 0*/, int16_t y /*= 0*/ )
{
  ae_individual* new_indiv = NULL;

  // ===========================================================================
  //  1) Copy parent
  // ===========================================================================
  #ifdef __NO_X
    #ifndef __REGUL
      new_indiv = new ae_individual( parent, index, parent->get_mut_prng(), parent->get_stoch_prng() );
    #else
      new_indiv = new ae_individual_R( dynamic_cast<ae_individual_R*>(parent), index, parent->get_mut_prng(), parent->get_stoch_prng() );
    #endif
  #elif defined __X11
    #ifndef __REGUL
      new_indiv = new ae_individual_X11( dynamic_cast<ae_individual_X11*>(parent), index, parent->get_mut_prng(), parent->get_stoch_prng() );
    #else
      new_indiv = new ae_individual_R_X11( dynamic_cast<ae_individual_R_X11*>(parent), index, parent->get_mut_prng(), parent->get_stoch_prng() );
    #endif
  #endif
  
  
  // ===========================================================================
  //  2) Set the new individual's location on the grid
  //     (needed if the population is structured)
  // ===========================================================================
  if ( _exp_m->is_spatially_structured() && (x != -1) )
  {
    new_indiv->set_grid_cell( _exp_m->get_spatial_structure()->get_grid_cell( x, y ) );
  }
  
  
  // ===========================================================================
  //  3) Perform transfer
  // ===========================================================================
  int32_t nb_indivs = _exp_m->get_pop()->get_nb_indivs();
  
  // -----------------------------------
  //  a) Insertion transfer
  // -----------------------------------
  if ( _prng->random() < _exp_m->get_HT_ins_rate() )
  {
    // Insertion transfer
    // Requirements:
    //    * A circular exogenote => an alignment on the donor chromosome
    //    * An alignment between the exogenote and the endogenote
    
    // 1) Draw a random donor (uniform drawing).
    // We use the rank because indivs are sorted by rank (1 for the worst, POP_SIZE for the best).
    ae_individual * donor = NULL;
    do donor = _exp_m->get_pop()->get_indiv_by_rank( _prng->random( nb_indivs ) + 1 );
    while ( donor == parent );
    
    // 2) Look for an alignment within the donor genome
    ae_vis_a_vis* alignment_1   = NULL;
    ae_dna*       donor_dna     = donor->get_genetic_unit( 0 )->get_dna();
    ae_dna*       new_indiv_dna = new_indiv->get_genetic_unit( 0 )->get_dna();
    int32_t       nb_pairs_1    = (int32_t)( ceil( donor_dna->get_length() * parent->get_neighbourhood_rate() ) );
    
    alignment_1 = donor_dna->search_alignment( donor_dna, nb_pairs_1, DIRECT );
    
    if ( alignment_1 != NULL )
    {
      // 3) Make a copy of the sequence to be transferred (the exogenote)
      ae_genetic_unit* exogenote = donor_dna->copy_into_new_GU( alignment_1->get_i_1(), alignment_1->get_i_2() );
      
      // 4) Look for an alignments between the exogenote and the endogenote
      ae_vis_a_vis* alignment_2 = NULL;
      int32_t       nb_pairs_2  = (int32_t)( ceil( new_indiv_dna->get_length() * parent->get_neighbourhood_rate() ) );
      
      alignment_2 = exogenote->get_dna()->search_alignment( new_indiv_dna, nb_pairs_2, BOTH_SENSES );
      
      if ( alignment_2 != NULL )
      {
        //~ FILE * logfile = fopen( "transfer.out", "w" );
        //~ fprintf( logfile, "DONOR:\n%s\n\n\n", donor_dna->get_data() );
        //~ fprintf( logfile, "EXOGENOTE:\n%s\n\n\n", exogenote->get_dna()->get_data() );
        //~ fprintf( logfile, "RECIPIENT:\n%s\n\n\n", new_indiv_dna->get_data() );
        //~ #ifdef DEBUG
          //~ ae_common::sim->get_logs()->flush();
          //~ new_indiv->assert_promoters();
          //~ new_indiv->assert_promoters_order();
        //~ #endif
        
        //~ printf( "new indiv length: %"PRId32" parent: %"PRId32" donor: %"PRId32" exogenote: %"PRId32"\n",
                //~ new_indiv_dna->get_length(), parent->get_genetic_unit( 0 )->get_dna()->get_length(), donor_dna->get_length(), exogenote->get_dna()->get_length() );
        
        //~ printf( "************************************************************\n" );
        //~ donor->get_genetic_unit( 0 )->print_rnas();
        //~ printf( "************************************************************\n" );
        //~ exogenote->print_rnas();
        //~ printf( "************************************************************\n" );
        
        int32_t genome_length_before  = new_indiv_dna->get_length();
        int32_t genome_length_after   = new_indiv_dna->get_length() + exogenote->get_dna()->get_length();
        
        if ( genome_length_after > new_indiv->get_max_genome_length() )
        {
          #warning LOG
          //~ if ( _exp_m->is_logged( LOG_BARRIER ) == true )
          //~ {
            //~ // Write an entry in the barrier log file
            //~ fprintf(  ae_common::rec_params->get_log( LOG_BARRIER ), "%"PRId32" %"PRId32" INS_TRANSFER %"PRId32" %"PRId32" %"PRId32"\n",
                      //~ _exp_m->get_num_gener(),
                      //~ new_indiv->get_index_in_population(),
                      //~ exogenote->get_dna()->get_length(),
                      //~ 0,
                      //~ genome_length_before );
          //~ }
        }
        else
        {
          new_indiv_dna->insert_GU( exogenote, alignment_2->get_i_2(), alignment_2->get_i_1(), (alignment_2->get_sense() == INDIRECT) );
          //~ fprintf( logfile, "RESULT:\n%s\n\n\n", new_indiv_dna->get_data() );
          //~ fflush( logfile );
          

          // Write a line in transfer logfile
          #warning LOG
          //~ if ( ae_common::rec_params->is_logged( LOG_TRANSFER ) == true )
          //~ {
            //~ fprintf(  ae_common::rec_params->get_log( LOG_TRANSFER ),
                      //~ "%"PRId32" %"PRId32" %"PRId32" %"PRId8" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId16" %"PRId32" %"PRId32" %"PRId16"\n",
                      //~ _exp_m->get_num_gener(),
                      //~ new_indiv->get_index_in_population(),
                      //~ donor->get_index_in_population(),
                      //~ 0, // Transfer type
                      //~ exogenote->get_dna()->get_length(),
                      //~ 0,
                      //~ genome_length_before,
                      //~ new_indiv_dna->get_length(),
                      //~ alignment_1->get_i_1(),
                      //~ alignment_1->get_i_2(),
                      //~ alignment_1->get_score(),
                      //~ alignment_2->get_i_1(),
                      //~ alignment_2->get_i_2(),
                      //~ alignment_2->get_score() );
          //~ }
    
          #ifdef BIG_DEBUG
            ae_common::sim->get_logs()->flush();
            new_indiv->assert_promoters();
            new_indiv->assert_promoters_order();
          #endif
        }
        
        delete alignment_2;
      }
      
      delete exogenote;
      delete alignment_1;
    }
  }
  
  // -----------------------------------
  //  b) Replacement transfer
  // -----------------------------------
  if ( _prng->random() < _exp_m->get_HT_repl_rate() )
  {
    // Replacement transfer
    // Requirements:
    //    * 2 distinct alignments between the (linear) exogenote and the endogenote
    
    // 1) Draw a random donor (uniform drawing).
    // We use the rank because indivs are sorted by rank (1 for the worst, POP_SIZE for the best).
    ae_individual * donor = NULL;
    do donor = _exp_m->get_pop()->get_indiv_by_rank( _prng->random( nb_indivs ) + 1 );
    while ( donor == parent );
    
    // 2) Look for an alignment between the parent genome and the donor genome
    ae_vis_a_vis* alignment_1   = NULL;
    ae_vis_a_vis* alignment_2   = NULL;
    ae_dna*       parent_dna    = parent->get_genetic_unit( 0 )->get_dna();
    ae_dna*       donor_dna     = donor->get_genetic_unit( 0 )->get_dna();
    ae_dna*       new_indiv_dna = new_indiv->get_genetic_unit( 0 )->get_dna();
    ae_sense      sense         = (_prng->random() < 0.5) ? DIRECT : INDIRECT;
    int32_t       nb_pairs_1    = (int32_t)( ceil( new_indiv_dna->get_length() * parent->get_neighbourhood_rate() ) );
    int32_t       nb_pairs_2    = (int32_t)( ceil( new_indiv_dna->get_length() * parent->get_neighbourhood_rate() ) );
    
    alignment_1 = parent_dna->search_alignment( donor_dna, nb_pairs_1, sense );
    
    if ( alignment_1 != NULL )
    {
      // Look for a second alignement between the parent and the donor (must be different from alignment_1)
      while ( alignment_2 == NULL && nb_pairs_2 > 0 )
      {
        alignment_2 = parent_dna->search_alignment( donor_dna, nb_pairs_2, sense );
        
        // Forbid the replacement of the whole genome of the parent
        if ( alignment_2 != NULL && alignment_2->get_i_1() == alignment_1->get_i_1() )
        {
          delete alignment_2;
          alignment_2 = NULL;
        }
      }
    
    
      // If both alignments were found, proceed to the transfer
      if ( alignment_2 != NULL )
      {
        int32_t genome_length_before  = new_indiv_dna->get_length();
        int32_t genome_length_after   = new_indiv_dna->get_length();
        int32_t exogenote_length      = ae_utils::mod( alignment_2->get_i_2() - alignment_1->get_i_2() - 1, donor_dna->get_length() ) + 1;
        int32_t replaced_seq_length   = ae_utils::mod( alignment_2->get_i_1() - alignment_1->get_i_1() - 1, genome_length_before ) + 1;
        
        if ( genome_length_after < new_indiv->get_min_genome_length() || genome_length_after > new_indiv->get_max_genome_length() )
        {
          #warning LOG
          //~ if ( ae_common::rec_params->is_logged( LOG_BARRIER ) == true )
          //~ {
            //~ // Write an entry in the barrier log file
            //~ fprintf(  ae_common::rec_params->get_log( LOG_BARRIER ), "%"PRId32" %"PRId32" REPL_TRANSFER %"PRId32" %"PRId32" %"PRId32"\n",
                      //~ _exp_m->get_num_gener(),
                      //~ new_indiv->get_index_in_population(),
                      //~ exogenote_length,
                      //~ replaced_seq_length,
                      //~ genome_length_before );
          //~ }
        }
        
        // 3) Make a copy of the sequence to be transferred (the exogenote)
        ae_genetic_unit* exogenote = NULL;
        if ( sense == DIRECT )
        {
          exogenote = donor_dna->copy_into_new_GU( alignment_1->get_i_2(), alignment_2->get_i_2() );
        }
        else
        {
          exogenote = donor_dna->copy_into_new_GU( alignment_2->get_i_2(), alignment_1->get_i_2() );
        }
        
        // Delete the sequence to be replaced
        new_indiv_dna->do_deletion( alignment_1->get_i_1(), alignment_2->get_i_1() );
        
        //~ printf( "intermediate genome length : %"PRId32"\n", new_indiv_dna->get_length() );
        
        // Insert the transfered sequence
        if ( alignment_1->get_i_1() < alignment_2->get_i_1() )
        {
          new_indiv_dna->insert_GU( exogenote, alignment_1->get_i_1(), 0, sense == INDIRECT );
        }
        else
        {
          new_indiv_dna->insert_GU( exogenote, 0, 0, sense == INDIRECT );
        }
        

        // Write a line in transfer logfile
        #warning LOG
        //~ if ( ae_common::rec_params->is_logged( LOG_TRANSFER ) == true )
        //~ {
          //~ fprintf(  ae_common::rec_params->get_log( LOG_TRANSFER ),
                    //~ "%"PRId32" %"PRId32" %"PRId32" %"PRId8" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId16" %"PRId32" %"PRId32" %"PRId16"\n",
                    //~ _exp_m->get_num_gener(),
                    //~ new_indiv->get_index_in_population(),
                    //~ donor->get_index_in_population(),
                    //~ 1, // Transfer type
                    //~ exogenote->get_dna()->get_length(),
                    //~ replaced_seq_length,
                    //~ genome_length_before,
                    //~ new_indiv_dna->get_length(),
                    //~ alignment_1->get_i_1(),
                    //~ alignment_1->get_i_2(),
                    //~ alignment_1->get_score(),
                    //~ alignment_2->get_i_1(),
                    //~ alignment_2->get_i_2(),
                    //~ alignment_2->get_score() );
        //~ }
        
        delete exogenote;
        delete alignment_2;
      }
      
      delete alignment_1;
    }
  }
  
  
  // ===========================================================================
  //  4) Perform rearrangements and mutations
  // ===========================================================================
  if ( ! new_indiv->get_allow_plasmids() )
  {
    ae_genetic_unit* chromosome = new_indiv->get_genetic_unit_list()->get_first()->get_obj();
    
    chromosome->get_dna()->perform_mutations();
    
    if ( new_indiv->get_replic_report() != NULL )
    {
      new_indiv->get_replic_report()->get_dna_replic_reports()->add( chromosome->get_dna()->get_replic_report() );
    }
  }
  else
  {
    // For each GU, apply mutations
    ae_list_node<ae_genetic_unit*>* gen_unit_node = NULL;
    ae_genetic_unit* gen_unit = NULL;
    
    // Randomly determine the order in which the GUs will undergo mutations
    bool inverse_order = (_prng->random((int32_t) 2) < 0.5);
    
    if ( ! inverse_order )
    // Apply mutations in normal GU order
    {
      gen_unit_node = new_indiv->get_genetic_unit_list()->get_first();
      
      while ( gen_unit_node != NULL )
      {
        gen_unit = gen_unit_node->get_obj();
        
        gen_unit->get_dna()->perform_mutations();
        
        if ( new_indiv->get_replic_report() != NULL )
        {
          new_indiv->get_replic_report()->get_dna_replic_reports()->add( gen_unit->get_dna()->get_replic_report() );
        }

        gen_unit_node = gen_unit_node->get_next();
      }
    }
    else
    // Apply mutations in inverse GU order
    {
      gen_unit_node = new_indiv->get_genetic_unit_list()->get_last();
      
      while ( gen_unit_node != NULL )
      {
        gen_unit = gen_unit_node->get_obj();
        
        gen_unit->get_dna()->perform_mutations();
        
        if ( new_indiv->get_replic_report() != NULL )
        {
          new_indiv->get_replic_report()->get_dna_replic_reports()->add( gen_unit->get_dna()->get_replic_report() );
        }

        gen_unit_node = gen_unit_node->get_prev();
      }
    }
  }


  
  // ===========================================================================
  //  5) Evaluate new individual
  // ===========================================================================
  new_indiv->evaluate();
  
  
  // ===========================================================================
  //  6) Compute statistics
  // ===========================================================================
  new_indiv->compute_statistical_data();
  
  
  #ifdef BIG_DEBUG
  ae_common::sim->get_logs()->flush();
  new_indiv->assert_promoters();
  #endif
  
  
  return new_indiv;
}

ae_individual* ae_selection::calculate_local_competition ( int16_t x, int16_t y )
{
  ae_spatial_structure* sp_struct = _exp_m->get_spatial_structure();
  
  int16_t neighborhood_size = 9;
  int16_t grid_width  = sp_struct->get_grid_width();
  int16_t grid_height = sp_struct->get_grid_height();
  int16_t cur_x;
  int16_t cur_y;

  // construct a temporarry local array of fitness values
  double *  local_fit_array   = new double[neighborhood_size];
  double *  sort_fit_array    = new double[neighborhood_size];
  int16_t * initial_location  = new int16_t[neighborhood_size];
  double *  probs             = new double[neighborhood_size];
  int16_t   count             = 0;
  double    sum_local_fit     = 0.0;
  
  for ( int8_t i = -1 ; i < 2 ; i++ )
  {
    for ( int8_t j = -1 ; j < 2 ; j++ )
    {
      cur_x = ( x + i + grid_width )  % grid_width;
      cur_y = ( y + j + grid_height ) % grid_height;
      local_fit_array[count]  = sp_struct->get_indiv_at( cur_x, cur_y )->get_fitness();
      sort_fit_array[count]   = local_fit_array[count];
      initial_location[count] = count;
      sum_local_fit += local_fit_array[count];
      count++;
    }
  }

  // Do the competitions between the individuals, based on one of the 4 methods: 
  // 1. Rank linear
  // 2. Rank exponential
  // 3. Fitness proportionate
  // 4. Fittest individual
  
  // Any rank based selection
  switch ( _selection_scheme )
  {
    case RANK_LINEAR :
    case RANK_EXPONENTIAL :
    case FITTEST :
    {
      // First we sort the local fitness values using bubble sort :
      // we sort by increasing order, so the first element will have the worst fitness.
      bool swaped = true;
      int16_t loop_length = 8; 
      double  tmp_holder;
      int16_t tmp_holder2;
      while ( swaped == true )
      {   
        swaped = false;
        for ( int16_t i = 0 ; i < loop_length ; i++ )
        {
          //if the first is higher than the second,  exchange them
          if ( sort_fit_array[i] > sort_fit_array[i+1] )
          {
            tmp_holder = sort_fit_array[i];
            sort_fit_array[i] = sort_fit_array[i+1];
            sort_fit_array[i+1] = tmp_holder;
            
            tmp_holder2 = initial_location[i];
            initial_location[i] = initial_location[i+1];
            initial_location[i+1] = tmp_holder2;
            
            swaped = true;
          }
        }
        
        loop_length = loop_length - 1;
      }
      
      
      // Then we use the already computed probabilities
      for ( int16_t i = 0 ; i < neighborhood_size ; i++ )
      {
        probs[initial_location[i]] = _prob_reprod[i];
      }
      
      break;
    }
    // Fitness proportionate selection
    case FITNESS_PROPORTIONATE :
    {
      for( int16_t i = 0 ; i < neighborhood_size ; i++ )
      {
        probs[i] = local_fit_array[i]/sum_local_fit;
      }
      
      break;
    }
    default :
    {
      printf( "ERROR, invalid selection scheme in file %s:%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    }
  }
  
  
  // pick one organism to reproduce, based on probs[] calculated above, using roulette selection
  int8_t found_org = _prng->roulette_random( probs, 9 );
  
  int16_t x_offset = ( found_org / 3 ) - 1;
  int16_t y_offset = ( found_org % 3 ) - 1;
  
  delete [] local_fit_array; 
  delete [] sort_fit_array; 
  delete [] initial_location;
  delete [] probs; 
  
  return sp_struct->get_indiv_at( (x+x_offset+grid_width)  % grid_width,
                                  (y+y_offset+grid_height) % grid_height );
}

// =================================================================
//                          Non inline accessors
// =================================================================
