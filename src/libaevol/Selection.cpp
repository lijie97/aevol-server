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
#include "Selection.h"

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __TBB
#include <tbb/task_group.h>
#include <tbb/parallel_sort.h>
#endif

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include "ExpManager.h"
#include "VisAVis.h"

#ifdef __NO_X
  #ifndef __REGUL
    #include "Individual.h"
  #else
    #include "raevol/Individual_R.h"
  #endif
#elif defined __X11
  #ifndef __REGUL
    #include "Individual_X11.h"
  #else
    #include "raevol/Individual_R_X11.h"
  #endif
#endif


namespace aevol {


//##############################################################################
//
//                              Class Selection
//
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Selection::Selection(ExpManager* exp_m) {
  exp_m_ = exp_m;

  // ----------------------------------------- Pseudo-random number generator
  prng_ = NULL;

  // -------------------------------------------------------------- Selection
  selection_scheme_   = RANK_EXPONENTIAL;
  selection_pressure_ = 0.998;

  // --------------------------- Probability of reproduction of each organism
  prob_reprod_ = NULL;
}

// =================================================================
//                             Destructors
// =================================================================
Selection::~Selection() {
  delete [] prob_reprod_;
}

// =================================================================
//                            Public Methods
// =================================================================
void Selection::step_to_next_generation() {
  // To create the new generation, we must create nb_indivs new individuals
  // (offspring) and "kill" the existing ones.
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

  if (prng_ == NULL) {
    printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  // -------------------------------------------------------------------------------
  // 1) Compute the probability of reproduction of each individual in the population
  // -------------------------------------------------------------------------------
  #ifndef FIXED_POPULATION_SIZE
    #error this method is not ready for variable population size
    compute_local_prob_reprod();
  #else
    // The function compute_local_prob_reprod creates and fills the array prob_reprod_, which is telling us the probability of being picked for reproduction according to the rank of an individual in its neighboorhood.
    // It is only usefull when selection is rank based. When selection scheme is FITNESS_PROPORTIONATE, we do not need to call it.
    // It shoud only be called once in the simulation and not at each generation. So if prob_reprod_ already exists we do not need to call it.
    if ((selection_scheme_ != FITNESS_PROPORTIONATE) && (prob_reprod_ == NULL)) {
      compute_local_prob_reprod();
    }
  #endif

  if (prng_ == NULL) {
    printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  // Create proxies
  World* world = exp_m_->world();
  int16_t grid_width  = world->width();
  int16_t grid_height = world->height();
  GridCell*** pop_grid = exp_m_->grid();

  // create a temporary grid to store the reproducers
  Individual *** reproducers = new Individual ** [grid_width];
  for (int16_t i = 0 ; i < grid_width ; i++) {
    reproducers[i] = new Individual* [grid_height];
  }


  // Do local competitions
  for (int16_t x = 0 ; x < grid_width ; x++) {
    for (int16_t y = 0 ; y < grid_height ; y++) {
      reproducers[x][y] = do_local_competition(x, y);
    }
  }


  // TODO : Why is that not *after* the creation of the new population ?
  // Add the compound secreted by the individuals
  if (exp_m_->with_secretion()) {
    for (int16_t x = 0 ; x < grid_width ; x++) {
      for (int16_t y = 0 ; y < grid_height ; y++) {
        pop_grid[x][y]->set_compound_amount(
            pop_grid[x][y]->compound_amount() +
            pop_grid[x][y]->individual()->fitness_by_feature(SECRETION));
      }
    }

    // Diffusion and degradation of compound in the habitat
    world->update_secretion_grid();
  }


  // Create the new generation
  std::list<Individual*> old_generation = exp_m_->indivs();;

  #ifndef __TBB
  std::list<Individual*> new_generation;
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(dynamic)
  #endif
  for (int16_t x = 0 ; x < grid_width ; x++)
    for (int16_t y = 0 ; y < grid_height ; y++)
      do_replication(reproducers[x][y], x * grid_height + y, x, y);

  for (int16_t x = 0 ; x < grid_width ; x++)
    for (int16_t y = 0 ; y < grid_height ; y++)
      new_generation.emplace_back(pop_grid[x][y]->individual());
  #else
  std::vector<Individual*> new_generation;
  tbb::task_group tgroup;

  for (int16_t x = 0 ; x < grid_width ; x++)
    for (int16_t y = 0 ; y < grid_height ; y++)
      tgroup.run([=] {do_replication(reproducers[x][y], x * grid_height + y, x, y);});
  tgroup.wait();

  for (int16_t x = 0 ; x < grid_width ; x++)
    for (int16_t y = 0 ; y < grid_height ; y++)
      new_generation.push_back(pop_grid[x][y]->get_individual());

  #endif

  // delete the temporary grid and the parental generation
  for (int16_t x = 0 ; x < grid_width ; x++) {
    delete [] reproducers[x];
  }
  delete [] reproducers;
  for (auto indiv : old_generation) {
    delete indiv;
  }

  // Compute the rank of each individual
  #ifndef __TBB
  new_generation.sort([](Individual* lhs, Individual* rhs) {
                       return lhs->fitness() < rhs->fitness();
                     });
  #else
  tbb::parallel_sort(new_generation.begin(),new_generation.end(), [](Individual* lhs, Individual* rhs) {
      return lhs->get_fitness() < rhs->get_fitness();
  });
  #endif
  int rank = 1;
  for (Individual* indiv : new_generation) {
    indiv->set_rank(rank++);
  }

  // randomly migrate some organisms, if necessary
  world->MixIndivs();

  PerformPlasmidTransfers();

  // Update the best individual
  exp_m_->update_best();

  // Notify observers of the end of the generation
  notifyObservers(END_GENERATION);
}

void Selection::PerformPlasmidTransfers() {
  if (exp_m_->with_plasmids() &&
      ((exp_m_->prob_plasmid_HT() != 0.0) ||
        (exp_m_->tune_donor_ability() != 0.0) ||
        (exp_m_->tune_recipient_ability() != 0.0))) {
    // Create proxies
    World* world = exp_m_->world();
    int16_t grid_width  = world->width();
    int16_t grid_height = world->height();

    int16_t x_offset, y_offset, new_x, new_y;

    // Shuffle the grid:
    int16_t total_size = ((grid_width)*(grid_height));
    int16_t** shuffled_table = new int16_t* [total_size];
    for (int16_t z = 0 ; z < total_size ; z++) {
      shuffled_table[z] = new int16_t[2];
      int16_t quotient = z / grid_width;
      int16_t remainder = z % grid_width;
      shuffled_table[z][0] = (int16_t) remainder;
      shuffled_table[z][1] = (int16_t) quotient;
    }

    for (int16_t z = 0 ;z < total_size - 1 ; z++) {
      int16_t rand_nb = prng_->random((int16_t) (total_size-z));
      int16_t* tmp=shuffled_table[z+rand_nb];
      shuffled_table[z+rand_nb]=shuffled_table[z];
      shuffled_table[z]=tmp;
    }


    // First transfer all the plasmids, but just add them at the end of the list of the GUs
    for (int16_t z = 0 ; z < total_size ; z++) { // for each individual x
      int16_t x=shuffled_table[z][0];
      int16_t y=shuffled_table[z][1];

      for (int16_t n = 0 ; n < 9 ; n++) { // for each neighbour n of x
        x_offset = (n / 3) - 1;
        y_offset = (n % 3) - 1;

        new_x = (x+x_offset+grid_width) % grid_width;
        new_y = (y+y_offset+grid_height) % grid_height;

        if ((new_x != x)||(new_y != y)) {
          double ptransfer = exp_m_->prob_plasmid_HT() + exp_m_->tune_donor_ability()
                            * world->indiv_at(x, y)->fitness_by_feature(DONOR)
                            +
            exp_m_->tune_recipient_ability() * world->indiv_at(new_x, new_y)->fitness_by_feature(RECIPIENT) ;
          if (prng_->random() < ptransfer) { // will x give a plasmid to n ?
            if (exp_m_->swap_GUs()) {
              world->indiv_at(new_x, new_y)->inject_2GUs(world->indiv_at(x, y));
            }
            else {
              world->indiv_at(new_x, new_y)->inject_GU(world->indiv_at(x, y));
            }
          }
        }
      }
    }

    for(int16_t z=0;z <total_size;z++) {
      delete [] shuffled_table[z];
    }
    delete [] shuffled_table;



    // If an individual has more than 2 GUs, we keep only the first (main chromosome) and the last one
    // and re-evaluate the individual
    for (int16_t x = 0 ; x < grid_width ; x++) {
      for (int16_t y = 0 ; y < grid_height ; y++) {
        bool reevaluate = (world->indiv_at(x, y)->nb_genetic_units() > 2);
        world->indiv_at(x, y)->drop_nested_genetic_units();
        if (reevaluate)
          world->indiv_at(x, y)->Reevaluate();
      }
    }
  }
}

/*!
*/
void Selection::write_setup_file(gzFile exp_setup_file) const {
  // ---------------------------------------------------- Selection Parameters
  int8_t tmp_sel_scheme = selection_scheme_;
  gzwrite(exp_setup_file, &tmp_sel_scheme,      sizeof(tmp_sel_scheme));
  gzwrite(exp_setup_file, &selection_pressure_, sizeof(selection_pressure_));
}

/*!
*/
void Selection::save(gzFile& backup_file) const {
  if (prng_ == NULL) {
    printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  // ----------------------------------------- Pseudo-random number generator
  prng_->save(backup_file);
}

void Selection::load(gzFile& exp_setup_file,
                     gzFile& backup_file,
                     bool verbose) {
  // ---------------------------------------------------- Selection parameters
  int8_t tmp_sel_scheme;
  gzread(exp_setup_file, &tmp_sel_scheme, sizeof(tmp_sel_scheme));
  selection_scheme_ = (SelectionScheme) tmp_sel_scheme;
  gzread(exp_setup_file, &selection_pressure_, sizeof(selection_pressure_));

  // ----------------------------------------- Pseudo-random number generator
#if __cplusplus == 201103L
  prng_ = make_unique<JumpingMT>(backup_file);
#else
  prng_ = std::make_unique<JumpingMT>(backup_file);
#endif
}


// =================================================================
//                           Protected Methods
// =================================================================
void Selection::compute_prob_reprod() { // non spatially structured only
  if (prob_reprod_ != NULL) { // TODO <david.parsons@inria.fr> remove
    delete [] prob_reprod_;
  }

  int32_t nb_indivs = exp_m_->nb_indivs();
  prob_reprod_ = new double[nb_indivs];

  if (selection_scheme_ == RANK_LINEAR) {
    // The probability of reproduction for an individual is given by
    // (2-SP + 2 * (SP-1) * (R-1)/(N-1)) / N
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

    double increment = (2 * (selection_pressure_-1)) / (nb_indivs * (nb_indivs-1));
    prob_reprod_[0]  = (2 - selection_pressure_) / nb_indivs;

    for (int32_t i = 1 ; i < nb_indivs ; i++) {
      prob_reprod_[i] = prob_reprod_[i-1] + increment;
    }

    // No need to normalize: The sum is always 1 for linear ranking
  }
  else if (selection_scheme_ == RANK_EXPONENTIAL) {
    // The probability of reproduction for an individual is given by
    // ((SP-1) * SP^(N-R)) / (SP^N - 1)
    // Which is equivalent to
    // ((SP-1) * SP^N) / ((SP^N - 1) * SP^R)
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

    double SP_N = pow(selection_pressure_, nb_indivs); // SP^N
    prob_reprod_[0] = ((selection_pressure_ - 1) * SP_N) /
                      ((SP_N - 1) * selection_pressure_);

    for (int32_t i = 1 ; i < nb_indivs ; i++) {
      prob_reprod_[i] = prob_reprod_[i-1] / selection_pressure_;
    }

    // No need to normalize: We don't allow ex-aequo
  }
  else if (selection_scheme_ == FITNESS_PROPORTIONATE) {
    // The probability of reproduction for an individual is given by
    // exp(-SP * gap) / sum of this measure on all individuals
    //    SP : selective pressure. Fitness proportionate allows values of SP in ]0, +inf[
    //                             The closer SP to 0, the closer the selection to being linear.

    double* fitnesses = new double[nb_indivs];
    double  sum       = 0;

    size_t i = 0;
    for (const auto& indiv: exp_m_->indivs()) {
      fitnesses[i] = indiv->fitness();
      sum += fitnesses[i];
      ++i;
    }

    for (int32_t i = 0 ; i < nb_indivs ; i++) {
      prob_reprod_[i] = fitnesses[i] / sum;
    }

    delete [] fitnesses;
  }
  else if (selection_scheme_ == FITTEST) {
    printf("ERROR, fittest selection scheme is meant to be used for spatially structured populations %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  else {
    printf("ERROR, invalid selection scheme in file %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

void Selection::compute_local_prob_reprod() {
  int16_t neighborhood_size = 9;

  if (prob_reprod_ != NULL) {
    printf ("Warning, already defined %s:%d\n", __FILE__, __LINE__);
    delete [] prob_reprod_;
  }

  prob_reprod_ = new double[neighborhood_size];

  if (selection_scheme_ == RANK_LINEAR) {
    double increment = (2 * (selection_pressure_-1)) / (neighborhood_size * (neighborhood_size-1));
    double init_prob = (2 - selection_pressure_) / neighborhood_size;

    for (int16_t i = 0 ; i < neighborhood_size ; i++) {
      prob_reprod_[i] = init_prob + increment * i;
    }
  }
  else if (selection_scheme_ == RANK_EXPONENTIAL) {
    double SP_N = pow(selection_pressure_, neighborhood_size);
    prob_reprod_[0] = ((selection_pressure_ - 1) * SP_N) /
    ((SP_N - 1) * selection_pressure_);

    for (int16_t i = 1 ; i < neighborhood_size ; i++) {
      prob_reprod_[i] =  prob_reprod_[i-1] /  selection_pressure_;
    }
  }
  else if (selection_scheme_ == FITTEST) {
    for (int16_t i = 0 ; i < neighborhood_size-1 ; i++) {
      prob_reprod_[i] = 0.;
    }
    prob_reprod_[neighborhood_size-1] = 1.;
  }
  else if (selection_scheme_ == FITNESS_PROPORTIONATE) {
    printf("ERROR, this function is not intented to be use with this selection scheme %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  else {
    printf("ERROR, invalid selection scheme in file %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

Individual* Selection::do_replication(Individual* parent, int32_t index) {
  //Individual* new_indiv = NULL;
  // ===========================================================================
  // Copy parent
  // ===========================================================================
  #ifdef __NO_X
    #ifndef __REGUL
      Individual* new_indiv = new Individual(parent, index, parent->mut_prng(), parent->stoch_prng() );
    #else
      Individual_R* new_indiv = new Individual_R(dynamic_cast<Individual_R*>(parent), index, parent->mut_prng(), parent->stoch_prng() );
    #endif
  #elif defined __X11
    #ifndef __REGUL
      Individual_X11* new_indiv = new Individual_X11(dynamic_cast<Individual_X11 *>(parent), index, parent->mut_prng(), parent->stoch_prng() );
    #else
      Individual_R_X11* new_indiv = new Individual_R_X11(dynamic_cast<Individual_R_X11*>(parent), index, parent->mut_prng(), parent->stoch_prng() );
    #endif
  #endif

  // Notify observers that a new individual was created from <parent>
  {
    Individual* msg[2] = {new_indiv, parent};
    notifyObservers(NEW_INDIV, msg);
  }

  // Perform transfer, rearrangements and mutations
  if (not new_indiv->allow_plasmids()) {
    const GeneticUnit* chromosome = &new_indiv->genetic_unit_list().front();

    chromosome->dna()->perform_mutations(parent->id());
  }
  else { // For each GU, apply mutations
    // Randomly determine the order in which the GUs will undergo mutations
    bool inverse_order = (prng_->random((int32_t) 2) < 0.5);

    if (not inverse_order) { // Apply mutations in normal GU order
      for (const auto& gen_unit: new_indiv->genetic_unit_list()) {
        gen_unit.dna()->perform_mutations(parent->id());
      }
    }
    else { // Apply mutations in inverse GU order
      const auto& gul = new_indiv->genetic_unit_list();
      for (auto gen_unit = gul.crbegin(); gen_unit != gul.crend(); ++gen_unit) {
        gen_unit->dna()->perform_mutations(parent->id());
      }
    }
  }

  return new_indiv;
}

Individual* Selection::do_replication(Individual* parent, int32_t index,
                                      int16_t x, int16_t y) {
  Individual* new_indiv = do_replication(parent, index);

  // Set the new individual's location on the grid
  exp_m_->world()->PlaceIndiv(new_indiv, x, y);

  // Evaluate new individual
  new_indiv->Evaluate();

  // Compute statistics
  new_indiv->compute_statistical_data();

  // Tell observers the replication is finished
  new_indiv->notifyObservers(END_REPLICATION, nullptr);

  return new_indiv;
}

Individual *Selection::do_local_competition (int16_t x, int16_t y) {
  // This function uses the array prob_reprod_ when selection scheme is
  // RANK_LINEAR, RANK_EXPONENTIAL, or FITTEST. For these selection schemes,
  // the function compute_local_prob_reprod (creating the array prob_reprod_)
  // must have been called before.
  // When selection scheme is FITNESS_PROPORTIONATE, this function only uses
  // the fitness values

  World* world = exp_m_->world();

  int16_t neighborhood_size = 9;
  int16_t grid_width  = world->width();
  int16_t grid_height = world->height();
  int16_t cur_x;
  int16_t cur_y;

  // Build a temporary local array of fitness values
  double *  local_fit_array   = new double[neighborhood_size];
  double *  sort_fit_array    = new double[neighborhood_size];
  int16_t * initial_location  = new int16_t[neighborhood_size];
  double *  probs             = new double[neighborhood_size];
  int16_t   count             = 0;
  double    sum_local_fit     = 0.0;

  for (int8_t i = -1 ; i < 2 ; i++) {
    for (int8_t j = -1 ; j < 2 ; j++) {
      cur_x = (x + i + grid_width)  % grid_width;
      cur_y = (y + j + grid_height) % grid_height;
      local_fit_array[count]  = world->indiv_at(cur_x, cur_y)->fitness();
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
  switch (selection_scheme_) {
    case RANK_LINEAR :
    case RANK_EXPONENTIAL :
    case FITTEST : {
      assert(prob_reprod_);
      // First we sort the local fitness values using bubble sort :
      // we sort by increasing order, so the first element will have the worst fitness.
      bool swaped = true;
      int16_t loop_length = 8;
      double  tmp_holder;
      int16_t tmp_holder2;
      while (swaped == true) {
        swaped = false;
        for (int16_t i = 0 ; i < loop_length ; i++) {
          //if the first is higher than the second,  exchange them
          if (sort_fit_array[i] > sort_fit_array[i+1]) {
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
      for (int16_t i = 0 ; i < neighborhood_size ; i++) {
        probs[initial_location[i]] = prob_reprod_[i];
      }

      break;
    }
    // Fitness proportionate selection
    case FITNESS_PROPORTIONATE : {
      for(int16_t i = 0 ; i < neighborhood_size ; i++) {
        probs[i] = local_fit_array[i]/sum_local_fit;
      }

      break;
    }
    default : {
      printf("ERROR, invalid selection scheme in file %s:%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }


  // pick one organism to reproduce, based on probs[] calculated above, using roulette selection
  int8_t found_org = prng_->roulette_random(probs, 9);

  int16_t x_offset = (found_org / 3) - 1;
  int16_t y_offset = (found_org % 3) - 1;

  delete [] local_fit_array;
  delete [] sort_fit_array;
  delete [] initial_location;
  delete [] probs;

  return world->indiv_at((x+x_offset+grid_width)  % grid_width,
                             (y+y_offset+grid_height) % grid_height);
}

#ifdef __TBB
Individual* Selection::do_replication_tbb(Individual* parent, int32_t index, int16_t x /*= -1 */, int16_t y /*= -1 */ )
{
  // ===========================================================================
  //  1) Copy parent
  // ===========================================================================
#ifdef __NO_X
    #ifndef __REGUL
      Individual* new_indiv = new Individual(parent, index, parent->get_mut_prng(), parent->get_stoch_prng() );
    #else
      Individual_R* new_indiv = new Individual_R(dynamic_cast<Individual_R*>(parent), index, parent->get_mut_prng(), parent->get_stoch_prng() );
    #endif
  #elif defined __X11
#ifndef __REGUL
      Individual_X11* new_indiv = new Individual_X11(dynamic_cast<Individual_X11 *>(parent), index, parent->get_mut_prng(), parent->get_stoch_prng() );
    #else
  Individual_R_X11* new_indiv = new Individual_R_X11(dynamic_cast<Individual_R_X11*>(parent), index, parent->get_mut_prng(), parent->get_stoch_prng() );
#endif
#endif

  // Notify observers that a new individual was created from <parent>
  {
    Individual* msg[2] = {new_indiv, parent};
    notifyObservers(NEW_INDIV, msg);
  }

  // Set the new individual's location on the grid
  _exp_m->world()->PlaceIndiv(new_indiv, x, y);

  // Perform transfer, rearrangements and mutations
  if (not new_indiv->get_allow_plasmids())
  {
    const GeneticUnit* chromosome = &new_indiv->get_genetic_unit_list().front();

    chromosome->get_dna()->perform_mutations(parent->get_id());
  }
  else
  { // For each GU, apply mutations
    // Randomly determine the order in which the GUs will undergo mutations
    bool inverse_order = (new_indiv->get_mut_prng()->random((int32_t) 2) < 0.5);

    if (not inverse_order) { // Apply mutations in normal GU order
      for (const auto& gen_unit: new_indiv->get_genetic_unit_list()) {
        gen_unit.get_dna()->perform_mutations(parent->get_id() );
      }
    }
    else { // Apply mutations in inverse GU order
      const auto& gul = new_indiv->get_genetic_unit_list();
      for (auto gen_unit = gul.crbegin(); gen_unit != gul.crend(); ++gen_unit) {
        gen_unit->get_dna()->perform_mutations(parent->get_id());
      }
    }
  }

  // Evaluate new individual
  new_indiv->Evaluate();

  // Compute statistics
  new_indiv->compute_statistical_data();

  // Tell observers the replication is finished
  new_indiv->notifyObservers(END_REPLICATION, nullptr);

  return new_indiv;
}
#endif

// =================================================================
//                          Non inline accessors
// =================================================================
} // namespace aevol
