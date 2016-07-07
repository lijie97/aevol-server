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
#include <err.h>
#include <cerrno>
#include <sys/stat.h>

#include <list>
#include <iostream>

#include <zlib.h>

#include "ExpManager.h"
#include "Individual.h"

#include "raevol/cuda_struct.h"

#ifdef __REGUL
#include "raevol/Individual_R.h"
#endif
#include<chrono>

#include <iostream>
#include <unordered_map>
using namespace std;
using namespace std::chrono;

#ifdef __TRACING__
#include "ae_logger.h"


unordered_map<int,unordered_multiset<string>> ae_logger::logMap;
string ae_logger::logFile = "logger_csv.log";
mutex ae_logger::loggerMtx;
#endif

using std::cout;
using std::endl;


namespace aevol {
//##############################################################################
//                                                                             #
//                             Class ExpManager                            #
//                                                                             #
//##############################################################################

// ===========================================================================
//                        Definition of static attributes
// ===========================================================================

// ===========================================================================
//                                 Constructors
// ===========================================================================
ExpManager::ExpManager()
{
  // ------------------------------------------------------ Experimental setup
  exp_s_ = new ExpSetup(this);

  // ------------------------------------------------------------------- World
  world_ = nullptr;

  // ---------------------------------------------------------- Output manager
  output_m_ = new OutputManager(this);


  // -------------------------------- Timestep up to which we want to simulate
  t_end_ = 0;


  // ------------------------------------------------------------- Quit signal
  quit_signal_received_ = false;
}

// ===========================================================================
//                                  Destructor
// ===========================================================================
ExpManager::~ExpManager() noexcept
{
  delete exp_s_;
  delete output_m_;
  delete world_;
}

// ===========================================================================
//                                   Algorithms
// ===========================================================================

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void ExpManager::InitializeWorld(int16_t grid_width,
                                     int16_t grid_height,
                                     std::shared_ptr<JumpingMT> prng,
                                     std::shared_ptr<JumpingMT> mut_prng,
                                     std::shared_ptr<JumpingMT> stoch_prng,
                                     Habitat& habitat,
                                     bool share_phenotypic_target)
{
  world_ = new World();
  world_->set_prng(prng);
  world_->set_mut_prng(mut_prng);
  world_->set_stoch_prng(stoch_prng);
  world_->InitGrid(grid_width, grid_height, habitat, share_phenotypic_target);
}

/*!
  \brief Save the experiment
*/
void ExpManager::Save() const
{
  WriteSetupFiles();
  output_m_->write_current_generation_outputs();
}

/*!
  \brief Save all the static data (the data that is constant throughout
  a simulation)

  Save
    * The experimental setup
    * The output profile

  These are written in the current directory in gzipped binary files
  as well as plain text files
  This data should be saved once and for all when the experiment is created.
  There is no need to save this data in the middle of a simulation since it
  is constant throughout the experiment.

  \see load(int64_t t0,
            char* exp_setup_file_name,
            char* out_prof_file_name,
            char* pop_file_name,
            char* sel_file_name,
            char* world_file_name,
            bool verbose)
  \see WriteDynamicFiles()
  \see save_copy(char* dir, int64_t t)
*/
void ExpManager::WriteSetupFiles() const
{
  // 1) Create missing directories
  create_missing_directories();

  // 2) Open setup files (experimental setup and output profile)
  gzFile exp_s_file, out_p_file;
  open_setup_files(exp_s_file, out_p_file, AeTime::time(), "w");

  // 4) Write setup data
  exp_s_->write_setup_file(exp_s_file);
  output_m_->WriteSetupFile(out_p_file);

  // 5) Close setup files
  close_setup_files(exp_s_file, out_p_file);
}

/*!
  \brief Save all the dynamic data (the data that changes over the simulation)

  Save the state of
    * The population
    * The world
    * The PRNG (random generator) used for the selection

  These are written in the current directory in gzipped binary files
  whose names are appended with the current generation number.

  WARNING: The experimental setup and output profile are *not* saved. These
  should be saved once and for all when the experiment is created.

  \see load(int64_t t0,
            char* exp_setup_file_name,
            char* out_prof_file_name,
            char* pop_file_name,
            char* sel_file_name,
            char* world_file_name,
            bool verbose)
  \see WriteSetupFiles()
  \see save_copy(char* dir, int64_t t)
*/
void ExpManager::WriteDynamicFiles() const
{
  // Create missing directories
  create_missing_directories();

  // Open backup files
  gzFile sel_file, world_file;
  open_backup_files(sel_file, world_file, AeTime::time(), "w");

  // Save experiment
  sel()->save(sel_file);
  world_->save(world_file);

  // Close backup files
  close_backup_files(sel_file, world_file);
}

/*!
  \brief Saves a complete copy of the experiment at the provided location.

  Save both the setup (constant) and the backup (dynamic) files to the <dir>
  directory. Dynamic file names will be appended with the <num_gener>

  Static files (saved both as gzipped binary and plain text)
    * The experimental setup
    * The output profile

  Dynamic files (saved as gzipped binary)
    * The world
    * The PRNG (random generator) used for the selection

  \see load(int32_t first_gener,
             char* exp_setup_file_name,
             char* out_prof_file_name,
             char* pop_file_name,
             char* sel_file_name,
             char* world_file_name,
             bool verbose)
  \see WriteSetupFiles()
  \see WriteDynamicFiles()
*/
void ExpManager::save_copy(char* dir, int64_t time) const
{
  // Set time to time
  AeTime::set_time(time);

  // Create missing directories
  create_missing_directories(dir);

  // Open setup files and backup files
  gzFile exp_s_file, out_p_file, sel_file, world_file;
  open_setup_files(exp_s_file, out_p_file, time, "w", dir);
  open_backup_files(sel_file, world_file, time, "w", dir);

  // Write setup data
  exp_s_->write_setup_file(exp_s_file);
  output_m_->WriteSetupFile(out_p_file);

  // Write the state of selection and world into the backups
  sel()->save(sel_file);
  world_->save(world_file);

  // Close setup and backup files
  close_setup_files(exp_s_file, out_p_file);
  close_backup_files(sel_file, world_file);

  // Copy stats
  output_m_->CopyStats(dir, time);

  // Write last gener file
  output_m_->WriteLastGenerFile(dir);
}

void ExpManager::step_to_next_generation() {
  // TODO <david.parsons@inria.fr> Apply phenotypic target variation and noise
  world_->ApplyHabitatVariation();

  // Take a step in time
  AeTime::plusplus();

  // Create the corresponding new generation
  exp_s_->step_to_next_generation();

  // Write statistical data and store phylogenetic data (tree)
  output_m_->write_current_generation_outputs();
}

/*!
  \brief Load an experiment with the provided files
 */
// TODO <david.parsons@inria.fr> check verbose (what does it do ?, is it consistent ?)
void ExpManager::load(gzFile& exp_s_file,
                          gzFile& exp_backup_file,
                          gzFile& world_file,
                          gzFile& out_p_file,
                          bool verbose,
                          bool to_be_run /*  = true */)
{
  // ---------------------------------------- Retrieve experimental setup data
  printf("  Loading experimental setup...");
  fflush(stdout);
  exp_s_->load(exp_s_file, exp_backup_file, verbose);
  printf(" OK\n");

  if (FuzzyFactory::fuzzyFactory == NULL)
    FuzzyFactory::fuzzyFactory = new FuzzyFactory(exp_s_);


  // ---------------------------------------------------------- Retrieve world
  printf("  Loading world...");
  fflush(stdout);
  world_ = new World();
  world_->load(world_file, this);
  printf(" OK\n");

  sel()->set_unique_id(  grid_height()*grid_width()+1 );

  // --------------------------------------------- Retrieve output profile data
  printf("  Loading output profile...");
  fflush(stdout);
  output_m_->load(out_p_file, verbose, to_be_run);
  printf(" OK\n");


  // -------------------------------------------- Link world and output profile
  if (record_tree()) {
    sel()->addObserver(tree(), NEW_INDIV);
    for (int16_t x = 0 ; x < grid_width() ; x++) {
      for (int16_t y = 0; y < grid_height(); y++) {
        world_->indiv_at(x,y)->addObserver(
            tree(),
            END_REPLICATION);
      }
    }
    sel()->addObserver(tree(), END_GENERATION);
  }

  // --------------------------------------------------- Recompute unsaved data
  world_->evaluate_individuals();
}


/*!
  \brief Load an experiment with default files from a given directory
 */
// TODO <david.parsons@inria.fr> check verbose (what does it do ?, is it consistent ?)
void ExpManager::load(const char* dir,
    int64_t t0, bool verbose, bool to_be_run /*  = true */)
{
  AeTime::set_time(t0);

  // -------------------------------------------------------------------------
  // Open setup files and backup files
  // -------------------------------------------------------------------------
  gzFile exp_s_file, out_p_file;
  gzFile exp_backup_file, world_file;
  open_setup_files(exp_s_file, out_p_file, t0, "r", dir);
  open_backup_files(exp_backup_file, world_file, t0, "r", dir);


  // -------------------------------------------------------------------------
  // Load data from backup and parameter files
  // -------------------------------------------------------------------------
  load(exp_s_file, exp_backup_file,
       world_file, out_p_file, verbose, to_be_run);


  // -------------------------------------------------------------------------
  // Close setup and backup files
  // -------------------------------------------------------------------------
  close_setup_files(exp_s_file, out_p_file);
  close_backup_files(exp_backup_file, world_file);


  if (FuzzyFactory::fuzzyFactory == NULL)
    FuzzyFactory::fuzzyFactory = new FuzzyFactory(exp_s_);
  printf("Factory flavor %d : %d\n",exp_s_->get_fuzzy_flavor(),FuzzyFactory::fuzzyFactory->get_fuzzy_flavor());
}


/**
 * \brief Load an experiment with the provided constitutive files
 */
// TODO <david.parsons@inria.fr> check verbose (what does it do ?, is it consistent ?)
void ExpManager::load(int64_t t0,
                          char* exp_setup_file_name,
                          char* exp_backup_file_name,
                          char* world_file_name,
                          char* out_prof_file_name,
                          bool verbose /*= false*/,
                          bool to_be_run /*= true*/)
{
  AeTime::set_time(t0);

  // ---------------------------------------------------------------------------
  // Open files and check them
  // ---------------------------------------------------------------------------
  gzFile exp_setup_file = gzopen(exp_setup_file_name, "r");
  gzFile out_prof_file = gzopen(out_prof_file_name, "r");
  gzFile exp_backup_file = gzopen(exp_backup_file_name, "r");
  gzFile world_file  = gzopen(world_file_name, "r");


  if (exp_setup_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, exp_setup_file_name);
    exit(EXIT_FAILURE);
  }
  if (out_prof_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, out_prof_file_name);
    exit(EXIT_FAILURE);
  }
  if (exp_backup_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, exp_backup_file_name);
    exit(EXIT_FAILURE);
  }
  if (world_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, world_file_name);
    exit(EXIT_FAILURE);
  }


  // ---------------------------------------------------------------------------
  // Actually load data
  // ---------------------------------------------------------------------------
  load(exp_setup_file, exp_backup_file, world_file, out_prof_file,
       verbose, to_be_run);


  // ---------------------------------------------------------------------------
  // Close setup and backup files
  // ---------------------------------------------------------------------------
  gzclose(exp_setup_file);
  gzclose(exp_backup_file);
  gzclose(world_file);
  gzclose(out_prof_file);
}



//int mutator = 0;

/**
 * Run the simulation
 */
void ExpManager::run_evolution() {
  // We are running a simulation.
  // Save the setup files to keep track of the setup history
  WriteSetupFiles();
#ifdef __TRACING__
  ae_logger::init("logger_csv.log");

  printf("Launching TRACING...\n");
#else
  printf("Launching NOT TRACING...\n");
#endif

#ifdef __TRACING__
  high_resolution_clock::time_point t_t1 = high_resolution_clock::now();
  high_resolution_clock::time_point t_t2,t1,t2;
#endif

  int max_protein = 0;
  int max_rna = 0;
  int max_influence = 0;

  const Habitat_R& hab = dynamic_cast<const Habitat_R&>(world()->grid(0,
                                                                      0)->habitat());

  int nb_signals = hab.signals().size();
  int life_time = exp_s()->get_nb_indiv_age();
  int nb_eval_ = exp_s()->get_nb_degradation_step();
  float selection_pressure = sel()->selection_pressure();
  World* world;
  int16_t grid_width;
  int16_t grid_height;
  int32_t pop_size;

  GridCell*** pop_grid;

  // create a temporary grid to store the reproducers

  int16_t x, y;
  int8_t what;
  high_resolution_clock::time_point t1, t2;

  std::unordered_map<unsigned long long, Individual*> unique_individual;

  std::list<Individual*> new_generation;
  std::vector<Individual*> to_evaluate;
  std::list<Individual*> old_generation;
  Individual*** reproducers = NULL;

    // For each generation
    while (AeTime::time() < t_end_ or !quit_signal_received()) { // termination condition is into the loop
#pragma omp parallel
      {
#pragma omp single
      {
      unique_individual.clear();

      new_generation.clear();
      to_evaluate.clear();
      old_generation.clear();

// Create proxies

//#pragma omp single
//      {

        printf(
            "============================== %" PRId64 " ==============================\n",
            AeTime::time());
        printf("  Best individual's distance to target (metabolic) : %f\n",
               best_indiv()->dist_to_target_by_feature(METABOLISM));


        int16_t nb_activators = 0;
        int16_t nb_operators = 0;

        Individual_R* test = dynamic_cast<Individual_R*>(best_indiv());

        //test->init_indiv(dynamic_cast<const Habitat_R&>(test->habitat()));

        int nb_protein = 0;
        for (const auto& rnax: test->_rna_list_coding) {
          Rna_R* rna = (Rna_R*) rnax;
          for (unsigned int i = 0; i < rna->nb_influences(); i++) {
            nb_protein++;
            if (rna->_enhancing_coef_list[i] > 0) {
              nb_activators++;
            }

            if (rna->_operating_coef_list[i] > 0) {
              nb_operators++;
            }
          }
        }

        printf("  Proteins %ld (%d) - RNA %ld - Link A %d - I %d\n",
               test->protein_list().size(), nb_protein,
               test->_rna_list_coding.size(),
               nb_activators, nb_operators);

        for (auto indiv : indivs()) {
          Individual_R* indiv_r = dynamic_cast<Individual_R*>(indiv);

          if (indiv_r->protein_list().size() > max_protein)
            max_protein = indiv_r->protein_list().size();

          if (indiv_r->_rna_list_coding.size() > max_rna)
            max_rna = indiv_r->_rna_list_coding.size();

          for (auto rna : indiv_r->_rna_list_coding)
            if (rna->_nb_influences > max_influence)
              max_influence = rna->_nb_influences;
        }
  //    }

      //if ()
      //  break;





//#pragma omp single
//      {
#ifdef __X11
        display();
#endif


#ifdef __TRACING__
        t1 = high_resolution_clock::now();
#endif
        // Take one step in the evolutionary loop
        /** Step to next generation **/

        // TODO <david.parsons@inria.fr> Apply phenotypic target variation and noise
        world_->ApplyHabitatVariation();

        // Take a step in time
        AeTime::plusplus();


      // Create the corresponding new generation
      /** Selection :: step_to_next_generation **/
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

      if (sel()->prng_ == NULL) {
        printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }

      // -------------------------------------------------------------------------------
      // 1) Compute the probability of reproduction of each individual in the population
      // -------------------------------------------------------------------------------
#ifndef FIXED_POPULATION_SIZE
#error this method is not ready for variable population size
      sel()->compute_local_prob_reprod();
#else
      // The function compute_local_prob_reprod creates and fills the array prob_reprod_, which is telling us the probability of being picked for reproduction according to the rank of an individual in its neighboorhood.
      // It is only usefull when selection is rank based. When selection scheme is FITNESS_PROPORTIONATE, we do not need to call it.
      // It shoud only be called once in the simulation and not at each generation. So if prob_reprod_ already exists we do not need to call it.
      if ((sel()->selection_scheme() != FITNESS_PROPORTIONATE) &&
          (sel()->prob_reprod() == NULL)) {
        sel()->compute_local_prob_reprod();
      }
#endif

      if (sel()->prng_ == NULL) {
        printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
  //    }


//#pragma omp single
      //{
        // Create proxies
        world = this->world();
        grid_width = world->width();
        grid_height = world->height();
        pop_size = grid_width * grid_height;

        pop_grid = grid();

        // create a temporary grid to store the reproducers
        reproducers = new Individual** [grid_width];


        for (int16_t i = 0; i < grid_width; i++) {
          reproducers[i] = new Individual* [grid_height];
        }
      }




      if (sel()->global) {
#pragma omp single
        {
          delete[] sel()->prob_reprod_;

          sel()->prob_reprod_ = new double[pop_size];

          double* fitnesses = new double[pop_size];
          double sum = 0;

          size_t i = 0;
          for (const auto& indiv: indivs()) {
            fitnesses[i] = indiv->fitness();
            sum += fitnesses[i];
            ++i;
          }

          for (int32_t i = 0; i < pop_size; i++) {
            sel()->prob_reprod_[i] = fitnesses[i] / sum;
          }

          delete[] fitnesses;

          int32_t* nb_offsprings = new int32_t[pop_size];
          sel()->prng_->multinomial_drawing(nb_offsprings, sel()->prob_reprod_,
                                            pop_size,
                                            pop_size);

          int index = 0;
          i = 0;

          for (const auto& indiv: indivs()) {
            for (int32_t j = 0; j < nb_offsprings[i]; j++) {
              x = index / grid_height;
              y = index % grid_height;

              reproducers[x][y] = indiv;

              index++;
            }
            i++;
          }
        }
        //printf("index %d\n",index);
      } else {
        // Do local competitions
#pragma omp for schedule(dynamic) private(x,y)
        for (int32_t index = 0; index < grid_width * grid_height; index++) {
          x = index / grid_height;
          y = index % grid_height;
          reproducers[x][y] = sel()->do_local_competition(x, y);
        }

#pragma omp barrier
      }

      //}

#pragma omp single
      {
        // TODO : Why is that not *after* the creation of the new population ?
        // Add the compound secreted by the individuals
        if (with_secretion()) {
          for (int16_t x = 0; x < grid_width; x++) {
            for (int16_t y = 0; y < grid_height; y++) {
              pop_grid[x][y]->set_compound_amount(
                  pop_grid[x][y]->compound_amount() +
                  pop_grid[x][y]->individual()->fitness_by_feature(SECRETION));
            }
          }

          // Diffusion and degradation of compound in the habitat
          world->update_secretion_grid();
        }

      }
        // Create the new generation

#pragma omp single
      {
        old_generation = indivs();
        //mutator = 0;
        for (auto indiv : old_generation) {
          indiv->number_of_clones_ = 0;
          unique_individual[indiv->id()] = indiv;
          (&indiv->genetic_unit_list().front())->dna()->set_hasMutate(false);
        }

      }


#pragma omp for schedule(dynamic) private(x,y,what)
        for (int32_t index = 0; index < grid_width * grid_height; index++) {
          x = index / grid_height;
          y = index % grid_height;
          sel()->do_replication(reproducers[x][y],
                                x * grid_height + y + pop_size * AeTime::time(),
                                what, x, y);
          if (what == 1 || what == 2) {
#pragma omp critical(updateindiv)
            {
              to_evaluate.push_back(pop_grid[x][y]->individual());
            }
          }
        }
#pragma omp barrier

#pragma omp single
        {
          t1 = high_resolution_clock::now();
        }

#pragma omp for schedule(dynamic)
        for (int i = 0; i < to_evaluate.size(); i++) {
          if ((dynamic_cast<PhenotypicTargetHandler_R*>(&to_evaluate[i]->grid_cell()->habitat().
              phenotypic_target_handler_nonconst())->hasChanged()) ||
              !to_evaluate[i]->evaluated_)
            sel()->run_life(dynamic_cast<Individual_R*>(to_evaluate[i]));
        }
#pragma omp barrier


#pragma omp single
        {
          for (int32_t index = 0; index < grid_width * grid_height; index++) {
            x = index / grid_height;
            y = index % grid_height;

            EndReplicationEvent* eindiv = new EndReplicationEvent(
                world->indiv_at(x, y), x, y);
            // Tell observers the replication is finished
            world->indiv_at(x, y)->notifyObservers(END_REPLICATION, eindiv);
          }

        t2 = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            t2 - t1).count();
        cout << "TIMER," << AeTime::time() << ",OLD," << duration << endl;

        for (int16_t x = 0; x < grid_width; x++)
          for (int16_t y = 0; y < grid_height; y++)
            new_generation.push_back(pop_grid[x][y]->individual());

        // delete the temporary grid and the parental generation
        for (int16_t x = 0; x < grid_width; x++) {
          for (int16_t y = 0; y < grid_height; y++) {
            reproducers[x][y] = nullptr;
          }
          delete[] reproducers[x];
        }
        delete[] reproducers;
        // Compute the rank of each individual
        new_generation.sort([](Individual* lhs, Individual* rhs) {
            return lhs->fitness() < rhs->fitness();
        });

        int rank = 1;
        for (Individual* indiv : new_generation) {
          indiv->set_rank(rank++);
        }

        // randomly migrate some organisms, if necessary
        world->MixIndivs();

        sel()->PerformPlasmidTransfers();

        // Update the best individual
        update_best();

        // Notify observers of the end of the generation
        sel()->notifyObservers(END_GENERATION);

        int number_of_clones = 0;
        for (auto iterator = unique_individual.begin();
             iterator != unique_individual.end(); iterator++) {
          if (iterator->second->number_of_clones_ == 0) {
            delete iterator->second;
          }
        }

        /** END of Selection :: step_to_next_generation **/

        // Write statistical data and store phylogenetic data (tree)
        output_m_->write_current_generation_outputs();

        /** END of step to next generation **/
#ifdef __TRACING__
        t2 = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        ae_logger::addLog(SELECTION,duration);
        ae_logger::flush(AeTime::get_time());
#endif
        }
      }
      }
  //  }
//  }
#ifdef __TRACING__
  t_t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t_t2 - t_t1 ).count();
  ae_logger::addLog(TOTAL,duration);
  ae_logger::flush(AeTime::get_time());
#endif


  output_m_->flush();
  printf("================================================================\n");
  printf("  The run is finished. \n");
  printf("  Printing the final best individual into " BEST_LAST_ORG_FNAME "\n");
  FILE* org_file = fopen(BEST_LAST_ORG_FNAME, "w");
  fputs(best_indiv()->genetic_unit_sequence(0), org_file);
  fclose(org_file);

}

void ExpManager::update_best()
{
  world_->update_best();
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================
void ExpManager::create_missing_directories(const char* dir /*= "."*/) const
{
  char cur_dir_name[255];
  int status;

  // Base directory
  status = mkdir(dir, 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }

  // Experimental setup
  sprintf(cur_dir_name, "%s/" EXP_S_DIR, dir);
  status = mkdir(cur_dir_name, 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }
  // Output profile
  sprintf(cur_dir_name, "%s/" OUT_P_DIR, dir);
  status = mkdir(cur_dir_name, 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }
  // World
  sprintf(cur_dir_name, "%s/" WORLD_DIR, dir);
  status = mkdir(cur_dir_name, 0755);
  if (status == -1 && errno != EEXIST)
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }
  // Stats
  sprintf(cur_dir_name, "%s/" STATS_DIR, dir);
  status = mkdir(cur_dir_name, 0755);
  if (status == -1 && errno != EEXIST)
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }
}

void ExpManager::open_backup_files(gzFile& exp_backup_file,
                                       gzFile& world_file,
                                       int64_t t,
                                       const char mode[3],
                                       const char* dir /*= "."*/) const
{
  assert(strcmp(mode, "w") == 0 or strcmp(mode, "r") == 0);

  // -------------------------------------------------------------------------
  // Generate backup file names for mandatory files.
  // -------------------------------------------------------------------------
  char exp_backup_file_name[255];
  char world_file_name[255];

  sprintf(exp_backup_file_name, "%s/" EXP_S_FNAME_FORMAT, dir, t);
  sprintf(world_file_name, "%s/" WORLD_FNAME_FORMAT, dir, t);


  // -------------------------------------------------------------------------
  // Open backup files
  // -------------------------------------------------------------------------
  exp_backup_file = gzopen(exp_backup_file_name, mode);
  world_file = gzopen(world_file_name, mode);


  // -------------------------------------------------------------------------
  // Check that files were correctly opened
  // -------------------------------------------------------------------------
  if (exp_backup_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, exp_backup_file_name);
    exit(EXIT_FAILURE);
  }
  if (world_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, world_file_name);
    exit(EXIT_FAILURE);
  }
}

void ExpManager::close_backup_files(gzFile& exp_backup_file,
                                        gzFile& world_file) const
{
  gzclose(exp_backup_file);
  gzclose(world_file);
}

void ExpManager::open_setup_files(
    gzFile& exp_s_file,
    gzFile& out_p_file,
    int64_t t,
    const char mode[3],
    const char* dir /*= "."*/) const
{
  // 1) Generate setup file names
  char exp_s_file_name[255];
  char out_p_file_name[255];

  sprintf(exp_s_file_name, "%s/" EXP_S_CONST_FNAME_FORMAT, dir);
  sprintf(out_p_file_name, "%s/" OUT_P_FNAME_FORMAT, dir);

  // 2) Open backup files
  exp_s_file = gzopen(exp_s_file_name, mode);
  out_p_file = gzopen(out_p_file_name, mode);

  // 3) Check that files were correctly opened
  if (exp_s_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, exp_s_file_name);
    exit(EXIT_FAILURE);
  }
  if (out_p_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, out_p_file_name);
    exit(EXIT_FAILURE);
  }
}

void ExpManager::close_setup_files(gzFile& exp_s_file,
                                       gzFile& out_p_file) const
{
  gzclose(exp_s_file);
  gzclose(out_p_file);
}

// ===========================================================================
//                               Non inline accessors
// ===========================================================================
Individual* ExpManager::indiv_by_id(int32_t id) const {
  return world_->indiv_by_id(id);
}

} // namespace aevol
