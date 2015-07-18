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
#include <err.h>
#include <cerrno>
#include <sys/stat.h>

#include <list>
#include <iostream>

#include <zlib.h>

#include "ExpManager.h"
#include "Individual.h"


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
  _exp_s = new ExpSetup(this);

  // ------------------------------------------------------------------- World
  world_ = nullptr;

  // ---------------------------------------------------------- Output manager
  _output_m = new OutputManager(this);


  // -------------------------------- Timestep up to which we want to simulate
  t_end = 0;


  // ------------------------------------------------------------- Quit signal
  _quit_signal_received = false;
}

// ===========================================================================
//                                  Destructor
// ===========================================================================
ExpManager::~ExpManager()
{
  delete _exp_s;
  delete _output_m;
  delete world_;
}

// ===========================================================================
//                                   Algorithms
// ===========================================================================
// void ExpManager::foreach_indiv(void (*processor)(Individual& indiv)) const
// {
//   for (const auto& indiv: get_indivs_std())
//     processor(*indiv);
// }

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void ExpManager::InitializeWorld(int16_t grid_width,
                                     int16_t grid_height,
                                     std::shared_ptr<JumpingMT> prng,
                                     const Habitat& habitat,
                                     bool share_phenotypic_target)
{
  world_ = new World();
  world_->InitGrid(grid_width, grid_height, habitat, share_phenotypic_target);
  world_->set_prng(prng);
}

/*!
  \brief Save the experiment
*/
void ExpManager::Save() const
{
  WriteSetupFiles();
  _output_m->write_current_generation_outputs();
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
  \see WriteDynamicFiles(void)
  \see save_copy(char* dir, int64_t t)
*/
void ExpManager::WriteSetupFiles() const
{
  // 1) Create missing directories
  create_missing_directories();

  // 2) Open setup files (experimental setup and output profile)
  gzFile exp_s_file, out_p_file;
  open_setup_files(exp_s_file, out_p_file, AeTime::get_time(), "w");

  // 4) Write setup data
  _exp_s->write_setup_file(exp_s_file);
  _output_m->WriteSetupFile(out_p_file);

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
  \see WriteSetupFiles(void)
  \see save_copy(char* dir, int64_t t)
*/
void ExpManager::WriteDynamicFiles(void) const
{
  // Create missing directories
  create_missing_directories();

  // Open backup files
  gzFile sel_file, world_file;
  open_backup_files(sel_file, world_file, AeTime::get_time(), "w");

  // Save experiment
  get_sel()->save(sel_file);
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
  \see WriteSetupFiles(void)
  \see WriteDynamicFiles(void)
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
  _exp_s->write_setup_file(exp_s_file);
  _output_m->WriteSetupFile(out_p_file);

  // Write the state of selection and world into the backups
  get_sel()->save(sel_file);
  world_->save(world_file);

  // Close setup and backup files
  close_setup_files(exp_s_file, out_p_file);
  close_backup_files(sel_file, world_file);

  // Copy stats
  _output_m->CopyStats(dir, time);

  // Write last gener file
  _output_m->WriteLastGenerFile(dir);
}

void ExpManager::step_to_next_generation() {
  // TODO <david.parsons@inria.fr> Apply phenotypic target variation and noise
  world_->ApplyHabitatVariation();

  // Take a step in time
  AeTime::plusplus();

  // Create the corresponding new generation
  _exp_s->step_to_next_generation();

  // Write statistical data and store phylogenetic data (tree)
  _output_m->write_current_generation_outputs();
}

/*!
  \brief Load an experiment with the provided files
 */
// TODO <david.parsons@inria.fr> check verbose (what doas it do ?, is it consistent ?)
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
  _exp_s->load(exp_s_file, exp_backup_file, verbose);
  printf(" OK\n");

  // ---------------------------------------------------------- Retrieve world
  printf("  Loading world...");
  fflush(stdout);
  world_ = new World();
  world_->load(world_file, this);
  printf(" OK\n");

  // --------------------------------------------- Retrieve output profile data
  printf("  Loading output profile...");
  fflush(stdout);
  _output_m->load(out_p_file, verbose, to_be_run);
  printf(" OK\n");

  // -------------------------------------------- Link world and output profile
  if (get_record_tree()) {
    get_sel()->addObserver(get_tree(), NEW_INDIV);
    for (auto indiv : world_->get_indivs())
      indiv->addObserver(
          get_tree()->get_report_by_index(AeTime::get_time(), indiv->get_id()),
          END_REPLICATION);
    get_sel()->addObserver(get_tree(), END_GENERATION);
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

/**
 * Run the simulation
 */
void ExpManager::run_evolution(void)
{
  // We are running a simulation.
  // Save the setup files to keep track of the setup history
  WriteSetupFiles();

  // For each generation
  while (true) { // termination condition is into the loop
    printf("============================== %" PRId64 " ==============================\n",
           AeTime::get_time());
    printf("  Best individual's distance to target (metabolic) : %f\n",
           get_best_indiv()->get_dist_to_target_by_feature(METABOLISM));

    if (AeTime::get_time() >= t_end or quit_signal_received())
      break;

#ifdef __X11
    display();
#endif

    // Take one step in the evolutionary loop
    step_to_next_generation();
  }

  _output_m->flush();
  printf("================================================================\n");
  printf("  The run is finished. \n");
  printf("  Printing the final best individual into " BEST_LAST_ORG_FNAME "\n");
  FILE* org_file = fopen(BEST_LAST_ORG_FNAME, "w");
  fputs(get_best_indiv()->get_genetic_unit_sequence(0), org_file);
  fclose(org_file);
}

void ExpManager::update_best(void)
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
Individual* ExpManager::get_indiv_by_id(int32_t id) const {
  return world_->get_indiv_by_id(id);
}
} // namespace aevol
