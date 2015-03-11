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
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "ae_exp_manager.h"
#include "ae_individual.h"



namespace aevol {
//##############################################################################
//                                                                             #
//                             Class ae_exp_manager                            #
//                                                                             #
//##############################################################################

// ===========================================================================
//                        Definition of static attributes
// ===========================================================================

// ===========================================================================
//                                 Constructors
// ===========================================================================
ae_exp_manager::ae_exp_manager(void)
{
  // ------------------------------------------------------ Experimental setup
  _exp_s = new ae_exp_setup(this);

  // ------------------------------------------------------------- Environment
  _env = new Environment();

  // ------------------------------------------------------------------- World
  world_ = NULL;

  // ---------------------------------------------------------- Output manager
  _output_m = new ae_output_manager(this);


  // -------------------------------- Timestep up to which we want to simulate
  t_end = 0;


  // ------------------------------------------------------------- Quit signal
  _quit_signal_received = false;
}

// ===========================================================================
//                                  Destructor
// ===========================================================================
ae_exp_manager::~ae_exp_manager(void)
{
  delete _exp_s;
  delete _output_m;
  delete _env;
  delete world_;
}

// ===========================================================================
//                                   Algorithms
// ===========================================================================
// void ae_exp_manager::foreach_indiv(void (*processor)(ae_individual& indiv)) const
// {
//   for (const auto& indiv: get_indivs_std())
//     processor(*indiv);
// }

// ===========================================================================
//                                 Public Methods
// ===========================================================================
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
            char* env_file_name,
            char* pop_file_name,
            char* sel_file_name,
            char* world_file_name,
            bool verbose)
  \see save(void)
  \see save_copy(char* dir, int64_t t)
*/
void ae_exp_manager::write_setup_files(void)
{
  // 1) Create missing directories
  create_missing_directories();

  // 2) Open setup files (experimental setup and output profile)
  gzFile exp_s_file, out_p_file;
  open_setup_files(exp_s_file, out_p_file, Time::get_time(), "w");

  // 4) Write setup data
  _exp_s->write_setup_file(exp_s_file);
  _output_m->write_setup_file(out_p_file);

  // 5) Close setup files
  close_setup_files(exp_s_file, out_p_file);
}

/*!
  \brief Save all the dynamic data (the data that changes over the simulation)

  Save the state of
    * The population
    * The environment iff it is changing (i.e. it has variation and/or noise)
    * The world
    * The PRNG (random generator) used for the selection

  These are written in the current directory in gzipped binary files
  whose names are appended with the current generation number.

  WARNING: The experimental setup and output profile are *not* saved. These
  should be saved once and for all when the experiment is created.

  \see load(int64_t t0,
            char* exp_setup_file_name,
            char* out_prof_file_name,
            char* env_file_name,
            char* pop_file_name,
            char* sel_file_name,
            char* world_file_name,
            bool verbose)
  \see write_setup_files(void)
  \see save_copy(char* dir, int64_t t)
*/
void ae_exp_manager::save(void) const
{
  // 1) Create missing directories
  create_missing_directories();

  // 2) Open backup files (environment, selection and world)
  gzFile env_file, sel_file, world_file;
  open_backup_files(env_file, sel_file, world_file,
      Time::get_time(), "w");

  // 3) Write the state of the environment and world into the backups
  get_env()->save(env_file);
  get_sel()->save(sel_file);
  world_->save(world_file);

  // 4) Close backup files
  close_backup_files(env_file, sel_file, world_file);
}

/*!
  \brief Saves a complete copy of the experiment at the provided location.

  Save both the setup (constant) and the backup (dynamic) files to the <dir>
  directory. Dynamic file names will be appended with the <num_gener>

  Static files (saved both as gzipped binary and plain text)
    * The experimental setup
    * The output profile

  Dynamic files (saved as gzipped binary)
    * The environment iff it is changing (i.e. it has variation and/or noise)
    * The world
    * The PRNG (random generator) used for the selection

  \see load(int32_t first_gener,
             char* exp_setup_file_name,
             char* out_prof_file_name,
             char* env_file_name,
             char* pop_file_name,
             char* sel_file_name,
             char* world_file_name,
             bool verbose)
  \see write_setup_files(void)
  \see save(void)
*/
void ae_exp_manager::save_copy(char* dir, int32_t num_gener /*= 0*/) const
{
  // 1) Create missing directories
  create_missing_directories(dir);

  // 2) Open setup files (experimental setup and output profile)
  //    and backup files (environment, selection and world)
  gzFile exp_s_file, out_p_file, env_file, sel_file, world_file;
  open_setup_files(exp_s_file, out_p_file, num_gener, "w", dir);
  open_backup_files(env_file, sel_file, world_file, num_gener, "w", dir);

  // 3) Write setup data
  _exp_s->write_setup_file(exp_s_file);
  _output_m->write_setup_file(out_p_file);

  // 4) Write the state of the environment, selection and world into the backups
  get_env()->save(env_file);
  get_sel()->save(sel_file);
  world_->save(world_file);

  // 4) Close setup and backup files
  close_setup_files(exp_s_file, out_p_file);
  close_backup_files(env_file, sel_file, world_file);
}


/*!
  \brief Load an experiment with the provided files
 */
void ae_exp_manager::load(gzFile& env_file,
                          gzFile& exp_s_file,
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

  // --------------------------------------------- Retrieve environmental data
  printf("  Loading environment...");
  fflush(stdout);
  //~ delete _env;
  //~ _env = new ae_environment();
  _env->load(env_file);
  printf(" OK\n");

  // --------------------------------------------- Retrieve output profile data
  printf("  Loading output profile...");
  fflush(stdout);
  _output_m->load(out_p_file, verbose, to_be_run);
  printf(" OK\n");
}


/*!
  \brief Load an experiment with default files from a given directory
 */
void ae_exp_manager::load(const char* dir,
    int64_t t0, bool verbose, bool to_be_run /*  = true */)
{
  Time::set_time(t0);

  // -------------------------------------------------------------------------
  // 1) Open setup files (experimental setup and output profile)
  //    and backup files (environment, selection and world)
  // -------------------------------------------------------------------------
  gzFile exp_s_file, out_p_file;
  gzFile env_file, exp_backup_file, world_file;
  open_setup_files(exp_s_file, out_p_file, t0, "r", dir);
  open_backup_files(env_file, exp_backup_file, world_file, t0, "r", dir);


  // -------------------------------------------------------------------------
  // 2) Load data from backup and parameter files
  // -------------------------------------------------------------------------
  load(env_file, exp_s_file, exp_backup_file,
       world_file, out_p_file, verbose, to_be_run);


  // -------------------------------------------------------------------------
  // 3) Close setup and backup files
  // -------------------------------------------------------------------------
  close_setup_files(exp_s_file, out_p_file);
  close_backup_files(env_file, exp_backup_file, world_file);


  // ---------------------------------------------------------------------------
  // 4) Recompute unsaved data
  // ---------------------------------------------------------------------------
  // Evaluate individuals
  world_->evaluate_individuals(get_env());
}


/*!
  \brief Load an experiment with the provided constitutive files
 */
void ae_exp_manager::load(int64_t t0,
                          char* env_file_name,
                          char* exp_setup_file_name,
                          char* exp_backup_file_name,
                          char* world_file_name,
                          char* out_prof_file_name,
                          bool verbose /*= false*/,
                          bool to_be_run /*= true*/)
{
  Time::set_time(t0);

  // ---------------------------------------------------------------------------
  // Open files and check them
  // ---------------------------------------------------------------------------
  gzFile exp_setup_file = gzopen(exp_setup_file_name, "r");
  gzFile out_prof_file = gzopen(out_prof_file_name, "r");
  gzFile env_file = gzopen(env_file_name, "r");
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
  if (env_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, env_file_name);
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
  load(env_file, exp_setup_file, exp_backup_file,
      world_file, out_prof_file, verbose, to_be_run);


  // ---------------------------------------------------------------------------
  // Close setup and backup files
  // ---------------------------------------------------------------------------
  gzclose(env_file);
  gzclose(exp_setup_file);
  gzclose(exp_backup_file);
  gzclose(world_file);
  gzclose(out_prof_file);


  // ---------------------------------------------------------------------------
  // 5) Recompute unsaved data
  // ---------------------------------------------------------------------------
  // Evaluate individuals
  world_->evaluate_individuals(get_env());
}

/*!
  Run the simulation
 */
void ae_exp_manager::run_evolution(void)
{
  // We are running a simulation.
  // Save the setup files to keep track of the setup history
  write_setup_files();

  // Dump the initial state of the population; useful for restarts
  _output_m->write_current_generation_outputs();

  while (Time::get_time() < t_end)
  {
    if (quit_signal_received()) break;

    printf("============================== %" PRId64 " ==============================\n",
        Time::get_time());
    printf("  Best individual's distance to target (metabolic) : %f\n",
        get_best_indiv()->get_dist_to_target_by_feature(METABOLISM));

    #ifdef __X11
      display();
    #endif

    // Take one step in the evolutionary loop
    step_to_next_generation();

    // Write statistical data and store phylogenetic data (tree)
    _output_m->write_current_generation_outputs();
  }

  _output_m->flush();


  printf("============================== %" PRId64 " ==============================\n",
      Time::get_time());
  printf("  Best individual's distance to target (metabolic) : %f\n",
      get_best_indiv()->get_dist_to_target_by_feature(METABOLISM));
  printf("===================================================================\n");
  printf("  The run is finished. \n");
  printf("  Printing the final best individual into " BEST_LAST_ORG_FNAME "\n");
  FILE* org_file = fopen(BEST_LAST_ORG_FNAME, "w");
  fputs(get_best_indiv()->get_genetic_unit_sequence(0), org_file);
  fclose(org_file);
}

void ae_exp_manager::update_best(void)
{
  world_->update_best();
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================
void ae_exp_manager::create_missing_directories(const char* dir /*= "."*/) const
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
  // Environment
  sprintf(cur_dir_name, "%s/" ENV_DIR, dir);
  status = mkdir(cur_dir_name, 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }
  // Population
  sprintf(cur_dir_name, "%s/" POP_DIR, dir);
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
}

void ae_exp_manager::open_backup_files(gzFile& env_file,
                                       gzFile& exp_backup_file,
                                       gzFile& world_file,
                                       int64_t t,
                                       const char mode[3],
                                       const char* dir /*= "."*/) const
{
  assert(strcmp(mode, "w") == 0 or strcmp(mode, "r") == 0);

  // -------------------------------------------------------------------------
  // 1) Generate backup file names for mandatory files.
  // -------------------------------------------------------------------------
  char env_file_name[255];
  char exp_backup_file_name[255];
  char world_file_name[255];

  sprintf(env_file_name, "%s/" ENV_FNAME_FORMAT, dir, t);
  sprintf(exp_backup_file_name, "%s/" EXP_S_FNAME_FORMAT, dir, t);
  sprintf(world_file_name, "%s/" WORLD_FNAME_FORMAT, dir, t);


  // -------------------------------------------------------------------------
  // 2) Open backup files (environment, population and selection)
  // -------------------------------------------------------------------------
  env_file = gzopen(env_file_name, mode);
  exp_backup_file = gzopen(exp_backup_file_name, mode);
  world_file = gzopen(world_file_name, mode);


  // -------------------------------------------------------------------------
  // 3) Check that files were correctly opened
  // -------------------------------------------------------------------------
  if (env_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, env_file_name);
    exit(EXIT_FAILURE);
  }
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

void ae_exp_manager::close_backup_files(gzFile& env_file,
                                        gzFile& exp_backup_file,
                                        gzFile& world_file) const
{
  gzclose(env_file);
  gzclose(exp_backup_file);
  gzclose(world_file);
}

void ae_exp_manager::open_setup_files(
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

  // 2) Open backup files (environment and world)
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

void ae_exp_manager::close_setup_files(gzFile& exp_s_file,
                                       gzFile& out_p_file) const
{
  gzclose(exp_s_file);
  gzclose(out_p_file);
}
// ===========================================================================
//                               Non inline accessors
// ===========================================================================
} // namespace aevol
