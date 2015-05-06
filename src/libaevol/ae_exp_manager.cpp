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
#include "ae_population.h"
#include "ae_individual.h"

#ifdef __X11
  #include "ae_population_X11.h"
#endif

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
ae_exp_manager::ae_exp_manager( void )
{
  // ------------------------------------------------------ Experimental setup
  _exp_s = new ae_exp_setup( this );

  // -------------------------------------------------------------- Population
  #ifdef __X11
    _pop = new ae_population_X11( this );
  #else
    _pop = new ae_population( this );
  #endif

  // ------------------------------------------------------------- Environment
  _env = new Environment();

  // ------------------------------------------------------- Spatial structure
  _spatial_structure = NULL;

  // ---------------------------------------------------------- Output manager
  _output_m = new ae_output_manager( this );


  // -------------------------------- Timestep up to which we want to simulate
  t_end = 0;


  // ------------------------------------------------------------- Quit signal
  _quit_signal_received = false;
}

// ===========================================================================
//                                  Destructor
// ===========================================================================
ae_exp_manager::~ae_exp_manager( void )
{
  delete _pop;
  delete _exp_s;
  delete _output_m;
  delete _env;
  if ( _spatial_structure != NULL ){delete _spatial_structure;}
}

// ===========================================================================
//                                   Algorithms
// ===========================================================================
void ae_exp_manager::foreach_indiv(void (*processor)(ae_individual& indiv)) const
{
  for (const auto& indiv: get_indivs_std())
    processor(*indiv);
}

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
            char* sp_struct_file_name,
            bool verbose)
  \see save(void)
  \see save_copy(char* dir, int64_t t)
*/
void ae_exp_manager::write_setup_files(void)
{
  // 1) Create missing directories
  create_missing_directories();

  // 2) Open setup files (experimental setup and output profile)
  gzFile exp_s_gzfile, out_p_gzfile;
  FILE*  exp_s_txtfile, * out_p_txtfile;
  open_setup_files(exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile,
      Time::get_time(), "w");

  // 4) Write setup data
  _exp_s->write_setup_file(exp_s_gzfile);
  _output_m->write_setup_file(out_p_gzfile);
  _exp_s->write_setup_file(exp_s_txtfile);
  _output_m->write_setup_file(out_p_txtfile);

  // 5) Close setup files
  close_setup_files(exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile);
}

/*!
  \brief Save all the dynamic data (the data that changes over the simulation)

  Save the state of
    * The population
    * The environment iff it is changing (i.e. it has variation and/or noise)
    * The spatial structure (if any)
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
            char* sp_struct_file_name,
            bool verbose)
  \see write_setup_files(void)
  \see save_copy(char* dir, int64_t t)
*/
void ae_exp_manager::save(void) const
{
  // 1) Create missing directories
  create_missing_directories();

  // 2) Open backup files (environment, population, selection and spatial structure)
  gzFile env_file, pop_file, sel_file, sp_struct_file;
  open_backup_files(env_file, pop_file, sel_file, sp_struct_file,
      Time::get_time(), "w");

  // 3) Write the state of the environment, population and spatial structure into the backups
  get_env()->save(env_file);
  get_pop()->save(pop_file);
  get_sel()->save(sel_file);
  if ( is_spatially_structured() )
  {
    get_spatial_structure()->save(sp_struct_file);
  }

  // 4) Close backup files
  close_backup_files( env_file, pop_file, sel_file, sp_struct_file );
}

/*!
  \brief Saves a complete copy of the experiment to the provided location.

  Save both the setup (constant) and the backup (dynamic) files to the <dir>
  directory. Dynamic file names will be appended with the <num_gener>

  Static files (saved both as gzipped binary and plain text)
    * The experimental setup
    * The output profile

  Dynamic files (saved as gzipped binary)
    * The population
    * The environment iff it is changing (i.e. it has variation and/or noise)
    * The spatial structure (if any)
    * The PRNG (random generator) used for the selection

  \see load( int32_t first_gener,
             char* exp_setup_file_name,
             char* out_prof_file_name,
             char* env_file_name,
             char* pop_file_name,
             char* sel_file_name,
             char* sp_struct_file_name,
             bool verbose )
  \see write_setup_files( void )
  \see save( void )
*/
void ae_exp_manager::save_copy( char* dir, int32_t num_gener /*= 0*/ ) const
{
  // 1) Create missing directories
  create_missing_directories(dir);

  // 2) Open setup files (experimental setup and output profile)
  //    and backup files (environment, population, selection and spacial structure)
  gzFile exp_s_gzfile, out_p_gzfile;
  FILE*  exp_s_txtfile, * out_p_txtfile;
  open_setup_files( exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile, num_gener, "w", dir );
  gzFile env_file, pop_file, sel_file, sp_struct_file;
  open_backup_files( env_file, pop_file, sel_file, sp_struct_file, num_gener, "w", dir );

  // 3) Write setup data
  _exp_s->write_setup_file( exp_s_gzfile );
  _output_m->write_setup_file( out_p_gzfile );

  // 4) Write the state of the environment, population, selection and spatial structure into the backups
  get_env()->save( env_file );
  get_pop()->save( pop_file );
  get_sel()->save( sel_file );
  if ( is_spatially_structured() )
  {
    get_spatial_structure()->save( sp_struct_file );
  }

  // 4) Close setup and backup files
  close_setup_files( exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile );
  close_backup_files( env_file, pop_file, sel_file, sp_struct_file );
}


/*!
  \brief Load an experiment with the provided files
 */
void ae_exp_manager::load( gzFile& pop_file,
                           gzFile& env_file,
                           gzFile& exp_s_gzfile,
                           FILE*&  exp_s_txtfile,
                           gzFile& exp_backup_file,
                           gzFile& sp_struct_file,
                           gzFile& out_p_gzfile,
                           FILE*& out_p_txtfile,
                           bool verbose,
                           bool to_be_run /*  = true */ )
{
  // ---------------------------------------- Retrieve experimental setup data
  printf( "  Loading experimental setup..." );
  fflush( stdout );
  if ( exp_s_gzfile != NULL )
  {
    _exp_s->load( exp_s_gzfile, exp_backup_file, verbose );
  }
  else
  {
    _exp_s->load( exp_s_txtfile, exp_backup_file, verbose );
  }
  printf( " OK\n" );

  // ----------------------------------------- Retrieve spatial structure data
  if ( sp_struct_file != NULL )
  {
    printf( "  Loading spatial structure..." );
    fflush( stdout );
    _spatial_structure = new ae_spatial_structure( sp_struct_file );
    printf( " OK\n" );
  }

  // ------------------------------------------------ Retrieve population data
  printf( "  Loading population..." );
  fflush( stdout );
  _pop->load( pop_file, verbose );
  printf( " OK\n" );

  // --------------------------------------------- Retrieve environmental data
  printf( "  Loading environment..." );
  fflush( stdout );
  //~ delete _env;
  //~ _env = new ae_environment();
  _env->load( env_file );
  printf( " OK\n" );

  // --------------------------------------------- Retrieve output profile data
  printf( "  Loading output profile..." );
  fflush( stdout );
  if ( out_p_gzfile != NULL )
  {
    _output_m->load( out_p_gzfile, verbose, to_be_run );
  }
  else
  {
    _output_m->load( out_p_txtfile, verbose, to_be_run );
  }
  printf( " OK\n" );
}


/*!
  \brief Load an experiment with default files from a given directory
 */
void ae_exp_manager::load( const char* dir,
    int64_t t0, bool use_text_files,
    bool verbose, bool to_be_run /*  = true */ )
{
  Time::set_time(t0);

  // -------------------------------------------------------------------------
  // 1) Open setup files (experimental setup and output profile)
  //    and backup files (environment, population, selection and sp structure)
  // -------------------------------------------------------------------------
  gzFile exp_s_gzfile, out_p_gzfile;
  FILE*  exp_s_txtfile, * out_p_txtfile;
  open_setup_files(exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile, t0, "r", dir );
  gzFile env_file, pop_file, exp_backup_file, sp_struct_file;
  open_backup_files(env_file, pop_file, exp_backup_file, sp_struct_file, t0, "r", dir );


  // -------------------------------------------------------------------------
  // 2) Load data from backup and parameter files
  // -------------------------------------------------------------------------
  load( pop_file, env_file, exp_s_gzfile, exp_s_txtfile, exp_backup_file,
        sp_struct_file, out_p_gzfile, out_p_txtfile, verbose, to_be_run );


  // -------------------------------------------------------------------------
  // 3) Close setup and backup files
  // -------------------------------------------------------------------------
  close_setup_files( exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile );
  close_backup_files( env_file, pop_file, exp_backup_file, sp_struct_file );


  // ---------------------------------------------------------------------------
  // 4) Recompute unsaved data
  // ---------------------------------------------------------------------------
  // Evaluate individuals
  _pop->evaluate_individuals( get_env() );

  // If the population is spatially structured, then the individuals are saved
  // and loaded in the order of the grid and not in increasing order of fitness
  // so we have to sort the individuals
  if ( is_spatially_structured() )
  {
    _pop->sort_individuals();
  }
}


/*!
  \brief Load an experiment with the provided constitutive files
 */
void ae_exp_manager::load(int64_t t0,
                          char* pop_file_name,
                          char* env_file_name,
                          char* exp_setup_file_name,
                          char* exp_backup_file_name,
                          char* sp_struct_file_name,
                          char* out_prof_file_name,
                          bool verbose ,
                          bool to_be_run /* = true */)
{
  Time::set_time(t0);

  // ---------------------------------------------------------------------------
  // 1) Determine whether the parameter files are in binary or plain text format
  //    and open in the corresponding "mode"
  // ---------------------------------------------------------------------------
  gzFile exp_setup_gzfile = NULL;
  gzFile out_prof_gzfile = NULL;
  FILE* exp_setup_txtfile = NULL;
  FILE* out_prof_txtfile = NULL;

  exp_setup_gzfile = gzopen( exp_setup_file_name, "r" );
  out_prof_gzfile = gzopen( out_prof_file_name, "r" );

  if ( exp_setup_gzfile == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, exp_setup_file_name );
    exit( EXIT_FAILURE );
  }
  if ( out_prof_gzfile == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, out_prof_file_name );
    exit( EXIT_FAILURE );
  }

  if (gzdirect(exp_setup_gzfile))
  {
    // The provided output profile file is not zipped, it is a plain text setup file
    // => Close gzfile and open as plain file instead
    gzclose(exp_setup_gzfile);
    exp_setup_gzfile = NULL;
    exp_setup_txtfile = fopen(exp_setup_file_name, "r");
  }

  if ( gzdirect( out_prof_gzfile ) )
  {
    // The provided output profile file is not zipped, it is a plain text setup file
    // => Close gzfile and open as plain file instead
    gzclose( out_prof_gzfile );
    out_prof_gzfile = NULL;
    out_prof_txtfile = fopen( out_prof_file_name, "r" );
  }


  // ---------------------------------------------------------------------------
  // 2) Open backup and parameter files
  // ---------------------------------------------------------------------------
  gzFile env_file = gzopen( env_file_name, "r" );
  gzFile pop_file = gzopen( pop_file_name, "r" );
  gzFile exp_backup_file = gzopen( exp_backup_file_name, "r" );
  gzFile sp_struct_file  = NULL;

  if ( env_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, env_file_name );
    exit( EXIT_FAILURE );
  }
  if ( pop_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, pop_file_name );
    exit( EXIT_FAILURE );
  }
  if ( exp_backup_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, exp_backup_file_name );
    exit( EXIT_FAILURE );
  }
  if ( sp_struct_file_name != NULL )
  {
    sp_struct_file = gzopen( sp_struct_file_name, "r" );
    if ( sp_struct_file == Z_NULL )
    {
      printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, sp_struct_file_name );
      exit( EXIT_FAILURE );
    }
  }


  // ---------------------------------------------------------------------------
  // 3) Load data from backup and parameter files
  // ---------------------------------------------------------------------------

  // ----------------------------------------- Retrieve spatial structure data
  if ( sp_struct_file != NULL )
  {
    printf( "  Loading spatial structure..." );
    fflush( stdout );
    _spatial_structure = new ae_spatial_structure( sp_struct_file );
    printf( "OK\n" );
  }

  // ------------------------------------------------ Retrieve population data
  printf( "  Loading population..." );
  fflush( stdout );
  _pop->load( pop_file, verbose );
  printf( "OK\n" );

  // --------------------------------------------- Retrieve environmental data
  printf( "  Loading environment..." );
  fflush( stdout );
  //~ delete _env;
  //~ _env = new ae_environment();
  _env->load( env_file );
  printf( "OK\n" );

  // ---------------------------------------- Retrieve experimental setup data
  printf( "  Loading experimental setup..." );
  fflush( stdout );
  if ( exp_setup_gzfile != NULL )
  {
    _exp_s->load( exp_setup_gzfile, exp_backup_file, verbose );
  }
  else
  {
    _exp_s->load( exp_setup_txtfile, exp_backup_file, verbose );
  }
  printf( "OK\n" );

  // --------------------------------------------- Retrieve output profile data
  if ( out_prof_gzfile != NULL )
  {
    _output_m->load( out_prof_gzfile, verbose, to_be_run );
  }
  else
  {
    _output_m->load( out_prof_txtfile, verbose, to_be_run );
  }


  // ---------------------------------------------------------------------------
  // 4) Close setup and backup files
  // ---------------------------------------------------------------------------
  gzclose(pop_file);
  gzclose(env_file);
  gzclose(exp_setup_gzfile);
  gzclose(exp_backup_file);
  gzclose(sp_struct_file);
  gzclose(out_prof_gzfile);


  // ---------------------------------------------------------------------------
  // 5) Recompute unsaved data
  // ---------------------------------------------------------------------------
  // Evaluate individuals
  _pop->evaluate_individuals( get_env() );

  // If the population is spatially structured, then the individuals are saved
  // and loaded in the order of the grid and not in increasing order of fitness
  // so we have to sort the individuals
  if ( is_spatially_structured() )
  {
    _pop->sort_individuals();
  }
}

/*!
  Run the simulation
 */
void ae_exp_manager::run_evolution( void )
{
  // We are running a simulation.
  // Save the setup files to keep track of the setup history
  write_setup_files();

  // Dump the initial state of the population; useful for restarts
  _output_m->write_current_generation_outputs();

  while (true) { // termination condition is into the loop
    printf("============================== %" PRId64 " ==============================\n",
           Time::get_time());
    printf("  Best individual's distance to target (metabolic) : %f\n",
           _pop->get_best()->get_dist_to_target_by_feature(METABOLISM));

    if (Time::get_time() >= t_end or quit_signal_received())
      break;

#ifdef __X11
    display();
#endif

    // Take one step in the evolutionary loop
    step_to_next_generation();

    // Write statistical data and store phylogenetic data (tree)
    _output_m->write_current_generation_outputs();
  }

  _output_m->flush();
  printf("================================================================\n");
  printf("  The run is finished. \n");
  printf("  Printing the final best individual into " BEST_LAST_ORG_FNAME "\n");
  FILE* org_file = fopen(BEST_LAST_ORG_FNAME, "w");
  fputs(_pop->get_best()->get_genetic_unit(0).get_dna()->get_data(), org_file);
  fclose(org_file);
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================
void ae_exp_manager::create_missing_directories( const char* dir /*= "."*/ ) const
{
  char cur_dir_name[255];
  int status;

  // Base directory
  status = mkdir( dir, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, cur_dir_name, errno );
  }

  // Experimental setup
  sprintf( cur_dir_name, "%s/" EXP_S_DIR, dir );
  status = mkdir( cur_dir_name, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, cur_dir_name, errno );
  }
  // Output profile
  sprintf( cur_dir_name, "%s/" OUT_P_DIR, dir );
  status = mkdir( cur_dir_name, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, cur_dir_name, errno );
  }
  // Environment
  sprintf( cur_dir_name, "%s/" ENV_DIR, dir );
  status = mkdir( cur_dir_name, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, cur_dir_name, errno );
  }
  // Population
  sprintf( cur_dir_name, "%s/" POP_DIR, dir );
  status = mkdir( cur_dir_name, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, cur_dir_name, errno );
  }
  // Spatial structure
  if ( is_spatially_structured() )
  {
    sprintf( cur_dir_name, "%s/" SP_STRUCT_DIR, dir );
    status = mkdir( cur_dir_name, 0755 );
    if ( status == -1 && errno != EEXIST )
    {
      err( EXIT_FAILURE, cur_dir_name, errno );
    }
  }
}

void ae_exp_manager::open_backup_files(gzFile& env_file,
                                       gzFile& pop_file,
                                       gzFile& exp_backup_file,
                                       gzFile& sp_struct_file,
                                       int64_t t,
                                       const char mode[3],
                                       const char* dir /*= "."*/ ) const
{
  assert( strcmp( mode, "w" ) == 0 or strcmp( mode, "r" ) == 0 );

  // -------------------------------------------------------------------------
  // 1) Generate backup file names for mandatory files.
  //    Optional files (e.g. spatial structure) will be handled separately
  // -------------------------------------------------------------------------
  char env_file_name[255];
  char pop_file_name[255];
  char exp_backup_file_name[255];

  sprintf(env_file_name, "%s/" ENV_FNAME_FORMAT, dir, t);
  sprintf(pop_file_name, "%s/" POP_FNAME_FORMAT, dir, t);
  sprintf(exp_backup_file_name, "%s/" EXP_S_FNAME_FORMAT, dir, t);


  // -------------------------------------------------------------------------
  // 2) Open backup files (environment, population and selection)
  // -------------------------------------------------------------------------
  env_file = gzopen(env_file_name, mode);
  pop_file = gzopen(pop_file_name, mode);
  exp_backup_file = gzopen(exp_backup_file_name, mode);


  // -------------------------------------------------------------------------
  // 3) Check that files were correctly opened
  // -------------------------------------------------------------------------
  if ( env_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, env_file_name );
    exit( EXIT_FAILURE );
  }
  if ( pop_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, pop_file_name );
    exit( EXIT_FAILURE );
  }
  if ( exp_backup_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, exp_backup_file_name );
    exit( EXIT_FAILURE );
  }


  // -------------------------------------------------------------------------
  // 4) Handle optional files
  // -------------------------------------------------------------------------
  char sp_struct_file_name[255];
  sprintf(sp_struct_file_name, "%s/" SP_STRUCT_FNAME_FORMAT, dir, t);

  // Check whether there should be spatial structure
  bool spatially_structured = false;
  if ( strcmp( mode, "w" ) == 0 )
  {
    if ( is_spatially_structured() )
      spatially_structured = true;
  }
  else if ( strcmp( mode, "r" ) == 0 )
  {
    // Check for the existence of the corresponding file in the file system.
    struct stat stat_buf;
    if ( stat( sp_struct_file_name, &stat_buf ) == -1 )
    {
      if ( errno != ENOENT )
      {
        printf( "%s:%d: error: unknown error.\n", __FILE__, __LINE__ );
        exit( EXIT_FAILURE );
      }
      // Else file does not exist, leave spatially_structured to false
    }
    else
    {
      spatially_structured = true;
    }
  }

  // If there should be spatial structure, open file and check
  if ( spatially_structured )
  {
    sp_struct_file = gzopen( sp_struct_file_name, mode );

    if ( sp_struct_file == Z_NULL )
    {
      printf( "%s:%d: error: could not open backup file %s\n",
              __FILE__, __LINE__, sp_struct_file_name );
      exit( EXIT_FAILURE );
    }
  }
  else
  {
    sp_struct_file = NULL;
  }
}

void ae_exp_manager::close_backup_files(  gzFile& env_file,
                                          gzFile& pop_file,
                                          gzFile& exp_backup_file,
                                          gzFile& sp_struct_file ) const
{
  gzclose( env_file );
  gzclose( pop_file );
  gzclose( exp_backup_file );
  if ( is_spatially_structured() )
  {
    gzclose( sp_struct_file );
  }
}

void ae_exp_manager::open_setup_files(
        gzFile& exp_s_gzfile, FILE*& exp_s_txtfile,
        gzFile& out_p_gzfile, FILE*& out_p_txtfile,
        int64_t t,
        const char mode[3],
        const char* dir /*= "."*/ ) const
{
  // 1) Generate setup file names
  char exp_s_gzfile_name[255];
  char exp_s_txtfile_name[255];
  char out_p_gzfile_name[255];
  char out_p_txtfile_name[255];

  sprintf( exp_s_gzfile_name,  "%s/" EXP_S_CONST_FNAME_FORMAT,     dir );
  sprintf( exp_s_txtfile_name, "%s/" EXP_S_TXT_CONST_FNAME_FORMAT, dir );
  sprintf( out_p_gzfile_name,  "%s/" OUT_P_FNAME_FORMAT,     dir );
  sprintf( out_p_txtfile_name, "%s/" OUT_P_TXT_FNAME_FORMAT, dir );

  // 2) Open backup files (environment, population and spacial structure)
  exp_s_gzfile  = gzopen( exp_s_gzfile_name, mode );
  exp_s_txtfile = fopen( exp_s_txtfile_name, mode );
  out_p_gzfile  = gzopen( out_p_gzfile_name, mode );
  out_p_txtfile = fopen( out_p_txtfile_name, mode );

  // 3) Check that files were correctly opened
  if ( exp_s_gzfile == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, exp_s_gzfile_name );
    exit( EXIT_FAILURE );
  }
  if ( exp_s_txtfile == NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, exp_s_txtfile_name );
    exit( EXIT_FAILURE );
  }
  if ( out_p_gzfile == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, out_p_gzfile_name );
    exit( EXIT_FAILURE );
  }
  if ( out_p_txtfile == NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, out_p_txtfile_name );
    exit( EXIT_FAILURE );
  }
}

void ae_exp_manager::close_setup_files(
        gzFile& exp_s_gzfile, FILE* exp_s_txtfile,
        gzFile& out_p_gzfile, FILE* out_p_txtfile ) const
{
  gzclose( exp_s_gzfile );
  fclose( exp_s_txtfile );
  gzclose( out_p_gzfile );
  fclose( out_p_txtfile );
}
// ===========================================================================
//                               Non inline accessors
// ===========================================================================
} // namespace aevol
