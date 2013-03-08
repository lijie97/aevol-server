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
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_exp_manager.h>
#include <ae_population.h>

#ifdef __X11
  #include <ae_population_X11.h>
#endif

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
  _num_gener    = 0;
  _first_gener  = 0;
  _last_gener   = 0;
  
  #ifdef __X11
    _pop = new ae_population_X11( this );
  #else
    _pop = new ae_population( this );
  #endif
  
  _exp_s    = new ae_exp_setup( this );
  _output_m = new ae_output_manager( this );
  
  _quit_signal_received = false;
}

// ===========================================================================
//                                  Destructors
// ===========================================================================
ae_exp_manager::~ae_exp_manager( void )
{
  delete _pop;
  delete _exp_s;
  delete _output_m;
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

  \see load( int32_t first_gener,
                           char* exp_setup_file_name,
                           char* out_prof_file_name,
                           char* env_file_name,
                           char* pop_file_name,
                           char* sp_struct_file_name,
                           bool verbose )
  \see save( void )
  \see save_copy( char* dir, int32_t num_gener )
*/
void ae_exp_manager::write_setup_files( void )
{
  // 1) Create missing directories
  create_missing_directories();
  
  // 2) Open setup files (experimental setup and output profile)
  gzFile exp_s_gzfile, out_p_gzfile;
  FILE*  exp_s_txtfile, * out_p_txtfile;
  open_setup_files( exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile, _num_gener );
  
  // 4) Write setup data
  _exp_s->write_setup_file( exp_s_gzfile );
  _output_m->write_setup_file( out_p_gzfile );
  _exp_s->write_setup_file( exp_s_txtfile );
  _output_m->write_setup_file( out_p_txtfile );
  
  // 5) Close setup files
  close_setup_files( exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile );
  
  //~ // 6) Link current setup files to those just written
  //~ unlink( EXP_S_CUR_FNAME );
  //~ unlink( OUT_P_CUR_FNAME );
  //~ unlink( EXP_S_TXT_CUR_FNAME );
  //~ unlink( OUT_P_TXT_CUR_FNAME );
  //~ symlink( exp_s_gzfile_name, EXP_S_CUR_FNAME );
  //~ symlink( out_p_gzfile_name, OUT_P_CUR_FNAME );
  //~ symlink( exp_s_txtfile_name, EXP_S_TXT_CUR_FNAME );
  //~ symlink( out_p_txtfile_name, OUT_P_TXT_CUR_FNAME );
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

  \see load( int32_t first_gener,
                           char* exp_setup_file_name,
                           char* out_prof_file_name,
                           char* env_file_name,
                           char* pop_file_name,
                           char* sp_struct_file_name,
                           bool verbose )
  \see write_setup_files( void )
  \see save_copy( char* dir, int32_t num_gener )
*/
void ae_exp_manager::save( void ) const
{
  // 1) Create missing directories
  create_missing_directories();
  
  // 2) Open backup files (environment, population and spatial structure)
  gzFile env_file, pop_file, sp_struct_file;
  open_backup_files( env_file, pop_file, sp_struct_file, _num_gener );
  
  // 3) Write the state of the environment, population and spatial structure into the backups
  get_env()->save( env_file );
  get_pop()->save( pop_file );
  if ( is_spatially_structured() )
  {
    get_spatial_structure()->save( sp_struct_file );
  }
  
  // 4) Close backup files
  close_backup_files( env_file, pop_file, sp_struct_file );
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
                           char* sp_struct_file_name,
                           bool verbose )
  \see write_setup_files( void )
  \see save( void )
*/
void ae_exp_manager::save_copy( char* dir, int32_t num_gener /*= 0*/ ) const
{
  // 1) Create missing directories
  create_missing_directories( dir );
  
  // 2) Open setup files (experimental setup and output profile)
  //    and backup files (environment, population and spacial structure)
  gzFile exp_s_gzfile, out_p_gzfile;
  FILE*  exp_s_txtfile, * out_p_txtfile;
  open_setup_files( exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile, num_gener, dir );
  gzFile env_file, pop_file, sp_struct_file;
  open_backup_files( env_file, pop_file, sp_struct_file, num_gener, dir );
  
  // 3) Write setup data
  _exp_s->write_setup_file( exp_s_gzfile );
  _output_m->write_setup_file( out_p_gzfile );
  
  // 4) Write the state of the environment, population and spatial structure into the backups
  get_env()->save( env_file );
  get_pop()->save( pop_file );
  if ( is_spatially_structured() )
  {
    get_spatial_structure()->save( sp_struct_file );
  }
  
  // 4) Close setup and backup files
  close_setup_files( exp_s_gzfile, exp_s_txtfile, out_p_gzfile, out_p_txtfile );
  close_backup_files( env_file, pop_file, sp_struct_file );
}


/*!
  \brief Load an experiment with the provided constitutive files
 */
void ae_exp_manager::load( int32_t first_gener,
                           char* exp_setup_file_name,
                           char* out_prof_file_name,
                           char* env_file_name,
                           char* pop_file_name,
                           char* sp_struct_file_name,
                           bool verbose )
{
  _first_gener = first_gener;
  
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
  
  if ( gzdirect( exp_setup_gzfile ) )
  {
    // The provided output profile file is not zipped, it is a plain text setup file
    // => Close gzfile and open as plain file instead
    gzclose( exp_setup_gzfile );
    exp_setup_gzfile = NULL;
    exp_setup_txtfile = fopen( exp_setup_file_name, "r" );
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
  // 2) Load data from backup and parameter files
  // ---------------------------------------------------------------------------
  // Open backup files
  gzFile env_file = gzopen( env_file_name, "r" );
  gzFile pop_file = gzopen( pop_file_name, "r" );
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
  if ( sp_struct_file_name != NULL )
  {
    sp_struct_file = gzopen( sp_struct_file_name, "r" );
    if ( sp_struct_file == Z_NULL )
    {
      printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, sp_struct_file_name );
      exit( EXIT_FAILURE );
    }
  }
  
  // Load experimental setup, population data and output profile from the backups
  if ( exp_setup_gzfile != NULL )
  {
    _exp_s->load( exp_setup_gzfile, env_file, sp_struct_file, verbose );
  }
  else
  {
    _exp_s->load( exp_setup_txtfile, env_file, sp_struct_file, verbose );
  }
  
  _pop->load( pop_file, verbose );
  
  if ( out_prof_gzfile != NULL )
  {
    _output_m->load( out_prof_gzfile, verbose );
  }
  else
  {
    _output_m->load( out_prof_txtfile, verbose );
  }
  
  // Close backup files
  gzclose( env_file );
  gzclose( pop_file );
  gzclose( sp_struct_file );
  gzclose( exp_setup_gzfile );
  gzclose( out_prof_gzfile );
  

  // ---------------------------------------------------------------------------
  // 3) Recompute unsaved data
  // ---------------------------------------------------------------------------
  // Evaluate individuals
  _pop->evaluate_individuals( _exp_s->get_env() );

  // If the population is spatially structured, then the individuals are saved
  // and loaded in the order of the grid and not in increasing order of fitness
  // so we have to sort the individuals
  if ( is_spatially_structured() )
  {
    _pop->sort_individuals();
  }


  // If the simulation is being continued (not just post-processed),
  // prepare output data accordingly
  //~ if ( to_be_run )
  //~ {
    //~ // Prepare stat files
    //~ _stats  = new ae_stats( _first_gener );
    
    //~ // Prepare tree
    //~ if ( ae_common::rec_params->get_record_tree() == true )
    //~ { 
      //~ mkdir( "tree", 0755 );
      //~ _tree  = new ae_tree(); 
    //~ }
    //~ else
    //~ {
      //~ _tree = NULL;
    //~ }
    
    //~ // Prepare dump
    //~ if ( ae_common::rec_params->get_dump_period() > 0 )
    //~ {
      //~ mkdir( "dump", 0755 );
      //~ _dump = new ae_dump();
    //~ }
    //~ else
    //~ {
      //~ _dump = NULL;
    //~ }
    
    //~ mkdir( "backup", 0755 );
    
    //~ if ( ae_common::rec_params->is_logged( LOG_LOADS ) == true )
    //~ {
      //~ // Write an entry in the LOADS log file
      //~ fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "GENERATION_OVERLOAD %"PRId32"\n", _num_gener );
      //~ if ( param_overloader->get_nb_overloaded() > 0 )
      //~ {
        //~ //fprintf( _logs->get_log( LOG_LOADS ), "  Overloaded parameters:\n" );
        //~ param_overloader->write_log( ae_common::rec_params->get_log( LOG_LOADS ) );
        //~ fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "\n" );
      //~ }
      //~ else
      //~ {
        //~ fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "  No overloaded parameters\n\n" );
      //~ }
    //~ }
  //~ }
  //~ else
  //~ {
    //~ // We just want to inspect the state of the simulation at this moment
    //~ ae_common::rec_params->init_logs( 0 );
  //~ }
  
  // Initialize display
  #ifdef __X11
    //~ ((ae_exp_setup_X11*) this)->initialize( ae_common::pop_structure, ae_common::params->get_allow_plasmids() );
  #endif // def __X11
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
  
  while ( _num_gener < _last_gener )
  {
    printf( "============================== %"PRId32" ==============================\n", _num_gener );
    printf( "  Best individual's distance to target (metabolic) : %f\n",
            _pop->get_best()->get_dist_to_target_by_feature( METABOLISM ) );

    if ( quit_signal_received() ) break;
    
    #ifdef __X11
      display();
    #endif
    
    // Take one step in the evolutionary loop
    _exp_s->step_to_next_generation();
    _num_gener++;

    // Write statistical data and store phylogenetic data (tree)
    _output_m->write_current_generation_outputs();
  }
  
  _output_m->flush();
  
  printf( "============================== %"PRId32" ==============================\n", _num_gener );
  printf( "  distance to target (metabolic) : %f\n", _pop->get_best()->get_dist_to_target_by_feature( METABOLISM ) );
  printf( "===================================================================\n");
  printf ("  The run is finished. \n"); 
  printf ("  Printing the final best individual into "BEST_LAST_ORG_FNAME"\n"); 
  FILE* org_file = fopen( BEST_LAST_ORG_FNAME, "w" );
  fputs( _pop->get_best()->get_genetic_unit(0)->get_dna()->get_data(), org_file ); 
  fclose ( org_file ); 
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================
void ae_exp_manager::create_missing_directories( const char* dir /*= "."*/ ) const
{
  char cur_dir_name[255];
  int status;
  
  // Experimental setup
  sprintf( cur_dir_name, "%s/"EXP_S_DIR, dir );
  status = mkdir( cur_dir_name, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, cur_dir_name, errno );
  }
  // Output profile
  sprintf( cur_dir_name, "%s/"OUT_P_DIR, dir );
  status = mkdir( cur_dir_name, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, cur_dir_name, errno );
  }
  // Environment
  sprintf( cur_dir_name, "%s/"ENV_DIR, dir );
  status = mkdir( cur_dir_name, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, cur_dir_name, errno );
  }
  // Population
  sprintf( cur_dir_name, "%s/"POP_DIR, dir );
  status = mkdir( cur_dir_name, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, cur_dir_name, errno );
  }
  // Spatial structure
  if ( is_spatially_structured() )
  {
    sprintf( cur_dir_name, "%s/"SP_STRUCT_DIR, dir );
    status = mkdir( cur_dir_name, 0755 );
    if ( status == -1 && errno != EEXIST )
    {
      err( EXIT_FAILURE, cur_dir_name, errno );
    }
  }
}

void ae_exp_manager::open_backup_files( gzFile& env_file,
                                        gzFile& pop_file,
                                        gzFile& sp_struct_file,
                                        int32_t num_gener,
                                        const char* dir /*= "."*/ ) const
{ 
  // 1) Generate backup file names
  char env_file_name[255];
  char pop_file_name[255];
  char sp_struct_file_name[255];
  
  sprintf( env_file_name, "%s/"ENV_FNAME_FORMAT, dir, num_gener );
  sprintf( pop_file_name, "%s/"POP_FNAME_FORMAT, dir, num_gener );
	
  if ( is_spatially_structured() )
  {
    sprintf( sp_struct_file_name, "%s/"SP_STRUCT_FNAME_FORMAT, dir, num_gener );
  }
  
  // 2) Open backup files (environment, population and spacial structure)
  env_file = gzopen( env_file_name, "w" );
  pop_file = gzopen( pop_file_name, "w" );
  if ( is_spatially_structured() )
  {
    sp_struct_file  = gzopen( sp_struct_file_name, "w" );
  }
  
  // 3) Check that files were correctly opened
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
  if ( is_spatially_structured() and sp_struct_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, sp_struct_file_name );
    exit( EXIT_FAILURE );
  }
}

void ae_exp_manager::close_backup_files(  gzFile& env_file,
                                          gzFile& pop_file,
                                          gzFile& sp_struct_file ) const
{
  gzclose( env_file );
  gzclose( pop_file );
  if ( is_spatially_structured() )
  {
    gzclose( sp_struct_file );
  }
}

void ae_exp_manager::open_setup_files(
        gzFile& exp_s_gzfile, FILE*& exp_s_txtfile,
        gzFile& out_p_gzfile, FILE*& out_p_txtfile,
        int32_t num_gener,
        const char* dir /*= "."*/ ) const
{
  // 1) Generate setup file names
  char exp_s_gzfile_name[255];
  char exp_s_txtfile_name[255];
  char out_p_gzfile_name[255];
  char out_p_txtfile_name[255];
  sprintf( exp_s_gzfile_name,  "%s/"EXP_S_FNAME_FORMAT,     dir, num_gener );
  sprintf( exp_s_txtfile_name, "%s/"EXP_S_TXT_FNAME_FORMAT, dir, num_gener );
  sprintf( out_p_gzfile_name,  "%s/"OUT_P_FNAME_FORMAT,     dir, num_gener );
  sprintf( out_p_txtfile_name, "%s/"OUT_P_TXT_FNAME_FORMAT, dir, num_gener );
  
  // 2) Open backup files (environment, population and spacial structure)
  exp_s_gzfile  = gzopen( exp_s_gzfile_name, "w" );
  exp_s_txtfile = fopen( exp_s_txtfile_name, "w" );
  out_p_gzfile  = gzopen( out_p_gzfile_name, "w" );
  out_p_txtfile = fopen( out_p_txtfile_name, "w" );
  
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
