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


/*! \class ae_exp_manager
    \brief This class allows to manage an experiment.
    
    An experiment manager allows to... manage an experiment.
    It owns a population and an experimental_setup that can be loaded from a pair of aevol binary files (pop and exp_setup)
*/
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <sys/stat.h>
#include <err.h>
#include <errno.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_exp_manager.h>

#ifdef __X11
  #include <ae_population_X11.h>
#endif

//##############################################################################
//                                                                             #
//                             Class ae_exp_manager                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
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

// =================================================================
//                             Destructors
// =================================================================
ae_exp_manager::~ae_exp_manager( void )
{
  delete _pop;
  delete _exp_s;
  delete _output_m;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_exp_manager::write_setup_files( void )
{
  // 1) Generate file names
  char* exp_s_gzfile_name = new char[63];
  char* out_m_gzfile_name  = new char[63];
  char* exp_s_txtfile_name = new char[63];
  char* out_m_txtfile_name  = new char[63];
  strcpy( exp_s_gzfile_name,  "exp_setup.ae" );
  strcpy( out_m_gzfile_name,  "output_profile.ae" );
  strcpy( exp_s_txtfile_name,  "exp_setup.in" );
  strcpy( out_m_txtfile_name,  "output_profile.in" );
  
  // 2) Open files
  gzFile* exp_s_gzfile  = (gzFile*) gzopen( exp_s_gzfile_name, "w" );
  gzFile* out_m_gzfile  = (gzFile*) gzopen( out_m_gzfile_name, "w" );
  FILE*   exp_s_txtfile = fopen( exp_s_txtfile_name, "w" );
  FILE*   out_m_txtfile = fopen( out_m_txtfile_name, "w" );
  
  // 4) Write data
  _exp_s->write_setup_file( exp_s_gzfile );
  _output_m->write_setup_file( out_m_gzfile );
  
  _exp_s->write_setup_file( exp_s_txtfile );
  _output_m->write_setup_file( out_m_txtfile );
  
  // 4) Close files
  gzclose( exp_s_gzfile );
  gzclose( out_m_gzfile );
  fclose( exp_s_txtfile );
  fclose( out_m_txtfile );
  
  // 5) Clean up
  delete [] exp_s_gzfile_name;
  delete [] out_m_gzfile_name;
  delete [] exp_s_txtfile_name;
  delete [] out_m_txtfile_name;
}

void ae_exp_manager::save_experiment( void ) const
{
  // 1) Generate backup file names
  char env_file_name[50];
  char pop_file_name[50];
  char sp_struct_file_name[50];
  
  sprintf( env_file_name,       ENV_FNAME_FORMAT,       _num_gener );
  sprintf( pop_file_name,       POP_FNAME_FORMAT,       _num_gener );
  if ( is_spatially_structured() )
  {
    sprintf( sp_struct_file_name, SP_STRUCT_FNAME_FORMAT, _num_gener );
  }
  
  
  // 2) Create missing backup directories
  int status;
  status = mkdir( ENV_DIR, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, ENV_DIR, errno );
  }
  status = mkdir( POP_DIR, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, POP_DIR, errno );
  }
  if ( is_spatially_structured() )
  {
    status = mkdir( SP_STRUCT_DIR, 0755 );
    if ( status == -1 && errno != EEXIST )
    {
      err( EXIT_FAILURE, SP_STRUCT_DIR, errno );
    }
  }
  
  
  // 3) Open backup files (population, exp_setup and output_man)
  gzFile* env_file        = (gzFile*) gzopen( env_file_name, "w" );
  gzFile* pop_file        = (gzFile*) gzopen( pop_file_name, "w" );
  gzFile* sp_struct_file;
  if ( is_spatially_structured() )
  {
    sp_struct_file  = (gzFile*) gzopen( sp_struct_file_name, "w" );
  }
  
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
  if ( is_spatially_structured() )
  {
    if ( sp_struct_file == Z_NULL )
    {
      printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, sp_struct_file_name );
      exit( EXIT_FAILURE );
    }
  }
  
  // 4) Write the state of the environment, population and spatial structure into the backups
  get_env()->save( env_file );
  get_pop()->save( pop_file );
  if ( is_spatially_structured() )
  {
    get_spatial_structure()->save( sp_struct_file );
  }
  
  // 5) Close backup files
  gzclose( env_file );
  gzclose( pop_file );
  if ( is_spatially_structured() )
  {
    gzclose( sp_struct_file );
  }
}

void ae_exp_manager::load_experiment( char* exp_setup_file_name,
                                      char* out_prof_file_name,
                                      char* env_file_name,
                                      char* pop_file_name,
                                      char* sp_struct_file_name,
                                      bool verbose )
{
  // ---------------------------------------------------------------------------
  // 1) Determine whether the parameter files are plain text and open in the 
  //    corresponding "mode"
  // ---------------------------------------------------------------------------
  gzFile* exp_setup_gzfile = NULL;
  gzFile* out_prof_gzfile = NULL;
  FILE* exp_setup_txtfile = NULL;
  FILE* out_prof_txtfile = NULL;
  
  exp_setup_gzfile = (gzFile*) gzopen( exp_setup_file_name, "r" );
  out_prof_gzfile = (gzFile*) gzopen( out_prof_file_name, "r" );
  
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
  gzFile* env_file = (gzFile*) gzopen( env_file_name, "r" );
  gzFile* pop_file = (gzFile*) gzopen( pop_file_name, "r" );
  gzFile* sp_struct_file  = NULL;
  
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
    sp_struct_file = (gzFile*) gzopen( sp_struct_file_name, "r" );
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

void ae_exp_manager::run_evolution( void )
{
  // dump the initial state of the population; useful for restarts
  _output_m->write_current_generation_outputs();
  
  while ( _num_gener < _last_gener )
  {
    printf( "============================== %"PRId32" ==============================\n", _num_gener );
    printf( "  distance to target (metabolic) : %f\n", ((ae_individual *) _pop->get_indivs()->get_last()->get_obj())->get_dist_to_target_by_feature( METABOLISM ) );

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
  printf( "  distance to target (metabolic) : %f\n", ((ae_individual *) _pop->get_indivs()->get_last()->get_obj())->get_dist_to_target_by_feature( METABOLISM ) );
  printf( "===================================================================\n");
  printf ("  The run is finished. \n"); 
  printf ("  Printing the final best individual into best_last_org.txt \n"); 
  const char* out_file_name = "best_last_org.txt"; 
  FILE* org_file = fopen( out_file_name, "w" );
  fputs( ((ae_individual *) _pop->get_indivs()->get_last()->get_obj())->get_genetic_unit(0)->get_dna()->get_data(), org_file); 
  fclose ( org_file ); 
}

// =================================================================
//                           Protected Methods
// =================================================================

// =================================================================
//                          Non inline accessors
// =================================================================
