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


/*! \class ae_simulation
    \brief This class contains all the data regarding a simulation (either main run or post-treatment)
*/


// =================================================================
//                              Libraries
// =================================================================
#include <stdio.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>




// =================================================================
//                            Project Files
// =================================================================
#ifdef __NO_X
#include <ae_simulation.h>
#include <ae_population.h>
#elif defined __X11
#include <ae_simulation_X11.h>
#include <ae_population_X11.h>
#endif

#include <ae_param_loader.h>
#include <ae_param_overloader.h>
#include <ae_simulation_X11.h>




//##############################################################################
//                                                                             #
//                             Class ae_simulation                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_simulation::ae_simulation( void )
{
  _num_gener    = 0;
  _first_gener  = 0;
  _last_gener   = 0;

  _env  = NULL;
  _pop  = NULL;
  
  _stats  = NULL;
  _tree   = NULL;
  _dump   = NULL;
}
  
/*ae_simulation::ae_simulation( ae_param_overloader* param_overloader )
{
  // Some calls in this constructor will need to access the simulation (i.e. <this>)
  // through <ae_common::sim> which would normaly be initialized at the end of the construction.
  //
  // In order to keep the access to the simulation homogeneous within the entire application,
  // we will force the affectation now.
  //
  // Note that this is done by obligation here, it must not be considered good practice
  
  ae_common::sim = this;  

  // Initialize "trivial" attributes
  _first_gener  = 0;
  _num_gener    = 0;

  // Load parameter file
  ae_param_loader* params = new ae_param_loader();
  params->load();
  if ( param_overloader != NULL )
  {
    param_overloader->overload_params();
  }
  
  // Create a pseudo-random number generator
  alea = new ae_rand_mt( ae_common::init_params->get_seed() );
  
  #ifdef __REGUL
    // Initialisation of the evaluation dates
    if ( ae_common::params->get_individual_evaluation_dates() == NULL )
    {
      ae_common::individual_evaluation_dates = new ae_array_short( ae_common::individual_evaluation_nbr );
      ae_common::individual_evaluation_dates->set_value( 0, 10 );
    }
    ae_common::init_binding_matrix();
  #endif
  
  delete params;

  // Create the environment
  _env = new ae_environment();

 // Create new population
  #ifdef __NO_X
    _pop = new ae_population();
  #elif defined __X11
    _pop = new ae_population_X11();
  #endif

  // Create statistics files
  _stats = new ae_stats();
  _stats->write_headers();
  _stats->write_current_generation_statistics();

  // Create backup directory and write backup 
  mkdir( "backup", 0755 );
  write_backup();
  
  // Create dump directory and init dump
  _dump = NULL;
  if ( ae_common::rec_params->get_dump_period() > 0 )
  {
    mkdir( "dump", 0755 );
    _dump = new ae_dump();
  }
  
  // Create tree directory and init tree
  _tree = NULL;
  if ( ae_common::rec_params->get_record_tree() == true )
  { 
    mkdir( "tree", 0755 );
    _tree = new ae_tree();
  }
  
  
  // Write an entry in the LOADS log file
  if ( ae_common::rec_params->is_logged( LOG_LOADS ) == true )
  {
    if ( param_overloader != NULL )
    {
      param_overloader->write_log( ae_common::rec_params->get_log( LOG_LOADS ) );
    }
    else
    {
      fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "  No overloaded parameters\n\n" );
    }
  }
}*/

/*ae_simulation::ae_simulation( char* organism_file_name, ae_param_overloader* param_overloader )
{
  // Some calls in this constructor will need to access the simulation (i.e. <this>)
  // through <ae_common::sim> which would normaly be initialized at the end of the construction.
  //
  // In order to keep the access to the simulation homogeneous within the entire application,
  // we will force the affectation now.
  //
  // Note that this is done by obligation here, it must not be considered good practice
  
  ae_common::sim = this;  

  // Initialize "trivial" attributes
  _first_gener  = 0;
  _num_gener    = 0;

  // Load parameter file
  ae_param_loader* params = new ae_param_loader();
  params->load();
  if ( param_overloader != NULL )
  {
    param_overloader->overload_params();
  }
  
  // Create a pseudo-random number generator
  alea = new ae_rand_mt( ae_common::init_params->get_seed() );
  
  #ifdef __REGUL
  // Initialisation of the evaluation dates
  if ( ae_common::individual_evaluation_dates == NULL )
  {
    ae_common::individual_evaluation_dates = new ae_array_short( ae_common::individual_evaluation_nbr );
    ae_common::individual_evaluation_dates->set_value( 0, 10 );
  }
  ae_common::init_binding_matrix();
  #endif
  
  delete params;

  // Create the environment
  _env = new ae_environment();

 // Create new population
  #ifdef __NO_X
    _pop = new ae_population( organism_file_name );
  #elif defined __X11
    printf( "ERROR : starting X11 runs from single organism is not yet implemented!" );
    _pop = new ae_population_X11();
  #endif

  // Create statistics files
  _stats = new ae_stats();
  _stats->write_headers();
  _stats->write_current_generation_statistics();

  // Create backup directory and write backup 
  mkdir( "backup", 0755 );
  write_backup();
  
  // Create dump directory and init dump
  _dump = NULL;
  if ( ae_common::rec_params->get_dump_period() > 0 )
  {
    mkdir( "dump", 0755 );
    _dump = new ae_dump();
  }
  
  // Create tree directory and init tree
  _tree = NULL;
  if ( ae_common::rec_params->get_record_tree() == true )
  { 
    mkdir( "tree", 0755 );
    _tree = new ae_tree();
  }
  
  
  // Write an entry in the LOADS log file
  if ( ae_common::rec_params->is_logged( LOG_LOADS ) == true )
  {
    fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "Simulation lauched\n" );
    if ( param_overloader != NULL )
    {
      fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "  Overloaded parameters:\n" );
      param_overloader->write_log( ae_common::rec_params->get_log( LOG_LOADS ) );
      fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "\n" );
    }
    else
    {
      fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "  No overloaded parameters\n\n" );
    }
  }
}*/


//ae_simulation::ae_simulation( char* backup_file_name, bool to_be_run /* = TRUE */, ae_param_overloader* param_overloader /* = NULL */ )
/*{
  // Some calls in this constructor will need to access the simulation (i.e. <this>)
  // through <ae_common::sim> which would normaly be initialized at the end of the construction.
  //
  // In order to keep the access to the simulation homogeneous within the entire application,
  // we will force the affectation now.
  //
  // Note that this is done by obligation here, it must not be considered good practice
  ae_common::sim = this;
  
  // Check whether the backup file is an aevol or a raevol backup
  int16_t bfn_len = strlen( backup_file_name );
  #ifdef __REGUL
    if ( strcmp( &backup_file_name[bfn_len-4], ".rae" ) != 0 )
    {
      printf( "ERROR : %s is not valid RAEVOL backup file.\n", backup_file_name );
      exit( EXIT_FAILURE );
    }
  #else
    if ( strcmp( &backup_file_name[bfn_len-3], ".ae" ) != 0 || strstr(backup_file_name,"best") !=NULL)
    {
      printf( "ERROR : %s is not valid AEVOL backup file.\n", backup_file_name );
      exit( EXIT_FAILURE );
    }
  #endif
  

  gzFile* backup_file = (gzFile*) gzopen( backup_file_name, "r" );
  
  if ( backup_file == Z_NULL )
  {
    printf( "ERROR : Could not read backup file %s\n", backup_file_name );
    exit( EXIT_FAILURE );
  }

  // Retrieve random generator state
  printf( "  Loading random generator\n" );
  alea = new ae_rand_mt( backup_file );

  // Retreive common data
  printf( "  Loading common data\n" );
  ae_common::read_from_backup( backup_file );
  if ( param_overloader != NULL )
  {
    int32_t seed_from_backup = ae_common::init_params->get_seed();
    param_overloader->overload_params();
    if ( ae_common::init_params->get_seed() != seed_from_backup )
    {
      // We need to reseed
      delete alea;
      alea = new ae_rand_mt( ae_common::init_params->get_seed() );
    }
  }

  // Retrieve current generation number
  printf( "  Loading simulation data\n" );
  gzread( backup_file, &_num_gener, sizeof(_num_gener) );
  printf( "    _num_gener : %d\n", _num_gener );
  _first_gener = _num_gener;

  // Retrieve environmental data
  printf( "  Loading environment\n" );
  _env = new ae_environment( backup_file );

  // Retrieve population data and individuals
  printf( "  Loading population\n" );
#ifdef __NO_X
  //~ printf( "new ae_population( backup_file );\n" );
  _pop = new ae_population( backup_file );
#elif defined __X11
  //~ printf( "new ae_population_X11( backup_file );\n" );
  _pop = new ae_population_X11( backup_file );
#endif

  gzclose( backup_file );

  // Evaluate individuals
  _pop->evaluate_individuals( _env );

  // if the population is spatially structured, then the individuals have been saved and loaded in the order of the grid and not in increasing order of fitness
  // so we have to sort the individuals

  if ( ae_common::pop_structure == true )
  {
    _pop->sort_individuals();
  }


  if ( to_be_run ) // We want to restart the run from here (not just do a post-processing)
  {
    // Prepare stat files
    _stats  = new ae_stats( _first_gener );
    
    // Prepare tree
    if ( ae_common::rec_params->get_record_tree() == true )
    { 
      mkdir( "tree", 0755 );
      _tree  = new ae_tree(); 
    }
    else
    {
      _tree = NULL;
    }
    
    // Prepare dump
    if ( ae_common::rec_params->get_dump_period() > 0 )
    {
      mkdir( "dump", 0755 );
      _dump = new ae_dump();
    }
    else
    {
      _dump = NULL;
    }
    
    mkdir( "backup", 0755 );
    
    if ( ae_common::rec_params->is_logged( LOG_LOADS ) == true )
    {
      // Write an entry in the LOADS log file
      fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "GENERATION_OVERLOAD %"PRId32"\n", _num_gener );
      if ( param_overloader->get_nb_overloaded() > 0 )
      {
        //fprintf( _logs->get_log( LOG_LOADS ), "  Overloaded parameters:\n" );
        param_overloader->write_log( ae_common::rec_params->get_log( LOG_LOADS ) );
        fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "\n" );
      }
      else
      {
        fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "  No overloaded parameters\n\n" );
      }
    }
  }
  else
  {
    // We just want to inspect the state of the simulation at this moment
    ae_common::rec_params->init_logs( 0 );
    _stats = NULL;
    _tree  = NULL; 
    _dump  = NULL;
  }
}*/

// =================================================================
//                             Destructors
// =================================================================
ae_simulation::~ae_simulation( void )
{
  if (ae_common::rec_params->get_record_tree() == true) { delete _tree; } 
  if (ae_common::rec_params->get_dump_period() > 0) { delete _dump; }

  delete _stats;
  delete _pop;
  delete _env;
  
  delete alea;
}

// =================================================================
//                            Public Methods
// =================================================================
/*!
  Create simulation according to parameters from a param file and pptionally parameter overloads
*/ 
void ae_simulation::load_params( ae_param_loader* param_loader, ae_param_overloader* param_overloader /*= NULL*/ )
{
  // Load parameters
  param_loader->load();
  if ( param_overloader != NULL )
  {
    param_overloader->overload_params();
  }
  
  // Create a pseudo-random number generator
  alea = new ae_rand_mt( ae_common::init_params->get_seed() );
  
  #ifdef __REGUL
    // Initialisation of the evaluation dates
    if ( ae_common::params->get_individual_evaluation_dates() == NULL )
    {
      ae_common::individual_evaluation_dates = new ae_array_short( ae_common::individual_evaluation_nbr );
      ae_common::individual_evaluation_dates->set_value( 0, 10 );
    }
    ae_common::init_binding_matrix();
  #endif

  // Create the environment
  _env = new ae_environment();

 // Create new population
  #ifdef __NO_X
    _pop = new ae_population();
  #elif defined __X11
    _pop = new ae_population_X11();
  #endif

  // Create statistics files
  _stats = new ae_stats();
  _stats->write_headers();
  _stats->write_current_generation_statistics();

  // Create backup directory and write backup 
  mkdir( "backup", 0755 );
  write_backup();
  
  // Create dump directory and init dump
  _dump = NULL;
  if ( ae_common::rec_params->get_dump_period() > 0 )
  {
    mkdir( "dump", 0755 );
    _dump = new ae_dump();
  }
  
  // Create tree directory and init tree
  _tree = NULL;
  if ( ae_common::rec_params->get_record_tree() == true )
  { 
    mkdir( "tree", 0755 );
    _tree = new ae_tree();
  }
  
  
  // Write an entry in the LOADS log file
  if ( ae_common::rec_params->is_logged( LOG_LOADS ) == true )
  {
    if ( param_overloader != NULL )
    {
      param_overloader->write_log( ae_common::rec_params->get_log( LOG_LOADS ) );
    }
    else
    {
      fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "  No overloaded parameters\n\n" );
    }
  }
  
  // Initialize display
  #ifdef __X11
    ((ae_simulation_X11*) this)->initialize( ae_common::pop_structure, ae_common::params->get_allow_plasmids() );
  #endif // def __X11
}

/*!
  Load the simulation from a backup. Optionally overload parameters and modify simulation accordingly
*/ 
void ae_simulation::load_backup( char* backup_file_name, bool to_be_run /*= true*/, ae_param_overloader* param_overloader /*= NULL*/ )
{
  // Open backup file
  gzFile* backup_file = (gzFile*) gzopen( backup_file_name, "r" );
  if ( backup_file == Z_NULL )
  {
    printf( "ERROR : Could not read backup file %s\n", backup_file_name );
    exit( EXIT_FAILURE );
  }
  
  // Retrieve data from the backup
  this->read_from_backup( backup_file );
  
  // Close backup file
  gzclose( backup_file );
  
  // Overload parameters
  if ( param_overloader != NULL )
  {
    // TODO if ( param_overloader->is_overloadded( SEED ) ) ... !!!!!
    int32_t seed_from_backup = ae_common::init_params->get_seed();
    param_overloader->overload_params();
    if ( ae_common::init_params->get_seed() != seed_from_backup )
    {
      // We need to reseed
      delete alea;
      alea = new ae_rand_mt( ae_common::init_params->get_seed() );
    }
  }

  // Evaluate individuals
  _pop->evaluate_individuals( _env );

  // If the population is spatially structured, then the individuals are saved and loaded in the order of the grid and not in increasing order of fitness
  // so we have to sort the individuals
  if ( ae_common::pop_structure == true )
  {
    _pop->sort_individuals();
  }


  // If the simulation is being continued (not just do a post-processed), prepare output data accordingly
  if ( to_be_run )
  {
    // Prepare stat files
    _stats  = new ae_stats( _first_gener );
    
    // Prepare tree
    if ( ae_common::rec_params->get_record_tree() == true )
    { 
      mkdir( "tree", 0755 );
      _tree  = new ae_tree(); 
    }
    else
    {
      _tree = NULL;
    }
    
    // Prepare dump
    if ( ae_common::rec_params->get_dump_period() > 0 )
    {
      mkdir( "dump", 0755 );
      _dump = new ae_dump();
    }
    else
    {
      _dump = NULL;
    }
    
    mkdir( "backup", 0755 );
    
    if ( ae_common::rec_params->is_logged( LOG_LOADS ) == true )
    {
      // Write an entry in the LOADS log file
      fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "GENERATION_OVERLOAD %"PRId32"\n", _num_gener );
      if ( param_overloader->get_nb_overloaded() > 0 )
      {
        //fprintf( _logs->get_log( LOG_LOADS ), "  Overloaded parameters:\n" );
        param_overloader->write_log( ae_common::rec_params->get_log( LOG_LOADS ) );
        fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "\n" );
      }
      else
      {
        fprintf( ae_common::rec_params->get_log( LOG_LOADS ), "  No overloaded parameters\n\n" );
      }
    }
  }
  else
  {
    // We just want to inspect the state of the simulation at this moment
    ae_common::rec_params->init_logs( 0 );
  }
  
  // Initialize display
  #ifdef __X11
    ((ae_simulation_X11*) this)->initialize( ae_common::pop_structure, ae_common::params->get_allow_plasmids() );
  #endif // def __X11
}

void ae_simulation::run( void )
{
  _last_gener = _first_gener + ae_common::nb_generations;

  // dump the initial state of the population; useful for restarts
  if( ae_common::rec_params->get_dump_period() > 0 ) 
  {
    _dump->write_current_generation_dump();
  }
  
  while ( _num_gener < _last_gener )
  {
    printf( "============================== %"PRId32" ==============================\n", _num_gener );
    printf( "  distance to target (metabolic) : %f\n", ((ae_individual *) _pop->get_indivs()->get_last()->get_obj())->get_dist_to_target_by_feature( METABOLISM ) );

    #ifdef __X11
      if ( ((ae_simulation_X11*)this)->quit_signal_received() ) break;
      ((ae_simulation_X11*)this)->display();
    #endif
    
    // Environmental variation
    _env->apply_variation();

    // Create next generation
    if ( ae_common::pop_structure == false )  { _pop->step_to_next_generation();      } 
    else                                      { _pop->step_to_next_generation_grid(); } 
    
    _num_gener++;

    // Write statistical data and store phylogenetic data (tree)
    _stats->write_current_generation_statistics();

    if ( ae_common::rec_params->get_record_tree() == true )
    { 
      _tree->fill_tree_with_cur_gener(); 
    }

    // Write backup and tree
    if ( _num_gener % ae_common::rec_params->get_tree_step() == 0 )    
    {
      if ( ae_common::rec_params->get_record_tree() == true  && ae_common::rec_params->get_tree_mode() == NORMAL ) 
      { 
        write_tree();
      }
    }
      
    if ( _num_gener % ae_common::rec_params->get_backup_step() == 0 )
    {
      _stats->flush();
      write_backup();
      
      // Update the last_gener.txt file
      FILE* last_gener_file = fopen( "last_gener.txt", "w" );
      if ( last_gener_file != NULL )
      {
        fprintf( last_gener_file, "%"PRId32"\n", _num_gener );
        fclose( last_gener_file );
      }
      else
      {
        printf( "Error : could not open file last_gener.txt\n" );
      }
      
      #ifdef __IN2P3
        if ( _num_gener % ae_common::big_backup_step == 0 )
        {
          system( "./make_big_backup.py" );
        }
      #endif
    }

    if( ae_common::rec_params->get_dump_period() > 0 ) 
    {
      if( _num_gener % ae_common::rec_params->get_dump_period() == 0 )
      {
        _dump->write_current_generation_dump();
      }
    }
  }
  
  _stats->flush();
  
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

void ae_simulation::set_total_generations( int32_t nb_generations )
{
  _last_gener = _first_gener + nb_generations;
}


void ae_simulation::write_backup( void )
{
  char backup_file_name[50];
  char best_indiv_file_name[50];
  
#ifdef __REGUL
  sprintf( backup_file_name, "backup/gen_%06"PRId32".rae", _num_gener );
  sprintf( best_indiv_file_name, "backup/best_%06"PRId32".rae", _num_gener );
#else
  sprintf( backup_file_name, "backup/gen_%06"PRId32".ae", _num_gener );
  sprintf( best_indiv_file_name, "backup/best_%06"PRId32".ae", _num_gener );
#endif
  
  gzFile* backup_file     = (gzFile*) gzopen( backup_file_name, "w" );
  gzFile* best_indiv_file = (gzFile*) gzopen( best_indiv_file_name, "w" );

  // Write random generator state
  alea->write_to_backup( backup_file );

  // Write common data
  ae_common::write_to_backup( backup_file );

  // Write simulation intrinsic data
  gzwrite( backup_file, &_num_gener, sizeof(_num_gener) );

  // Write environmental data
  _env->write_to_backup( backup_file );

  // Write population data and individuals
  _pop->write_to_backup( backup_file );
  
  // Write common data and best individual in best_indiv_file
  ae_common::write_to_backup( best_indiv_file );
  _pop->get_best()->write_to_backup( best_indiv_file );
  
  gzclose( backup_file );
  gzclose( best_indiv_file );
}

void ae_simulation::write_tree( void )
{
  char tree_file_name[50];
  
#ifdef __REGUL
  sprintf( tree_file_name, "tree/tree_%06"PRId32".rae", _num_gener );
#else
  sprintf( tree_file_name, "tree/tree_%06"PRId32".ae", _num_gener );
#endif
  
  gzFile* tree_file = (gzFile*) gzopen( tree_file_name, "w" );
  
  // Write phylogenetic data (tree)
  _tree->write_to_tree_file( tree_file );
  
  gzclose( tree_file );
}




// =================================================================
//                           Protected Methods
// =================================================================
void ae_simulation::read_from_backup( gzFile* backup_file )
{
  // Retrieve random generator state
  printf( "  Loading random generator..." );
  fflush( stdout );
  alea = new ae_rand_mt( backup_file );
  printf( "OK\n" );

  // Retreive common data
  printf( "  Loading common data..." );
  fflush( stdout );
  ae_common::read_from_backup( backup_file );
  printf( "OK\n" );

  // Retrieve current generation number
  printf( "  Loading simulation data..." );
  fflush( stdout );
  gzread( backup_file, &_num_gener, sizeof(_num_gener) );
  _first_gener = _num_gener;
  printf( "OK\n" );

  // Retrieve environmental data
  printf( "  Loading environment..." );
  fflush( stdout );
  _env = new ae_environment( backup_file );
  printf( "OK\n" );

  // Retrieve population data and individuals
  printf( "  Loading population..." );
  fflush( stdout );
  #ifdef __NO_X
    //~ printf( "new ae_population( backup_file );\n" );
    _pop = new ae_population( backup_file );
  #elif defined __X11
    //~ printf( "new ae_population_X11( backup_file );\n" );
    _pop = new ae_population_X11( backup_file );
  #endif
  printf( "OK\n" );
}
