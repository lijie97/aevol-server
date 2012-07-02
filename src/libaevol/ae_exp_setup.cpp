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


/*! \class ae_exp_setup
    \brief This class contains the experiment setup among which the environment
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
#include <ae_exp_setup.h>
#include <ae_population.h>
#elif defined __X11
#include <ae_exp_setup_X11.h>
#include <ae_population_X11.h>
#endif




//##############################################################################
//                                                                             #
//                              Class ae_exp_setup                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_exp_setup::ae_exp_setup( ae_exp_manager* exp_m )
{
  _exp_m  = exp_m;
  _env    = new ae_environment();
  _sel    = new ae_selection( exp_m );
}
  
/*ae_exp_setup::ae_exp_setup( ae_param_overloader* param_overloader )
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

/*ae_exp_setup::ae_exp_setup( char* organism_file_name, ae_param_overloader* param_overloader )
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


//ae_exp_setup::ae_exp_setup( char* backup_file_name, bool to_be_run /* = TRUE */, ae_param_overloader* param_overloader /* = NULL */ )
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
ae_exp_setup::~ae_exp_setup( void )
{
  delete _env;
}

// =================================================================
//                            Public Methods
// =================================================================
/*!
  Create simulation according to parameters from a param file and pptionally parameter overloads
*/ 
//void ae_exp_setup::create_from_param_file( ae_param_loader* param_loader, ae_param_overloader* param_overloader /*= NULL*/ )
/*{ 
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
    ((ae_exp_setup_X11*) this)->initialize( ae_common::pop_structure, ae_common::params->get_allow_plasmids() );
  #endif // def __X11
}*/

/*!
  Load the simulation from a backup. Optionally overload parameters and modify simulation accordingly
*/ 











// =================================================================
//                           Protected Methods
// =================================================================
void ae_exp_setup::write_to_backup( gzFile* backup_file ) const
{
  // Write environmental data
  _env->write_to_backup( backup_file );
}

void ae_exp_setup::read_from_backup( gzFile* backup_file, bool verbose )
{
  //~ // Retreive common data
  //~ printf( "  Loading common data..." );
  //~ fflush( stdout );
  //~ ae_common::read_from_backup( backup_file );
  //~ printf( "OK\n" );

  //~ // Retrieve current generation number
  //~ printf( "  Loading simulation data..." );
  //~ fflush( stdout );
  //~ gzread( backup_file, &_num_gener, sizeof(_num_gener) );
  //~ _first_gener = _num_gener;
  //~ printf( "OK\n" );

  // Retrieve environmental data
  printf( "  Loading environment..." );
  fflush( stdout );
  _env = new ae_environment();
  _env->read_from_backup( backup_file );
  printf( "OK\n" );
}
