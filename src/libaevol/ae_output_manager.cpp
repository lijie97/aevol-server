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


/*! \class ae_output_manager
    \brief
*/
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <sys/stat.h>
#include <err.h>
#include <errno.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_output_manager.h>
#include <ae_exp_manager.h>




//##############################################################################
//                                                                             #
//                           Class ae_output_manager                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_output_manager::ae_output_manager( ae_exp_manager* exp_m )
{
  _exp_m  = exp_m;
  _stats  = NULL;
  _tree   = NULL;
  _dump   = NULL;
}

// =================================================================
//                             Destructors
// =================================================================
ae_output_manager::~ae_output_manager( void )
{
  delete _stats;
  delete _tree;
  delete _dump;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_output_manager::save_experiment( void ) const
{
  // Generate backup file names
  char exp_setup_file_name[50];
  char pop_file_name[50];
  char out_man_file_name[50];
  
  int32_t num_gener = _exp_m->get_num_gener();
  
  sprintf( exp_setup_file_name, EXP_SETUP_BACKUP_FNAME_FORMAT,  num_gener );
  sprintf( pop_file_name,       POP_BACKUP_FNAME_FORMAT,        num_gener );
  sprintf( out_man_file_name,   OUT_PROF_BACKUP_FNAME_FORMAT,   num_gener );
  
  
  // Create missing backup directories
  int status = mkdir( POP_BACKUP_DIR, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, POP_BACKUP_DIR, errno );
  }
  status = mkdir( EXP_SETUP_BACKUP_DIR, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, EXP_SETUP_BACKUP_DIR, errno );
  }
  status = mkdir( OUT_PROF_BACKUP_DIR, 0755 );
  if ( status == -1 && errno != EEXIST )
  {
    err( EXIT_FAILURE, OUT_PROF_BACKUP_DIR, errno );
  }
  
  
  // Open backup files (population, exp_setup and output_man)
  gzFile* exp_setup_file  = (gzFile*) gzopen( exp_setup_file_name, "w" );
  gzFile* pop_file        = (gzFile*) gzopen( pop_file_name, "w" );
  gzFile* out_man_file    = (gzFile*) gzopen( out_man_file_name, "w" );
  if ( exp_setup_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, exp_setup_file_name );
    exit( EXIT_FAILURE );
  }
  if ( pop_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, pop_file_name );
    exit( EXIT_FAILURE );
  }
  if ( out_man_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, out_man_file_name );
    exit( EXIT_FAILURE );
  }
  
  // Write experimental setup, population data and output profile into the backups
  _exp_m->get_exp_s()->write_to_backup( exp_setup_file );
  _exp_m->get_pop()->write_to_backup( pop_file );
  this->write_to_backup( out_man_file );
  
  // Close backup files
  gzclose( exp_setup_file );
  gzclose( pop_file );
  gzclose( out_man_file );
}

void ae_output_manager::load_experiment(  char* exp_setup_file_name /*= NULL*/,
                                          char* pop_file_name /*= NULL*/,
                                          char* out_man_file_name /*= NULL*/,
                                          bool verbose /*= true*/ )
{
  //~ int32_t num_gener = _exp_m->get_num_gener();
  
  //~ // Generate default backup file names where needed
  //~ if ( exp_setup_file_name == NULL )
  //~ {
    //~ exp_setup_file_name = new char[50];
    //~ #ifdef __REGUL
      //~ sprintf( exp_setup_file_name, "backup/exp_setup_%06"PRId32".rae", num_gener );
    //~ #else
      //~ sprintf( exp_setup_file_name, "backup/exp_setup_%06"PRId32".ae", num_gener );
    //~ #endif
  //~ }
  //~ if ( exp_setup_file_name == NULL )
  //~ {
    //~ exp_setup_file_name = new char[50];
    //~ #ifdef __REGUL
      //~ sprintf( pop_file_name, "backup/pop_%06"PRId32".rae", num_gener );
    //~ #else
      //~ sprintf( pop_file_name, "backup/pop_%06"PRId32".ae", num_gener );
    //~ #endif
  //~ }
  //~ if ( exp_setup_file_name == NULL )
  //~ {
    //~ exp_setup_file_name = new char[50];
    //~ #ifdef __REGUL
      //~ sprintf( out_man_file_name, "backup/out_man_%06"PRId32".rae", num_gener );
    //~ #else
      //~ sprintf( out_man_file_name, "backup/out_man_%06"PRId32".ae", num_gener );
    //~ #endif
  //~ }
  
  // Open backup files (population, exp_setup and output_man)
  gzFile* exp_setup_file  = (gzFile*) gzopen( exp_setup_file_name, "r" );
  gzFile* pop_file        = (gzFile*) gzopen( pop_file_name, "r" );
  gzFile* out_man_file    = (gzFile*) gzopen( out_man_file_name, "r" );
  if ( exp_setup_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, exp_setup_file_name );
    exit( EXIT_FAILURE );
  }
  if ( pop_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, pop_file_name );
    exit( EXIT_FAILURE );
  }
  if ( out_man_file == Z_NULL )
  {
    printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, out_man_file_name );
    exit( EXIT_FAILURE );
  }
  
  // Load experimental setup, population data and output profile from the backups
  _exp_m->get_exp_s()->read_from_backup( exp_setup_file, verbose );
  _exp_m->get_pop()->read_from_backup( pop_file, verbose );
  this->read_from_backup( out_man_file, verbose );
  
  // Close backup files
  gzclose( exp_setup_file );
  gzclose( pop_file );
  gzclose( out_man_file );
}

void ae_output_manager::write_to_backup( gzFile* backup_file ) const
{
  // Write the backup steps
  gzwrite( backup_file, &_backup_step,      sizeof(_backup_step) );
  gzwrite( backup_file, &_big_backup_step,  sizeof(_big_backup_step) );
  
  // Stats
  
  // Tree
  int8_t record_tree = _record_tree;
  gzwrite( backup_file, &record_tree, sizeof(record_tree) );
  
  // Dumps
  int8_t make_dumps = _make_dumps;
  gzwrite( backup_file, &make_dumps,    sizeof(make_dumps) );
  gzwrite( backup_file, &_dump_period,  sizeof(_dump_period) );
}

void ae_output_manager::read_from_backup( gzFile* backup_file, bool verbose )
{
  // Write the backup steps
  gzread( backup_file, &_backup_step,      sizeof(_backup_step) );
  gzread( backup_file, &_big_backup_step,  sizeof(_big_backup_step) );
  
  // Stats
  
  // Tree
  int8_t record_tree;
  gzread( backup_file, &record_tree, sizeof(record_tree) );
  _record_tree = record_tree;
  
  // Dumps
  int8_t make_dumps;
  gzread( backup_file, &make_dumps, sizeof(make_dumps) );
  _make_dumps = make_dumps;
  gzread( backup_file, &_dump_period,  sizeof(_dump_period) );
}

void ae_output_manager::write_current_generation_outputs( void ) const
{
  int32_t num_gener = _exp_m->get_num_gener();
  
  _stats->write_current_generation_statistics();
  _dump->write_current_generation_dump();
  
  if ( _record_tree )
  { 
    _tree->fill_tree_with_cur_gener(); 
  }

  // Write backup and tree
  if ( num_gener % _tree->get_tree_step() == 0 )    
  {
    if ( _record_tree  && _tree->get_tree_mode() == NORMAL ) 
    { 
      write_tree();
    }
  }
    
  if ( num_gener % _backup_step == 0 )
  {
    _stats->flush();
    save_experiment();
    
    // Update the last_gener.txt file
    FILE* last_gener_file = fopen( "last_gener.txt", "w" );
    if ( last_gener_file != NULL )
    {
      fprintf( last_gener_file, "%"PRId32"\n", num_gener );
      fclose( last_gener_file );
    }
    else
    {
      printf( "Error : could not open file last_gener.txt\n" );
    }
    
    #ifdef __IN2P3
      if ( num_gener % _big_backup_step == 0 )
      {
        system( "./make_big_backup.py" );
      }
    #endif
  }

  if ( _make_dumps ) 
  {
    if( num_gener % _dump_period == 0 )
    {
      _dump->write_current_generation_dump();
    }
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
void ae_output_manager::write_tree( void ) const
{
  char tree_file_name[50];
  
#ifdef __REGUL
  sprintf( tree_file_name, "tree/tree_%06"PRId32".rae", _exp_m->get_num_gener() );
#else
  sprintf( tree_file_name, "tree/tree_%06"PRId32".ae", _exp_m->get_num_gener() );
#endif
  
  gzFile* tree_file = (gzFile*) gzopen( tree_file_name, "w" );
  
  // Write phylogenetic data (tree)
  _tree->write_to_tree_file( tree_file );
  
  gzclose( tree_file );
}

// =================================================================
//                          Non inline accessors
// =================================================================
