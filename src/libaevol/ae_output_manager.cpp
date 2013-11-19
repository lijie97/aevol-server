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
  _compute_phen_contrib_by_GU = false;
  _record_tree = false;
  _make_dumps = false;
  _dump_step = 0;
  _logs  = new ae_logs();
}

// =================================================================
//                             Destructors
// =================================================================
ae_output_manager::~ae_output_manager( void )
{
  delete _stats;
  delete _tree;
  delete _dump;
  delete _logs;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_output_manager::write_setup_file( gzFile setup_file ) const
{
  // Write the backup steps
  gzwrite( setup_file, &_backup_step,      sizeof(_backup_step) );
  gzwrite( setup_file, &_big_backup_step,  sizeof(_big_backup_step) );
  
  // Stats
  gzwrite( setup_file, &_compute_phen_contrib_by_GU,  sizeof(_compute_phen_contrib_by_GU) );
  
  // Tree
  int8_t record_tree = _record_tree;
  gzwrite( setup_file, &record_tree, sizeof(record_tree) );
  if ( _record_tree )
  {
    int32_t tmp_tree_step = _tree->get_tree_step();
    gzwrite( setup_file, &tmp_tree_step, sizeof(tmp_tree_step) );
    int8_t tmp_tree_mode = (int8_t) _tree->get_tree_mode();
    gzwrite( setup_file, &tmp_tree_mode, sizeof(tmp_tree_mode) );
  }
  
  // Dumps
  int8_t make_dumps = _make_dumps;
  gzwrite( setup_file, &make_dumps,  sizeof(make_dumps) );
  gzwrite( setup_file, &_dump_step,  sizeof(_dump_step) );
  
  // Logs
  int8_t logs = _logs->get_logs();
  gzwrite( setup_file, &logs,  sizeof(logs) );
}

void ae_output_manager::write_setup_file( FILE* setup_file ) const
{
  // Write the backup steps
  fprintf( setup_file, "BACKUP_STEP %"PRId32"\n", _backup_step );
  fprintf( setup_file, "BIG_BACKUP_STEP %"PRId32"\n", _big_backup_step );
  
  // Stats
  fprintf( setup_file, "COMPUTE_PHENOTYPIC_CONTRIBUTION_BY_GU %"PRId8"\n", _compute_phen_contrib_by_GU );
  
  // Tree
  fprintf( setup_file, "RECORD_TREE %s\n", _record_tree ? "true" : "false" );
  if ( _record_tree )
  {
    fprintf( setup_file, "TREE_STEP %"PRId32"\n", _tree->get_tree_step() );
    
    if ( _tree->get_tree_mode() == LIGHT )
    {
      fprintf( setup_file, "TREE_MODE LIGHT\n" );
    }
    else if ( _tree->get_tree_mode() == NORMAL )
    {
      fprintf( setup_file, "TREE_MODE NORMAL\n" );
    }
    else
    {
      fprintf( setup_file, "TREE_MODE UNKNOWN\n" );
    }
  }
  
  // Dumps
  fprintf( setup_file, "MAKE_DUMPS %s\n", _make_dumps ? "true" : "false" );
  fprintf( setup_file, "DUMP_STEP %"PRId32"\n", _dump_step );
  
  // Logs
  int8_t logs = _logs->get_logs();
  fprintf( setup_file, "LOGS %"PRId8"\n", logs );
}

void ae_output_manager::load( gzFile setup_file, bool verbose, bool to_be_run  )
{
  // Write the backup steps
  gzread( setup_file, &_backup_step,      sizeof(_backup_step) );
  gzread( setup_file, &_big_backup_step,  sizeof(_big_backup_step) );
  
  // Stats
  int32_t num_gener = _exp_m->get_num_gener();
  if( to_be_run)
  {
    if ( num_gener > 0 )
    {
      _stats = new ae_stats( _exp_m, num_gener );
    }
    else
    {
      _stats = new ae_stats( _exp_m );
    }
  }
  gzread( setup_file, &_compute_phen_contrib_by_GU,  sizeof(_compute_phen_contrib_by_GU) );
  
  // Tree
  int8_t record_tree;
  gzread( setup_file, &record_tree, sizeof(record_tree) );
  _record_tree = record_tree;
  if ( _record_tree )
  {
    int32_t tmp_tree_step;
    gzread( setup_file, &tmp_tree_step, sizeof(tmp_tree_step) );
    int8_t tmp_tree_mode;
    gzread( setup_file, &tmp_tree_mode, sizeof(tmp_tree_mode) );
    if ( (ae_tree_mode)tmp_tree_mode != LIGHT && (ae_tree_mode)tmp_tree_mode != NORMAL)
    {
        printf( "%s:%d: error: invalid tree mode\n", __FILE__, __LINE__ );
        exit( EXIT_FAILURE );
    }
    
    _tree = new ae_tree( _exp_m, (ae_tree_mode) tmp_tree_mode, tmp_tree_step );
  }
  
  // Dumps
  int8_t make_dumps;
  gzread( setup_file, &make_dumps, sizeof(make_dumps) );
  _make_dumps = make_dumps;
  gzread( setup_file, &_dump_step,  sizeof(_dump_step) );
  if( _make_dumps == true)
  {
    _dump = new ae_dump(_exp_m);
  }
  
  // Logs
  int8_t logs;
  gzread( setup_file, &logs, sizeof(logs) );
  if( to_be_run)
  {
    _logs->load(logs, num_gener);
  }
}

void ae_output_manager::load( FILE* setup_file, bool verbose, bool to_be_run  )
{
  // Write the backup steps
  fscanf( setup_file, "BACKUP_STEP %"PRId32"\n", &_backup_step );
  fscanf( setup_file, "BIG_BACKUP_STEP %"PRId32"\n", &_big_backup_step );
  
  // Stats
  int32_t num_gener = _exp_m->get_num_gener();
  if( to_be_run)
  {
    if ( num_gener > 0 )
    {
      _stats = new ae_stats( _exp_m, num_gener );
    }
    else
    {
      _stats = new ae_stats( _exp_m );
    }
  }
  fscanf( setup_file, "COMPUTE_PHENOTYPIC_CONTRIBUTION_BY_GU %"PRId8"\n", &_compute_phen_contrib_by_GU );
  
  char tmp[10];
  
  // Tree
  fscanf( setup_file, "RECORD_TREE %s\n", tmp );
  _record_tree = ! strcmp( tmp, "true" );
  if ( _record_tree )
  {
    int32_t tmp_tree_step;
    fscanf( setup_file, "TREE_STEP %"PRId32"\n", &tmp_tree_step );
    int8_t tmp_tree_mode;
    fscanf( setup_file, "TREE_MODE %"PRId8"\n", &tmp_tree_mode );
    if ( (ae_tree_mode)tmp_tree_mode != LIGHT && (ae_tree_mode)tmp_tree_mode != NORMAL)
    {
        printf( "%s:%d: error: invalid tree mode\n", __FILE__, __LINE__ );
        assert( false );
        exit( EXIT_FAILURE );
    }
    
    _tree = new ae_tree( _exp_m, (ae_tree_mode) tmp_tree_mode, tmp_tree_step );
  }
  
  // Dumps
  fscanf( setup_file, "MAKE_DUMPS %s\n", tmp );
  _make_dumps = ! strcmp( tmp, "true" );
  fscanf( setup_file, "DUMP_STEP %"PRId32"\n", &_dump_step );
  if( _make_dumps == true)
  {
    _dump = new ae_dump(_exp_m);
  }
  
  // Logs
  int8_t logs;
  fscanf( setup_file, "LOGS %"PRId8"\n", &logs );
  _logs->load(logs, num_gener);
}

void ae_output_manager::write_current_generation_outputs( void ) const
{
  int32_t num_gener = _exp_m->get_num_gener();
  
  _stats->write_current_generation_statistics();
  
  if ( _record_tree )
  {
    if ( num_gener > 0 )
    {
      _tree->fill_tree_with_cur_gener(); 
    }
  }

  // Write backup and tree
  if ( _record_tree && (num_gener != _exp_m->get_first_gener()) && (num_gener % _tree->get_tree_step() == 0) )    
  {
    if ( _tree->get_tree_mode() == NORMAL ) 
    { 
      write_tree();
    }
  }
  
  if ( num_gener % _backup_step == 0 )
  {
    _stats->flush();
    _exp_m->save();
    
    // Update the LAST_GENER file
    FILE* last_gener_file = fopen( LAST_GENER_FNAME, "w" );
    if ( last_gener_file != NULL )
    {
      fprintf( last_gener_file, "%"PRId32"\n", num_gener );
      fclose( last_gener_file );
    }
    else
    {
      printf( "Error : could not open file "LAST_GENER_FNAME"\n" );
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
    if( num_gener % _dump_step == 0 )
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
  // Create the tree directory if it doesn't exist
  int status;
  status = mkdir( TREE_DIR, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, TREE_DIR, errno );
  }
  
  char tree_file_name[50];
  
#ifdef __REGUL
  sprintf( tree_file_name, "tree/tree_%06"PRId32".rae", _exp_m->get_num_gener() );
#else
  sprintf( tree_file_name, "tree/tree_%06"PRId32".ae", _exp_m->get_num_gener() );
#endif
  
  gzFile tree_file = gzopen( tree_file_name, "w" );
  
  // Write phylogenetic data (tree)
  _tree->write_to_tree_file( tree_file );
  
  gzclose( tree_file );
}

// =================================================================
//                          Non inline accessors
// =================================================================
