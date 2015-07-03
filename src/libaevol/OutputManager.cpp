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
#include "OutputManager.h"
#include "ExpManager.h"
#include "Time.h"

namespace aevol {



//##############################################################################
//                                                                             #
//                           Class OutputManager                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
OutputManager::OutputManager(ExpManager * exp_m)
{
  _exp_m  = exp_m;
  _stats  = new Stats(exp_m);;
  _tree   = NULL;
  _dump   = NULL;
  _compute_phen_contrib_by_GU = false;
  _record_tree = false;
  _make_dumps = false;
  _dump_step = 0;
  _logs  = new Logging();
}

// =================================================================
//                             Destructors
// =================================================================
OutputManager::~OutputManager( void )
{
  delete _stats;
  delete _tree;
  delete _dump;
  delete _logs;
}

// =================================================================
//                            Public Methods
// =================================================================
void OutputManager::write_setup_file( gzFile setup_file ) const
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
    auto tmp_tree_step = _tree->get_tree_step();
    gzwrite( setup_file, &tmp_tree_step, sizeof(tmp_tree_step) );
  }
  
  // Dumps
  int8_t make_dumps = _make_dumps;
  gzwrite( setup_file, &make_dumps,  sizeof(make_dumps) );
  gzwrite( setup_file, &_dump_step,  sizeof(_dump_step) );
  
  // Logs
  int8_t logs = _logs->get_logs();
  gzwrite( setup_file, &logs,  sizeof(logs) );
}

void OutputManager::write_setup_file(FILE* setup_file) const
{
  // Write the backup steps
  fprintf(setup_file, "BACKUP_STEP %" PRId64 "\n", _backup_step);
  fprintf(setup_file, "BIG_BACKUP_STEP %" PRId64 "\n", _big_backup_step);
  
  // Stats
  fprintf(setup_file, "COMPUTE_PHENOTYPIC_CONTRIBUTION_BY_GU %" PRId8 "\n", (int8_t) _compute_phen_contrib_by_GU );
  
  // Tree
  fprintf(setup_file, "RECORD_TREE %s\n", _record_tree ? "true" : "false");
  if (_record_tree)
    fprintf(setup_file, "TREE_STEP %" PRId64 "\n", _tree->get_tree_step());
  
  // Dumps
  fprintf(setup_file, "MAKE_DUMPS %s\n", _make_dumps ? "true" : "false");
  fprintf(setup_file, "DUMP_STEP %" PRId64 "\n", _dump_step);
  
  // Logs
  int8_t logs = _logs->get_logs();
  fprintf(setup_file, "LOGS %" PRId8 "\n", logs);
}

void OutputManager::load( gzFile setup_file, bool verbose, bool to_be_run  )
{
  // Write the backup steps
  gzread( setup_file, &_backup_step,      sizeof(_backup_step) );
  gzread( setup_file, &_big_backup_step,  sizeof(_big_backup_step) );
  
  // Stats
  if (to_be_run)
  {
    if (Time::get_time() > 0)
    {
      _stats = new Stats(_exp_m, Time::get_time());
    }
    else
    {
      _stats = new Stats(_exp_m);
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
    
    _tree = new Tree( _exp_m, tmp_tree_step );
  }
  
  // Dumps
  int8_t make_dumps;
  gzread( setup_file, &make_dumps, sizeof(make_dumps) );
  _make_dumps = make_dumps;
  gzread( setup_file, &_dump_step,  sizeof(_dump_step) );
  if( _make_dumps == true)
  {
    _dump = new Dump(_exp_m);
  }
  
  // Logs
  int8_t logs;
  gzread(setup_file, &logs, sizeof(logs));
  if (to_be_run)
  {
    _logs->load(logs, Time::get_time());
  }
}

void OutputManager::load(FILE* setup_file, bool verbose, bool to_be_run)
{
  int ret;
  // Write the backup steps
  ret = fscanf(setup_file, "BACKUP_STEP %" SCNd64 "\n", &_backup_step);
  ret = fscanf(setup_file, "BIG_BACKUP_STEP %" SCNd64 "\n", &_big_backup_step);
  
  // Stats
  if(to_be_run)
  {
    if (Time::get_time() > 0)
    {
      _stats = new Stats(_exp_m, Time::get_time());
    }
    else
    {
      _stats = new Stats(_exp_m);
    }
  }
  {
    int tmp;
    ret = fscanf( setup_file, "COMPUTE_PHENOTYPIC_CONTRIBUTION_BY_GU %d\n", &tmp);
    _compute_phen_contrib_by_GU = tmp;
  }
  
  char tmp[10];
  
  // Tree
  ret = fscanf(setup_file, "RECORD_TREE %s\n", tmp);
  _record_tree = ! strcmp( tmp, "true" );
  if ( _record_tree )
  {
    int64_t tmp_tree_step;

    ret = fscanf( setup_file, "TREE_STEP %" SCNd64 "\n", &tmp_tree_step );

    
    _tree = new Tree( _exp_m, tmp_tree_step );
  }
  
  // Dumps
  ret = fscanf( setup_file, "MAKE_DUMPS %s\n", tmp );
  _make_dumps = ! strcmp( tmp, "true" );
  ret = fscanf( setup_file, "DUMP_STEP %" SCNd64 "\n", &_dump_step );
  if( _make_dumps == true)
  {
    _dump = new Dump(_exp_m);
  }
  
  // Logs
  int8_t logs;
  ret = fscanf(setup_file, "LOGS %" SCNd8 "\n", &logs);
  _logs->load(logs, Time::get_time());
}

void OutputManager::write_current_generation_outputs( void ) const
{
  // Write stats
  _stats->write_current_generation_statistics();

  // Manage tree
  if (_record_tree &&
      Time::get_time() > 0 &&
      (Time::get_time() % _tree->get_tree_step() == 0)) {
    write_tree();
  }

  // Write backup
  if (Time::get_time() % _backup_step == 0) {
    _stats->flush();
    _exp_m->save();
    
    // Update the LAST_GENER file
    FILE* last_gener_file = fopen(LAST_GENER_FNAME, "w");
    if (last_gener_file != NULL) {
      fprintf(last_gener_file, "%" PRId64 "\n", Time::get_time());
      fclose(last_gener_file);
    }
    else {
      printf( "Error : could not open file " LAST_GENER_FNAME "\n" );
    }
  }

  // Write dumps
  if (_make_dumps) {
    if(Time::get_time() % _dump_step == 0) {
      _dump->write_current_generation_dump();
    }
  }
}

// TODO <david.parsons@inria.fr> we need an input_dir attribute in this class !
int32_t OutputManager::get_last_gener() {
  int32_t num_gener;
  FILE* lg_file = fopen(LAST_GENER_FNAME, "r");
  if (lg_file != NULL) {
    if (fscanf(lg_file, "%" PRId32 "\n", &num_gener) == EOF) {
      Utils::ExitWithMsg("failed to read last generation", __FILE__, __LINE__);
    }
    fclose(lg_file);
  }
  return num_gener;
}

// =================================================================
//                           Protected Methods
// =================================================================
void OutputManager::write_tree( void ) const
{
  // Create the tree directory if it doesn't exist
  int status;
  status = mkdir( TREE_DIR, 0755 );
  if ( (status == -1) && (errno != EEXIST) )
  {
    err( EXIT_FAILURE, "Impossible to create the directory %s", TREE_DIR );
  }
  
  char tree_file_name[50];
  
  sprintf(tree_file_name, "tree/tree_%06" PRId64 ".ae", Time::get_time());

  
  gzFile tree_file = gzopen( tree_file_name, "w" );
  
  // Write phylogenetic data (tree)
  _tree->write_to_tree_file( tree_file );
  
  gzclose( tree_file );
}

// =================================================================
//                          Non inline accessors
// =================================================================
} // namespace aevol
