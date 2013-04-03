//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: http://www.aevol.fr/
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


/*! \class
    \brief
*/
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_params_record.h>




//##############################################################################
//                                                                             #
//                           Class ae_params_record                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_params_record::ae_params_record( void )
{
  // Backups
  _backup_step     = 500;
  _big_backup_step = 10000;
  
  // Tree
  _tree_step   = 100;
  _record_tree = true;
  _tree_mode   = NORMAL;
  
  // Logs
  _logs = new ae_logs();
  
  // Other
  _more_stats  = false;
  _dump_period = 0;
}

// =================================================================
//                             Destructors
// =================================================================
ae_params_record::~ae_params_record( void )
{
  if ( _logs != NULL ) delete _logs;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_params_record::save( gzFile backup_file ) const
{
  // Backups
  gzwrite( backup_file, &_backup_step,     sizeof(_backup_step)     );
  gzwrite( backup_file, &_big_backup_step, sizeof(_big_backup_step) );
  
  // Tree
  gzwrite( backup_file, &_tree_step, sizeof(_tree_step) );
  int8_t tmp_record_tree = _record_tree? 1 : 0;
  gzwrite( backup_file, &tmp_record_tree, sizeof(tmp_record_tree) );
  int8_t tmp_tree_mode = _tree_mode;
  gzwrite( backup_file, &tmp_tree_mode, sizeof(tmp_tree_mode) );
  
  // Logs
  int8_t logs = _logs->get_logs();
  gzwrite( backup_file, &logs,  sizeof(logs) );
  
  // Other
  int8_t tmp_more_stats = _more_stats? 1 : 0;
  gzwrite( backup_file, &tmp_more_stats, sizeof(tmp_more_stats) );
  gzwrite( backup_file, &_dump_period,   sizeof(_dump_period)   );
}

void ae_params_record::load( gzFile backup_file, bool verbose )
{
  // Backups
  gzread( backup_file, &_backup_step,     sizeof(_backup_step)     );
  gzread( backup_file, &_big_backup_step, sizeof(_big_backup_step) );
  
  // Tree
  gzread( backup_file, &_tree_step,       sizeof(_tree_step)       );
  int8_t tmp_record_tree;
  gzread( backup_file, &tmp_record_tree, sizeof(tmp_record_tree) );
  _record_tree = (tmp_record_tree != 0);
  int8_t tmp_tree_mode;
  gzread( backup_file, &tmp_tree_mode, sizeof(tmp_tree_mode) );
  _tree_mode = (ae_tree_mode) tmp_tree_mode;
  
  // Logs
  int8_t logs;
  gzread( backup_file, &logs, sizeof(logs) );
  _logs->load(logs, 0);
  
  // Other
  int8_t tmp_more_stats;
  gzread( backup_file, &tmp_more_stats, sizeof(tmp_more_stats) );
  _more_stats = (tmp_more_stats != 0);
  gzread( backup_file, &_dump_period, sizeof(_dump_period) );
}

void ae_params_record::print_to_file( FILE* file ) const
{
  // Backups
  fprintf( file, "\nRecording parameters ----------------------------------------\n" );
  fprintf( file, "backup_step :                %"PRId32"\n", _backup_step );
  fprintf( file, "big_backup_step              %"PRId32"\n", _big_backup_step);
  
  // Tree
  fprintf( file, "tree_step :                  %"PRId32"\n", _tree_step );
  fprintf( file, "record_tree :                %s\n",  _record_tree? "true" : "false" );
  switch ( _tree_mode )
  {
    case LIGHT :
    {
      fprintf( file, "tree_mode   :                LIGHT\n" );
      break;
    }
    case NORMAL :
    {
      fprintf( file, "tree_mode   :                NORMAL\n" );
      break;
    }
    default :
    {
      fprintf( file, "tree_mode   :                UNKNOWN\n" );
      break;
    }
  }
  
  // Logs
  _logs->print_to_file( file );
  
  // Other
  fprintf( file, "more_stats :                 %s\n",  _more_stats? "true" : "false"  );
  fprintf( file, "dump_period :                %"PRId32"\n",  _dump_period            );
}

// =================================================================
//                           Protected Methods
// =================================================================

// =================================================================
//                          Non inline accessors
// =================================================================
