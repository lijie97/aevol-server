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


/** \class
 *  \brief
 */
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_logs.h>
#include <ae_common.h>




//##############################################################################
//                                                                             #
//                                Class ae_logs                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_logs::ae_logs( void )
{
  _rear_log             = NULL;
  _barrier_log          = NULL;
  _load_from_backup_log = NULL;  
  
  // Open required log files
  if ( ae_common::logs & REAR )
  {
    _rear_log = fopen( "log_rear.out", "w" );
    if ( _rear_log == NULL )
    {
      printf( "Error: Failed to open file \"log_rear.out\"\n" );
      exit( EXIT_FAILURE );
    }
  }
  if ( ae_common::logs & BARRIER )
  {
    _barrier_log = fopen( "log_barrier.out", "w" );
    
    if ( _barrier_log == NULL )
    {
      printf( "Error: Failed to open file \"log_barrier.out\"\n" );
      exit( EXIT_FAILURE );
    }
  }
  if ( ae_common::logs & LOADS )
  {
    _load_from_backup_log = fopen( "log_load_from_backup.out", "w" );
    
    if ( _load_from_backup_log == NULL )
    {
      printf( "Error: Failed to open file \"log_load_from_backup.out\"\n" );
      exit( EXIT_FAILURE );
    }
  }
}

ae_logs::ae_logs( int32_t num_gener )
{
  char* line = new char[500];
  char* ret;
  
  _rear_log             = NULL;
  _barrier_log          = NULL;
  _load_from_backup_log = NULL;
  
  // Prepare required log files
  if ( ae_common::logs & REAR )
  {
    rename ( "log_rear.out", "log_rear.out.old" );
    FILE* old_rear_log = fopen( "log_rear.out.old", "r" );
    if ( old_rear_log == NULL )
    {
      printf( "Error: Failed to open file \"log_rear.out.old\"\n" );
      exit( EXIT_FAILURE );
    }
    
    _rear_log = fopen( "log_rear.out", "w" );
    if ( _rear_log == NULL )
    {
      printf( "Error: Failed to open file \"log_rear.out\"\n" );
      exit( EXIT_FAILURE );
    }
    
    // Copy file headers
    ret = fgets( line, 500, old_rear_log );
    while ( !feof( old_rear_log ) && line[0] == '#' )
    {
      fputs( line, _rear_log );
      ret = fgets( line, 500, old_rear_log );
    }
    // This is the empty line between the header and the values
    fputs( line, _rear_log );
    
    // Copy log entries until num_gener (included)
    ret = fgets( line, 500, old_rear_log );
    while ( (int32_t)atol(line) <= num_gener && !feof(old_rear_log) )
    {
      fputs( line, _rear_log );
      ret = fgets( line, 500, old_rear_log );
    }
    
    fclose( old_rear_log );
  }
  if ( ae_common::logs & BARRIER )
  {
    rename ( "log_barrier.out", "log_barrier.out.old" );
    FILE* old_barrier_log = fopen( "log_barrier.out.old", "r" );
    if ( old_barrier_log == NULL )
    {
      printf( "Error: Failed to open file \"log_barrier.out.old\"\n" );
      exit( EXIT_FAILURE );
    }
    
    _barrier_log = fopen( "log_barrier.out", "w" );
    if ( _barrier_log == NULL )
    {
      printf( "Error: Failed to open file \"log_barrier.out\"\n" );
      exit( EXIT_FAILURE );
    }
    
    // Copy file headers
    ret = fgets( line, 500, old_barrier_log );
    while ( !feof( old_barrier_log ) && line[0] == '#' )
    {
      fputs( line, _barrier_log );
      ret = fgets( line, 500, old_barrier_log );
    }
    // This is the empty line between the header and the values
    fputs( line, _barrier_log );
    
    // Copy log entries until num_gener (included)
    ret = fgets( line, 500, old_barrier_log );
    while ( (int32_t)atol(line) <= num_gener && !feof(old_barrier_log) )
    {
      fputs( line, _barrier_log );
      ret = fgets( line, 500, old_barrier_log );
    }
    
    fclose( old_barrier_log );
  }
  if ( ae_common::logs & LOADS )
  {
    rename ( "log_load_from_backup.out", "log_load_from_backup.out.old" );
    FILE* old_load_from_backup_log = fopen( "log_load_from_backup.out.old", "r" );
    if ( old_load_from_backup_log == NULL )
    {
      printf( "Error: Failed to open file \"log_load_from_backup.out.old\"\n" );
      exit( EXIT_FAILURE );
    }
    
    _load_from_backup_log = fopen( "log_load_from_backup.out", "w" );
    if ( _load_from_backup_log == NULL )
    {
      printf( "Error: Failed to open file \"log_load_from_backup.out\"\n" );
      exit( EXIT_FAILURE );
    }
    
    // Copy file headers
    ret = fgets( line, 500, old_load_from_backup_log );
    while ( !feof( old_load_from_backup_log ) && line[0] == '#' )
    {
      fputs( line, _load_from_backup_log );
      ret = fgets( line, 500, old_load_from_backup_log );
    }
    // This is the empty line between the header and the values
    fputs( line, _load_from_backup_log );
    
    // Copy log entries until num_gener (included)
    ret = fgets( line, 500, old_load_from_backup_log );
    while ( (int32_t)atol(line) <= num_gener && !feof(old_load_from_backup_log) )
    {
      fputs( line, _load_from_backup_log );
      ret = fgets( line, 500, old_load_from_backup_log );
    }
    
    fclose( old_load_from_backup_log );
  }
  
  delete [] line;
}

// =================================================================
//                             Destructors
// =================================================================
ae_logs::~ae_logs( void )
{
  if ( ae_common::logs & REAR )
  {
    fclose( _rear_log );
  }
  if ( ae_common::logs & BARRIER )
  {
    fclose( _barrier_log );
  }
  if ( ae_common::logs & LOADS )
  {
    fclose( _load_from_backup_log );
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_logs::write_headers( void ) const
{
  // ========== REAR LOG ==========
  fprintf( _rear_log, "######################################################################\n" );
  fprintf( _rear_log, "#                 Chromosomal rearrangement log\n" );
  fprintf( _rear_log, "#\n" );
  fprintf( _rear_log, "# Log of every rearrangement that occured during the simulation\n" );
  fprintf( _rear_log, "# (not just one lineage)\n" );
  fprintf( _rear_log, "#\n" );
  fprintf( _rear_log, "# 1. Generation\n" );
  fprintf( _rear_log, "# 2. Index of the individual that has undergone the rearrangement\n" );
  fprintf( _rear_log, "# 3. Type of rearrangement\n" );
  fprintf( _rear_log, "# 4. Length of the rearranged segment\n" );
  fprintf( _rear_log, "# 5. Size of the genome before the rearrangement\n" );
  fprintf( _rear_log, "# 6. Alignment score that was needed for this rearrangement to occur\n" );
  fprintf( _rear_log, "# 7. Second alignment score (translocation only)\n" );
  fprintf( _rear_log, "#\n" );
  fprintf( _rear_log, "######################################################################\n" );
  fprintf( _rear_log, "#\n" );
  fprintf( _rear_log, "# Header for R\n" );
  fprintf( _rear_log, "gener indiv r_type seg_len genome_size score1 score2\n" );
  fprintf( _rear_log, "#\n" );
  
  // ========== BARRIER LOG ==========
  fprintf( _barrier_log, "######################################################################\n" );
  fprintf( _barrier_log, "#                     Genome size limits log\n" );
  fprintf( _barrier_log, "#\n" );
  fprintf( _barrier_log, "# An entry is written whenever a mutation would have produced a\n" );
  fprintf( _barrier_log, "# genome whose size wouldn't lie in [min, max].\n" );
  fprintf( _barrier_log, "# The corresponding mutation is \"cancelled\"\n" );
  fprintf( _barrier_log, "#\n" );
  fprintf( _barrier_log, "# 1. Generation\n" );
  fprintf( _barrier_log, "# 2. Index of the individual\n" );
  fprintf( _barrier_log, "# 3. Type of mutation\n" );
  fprintf( _barrier_log, "# 4. Segment length\n" );
  fprintf( _barrier_log, "# 5. Genome size (before the mutation)\n" );
  fprintf( _barrier_log, "#\n" );
  fprintf( _barrier_log, "######################################################################\n" );
  
  
  // ========== LOADS LOG ==========
  fprintf( _load_from_backup_log, "######################################################################\n" );
  fprintf( _load_from_backup_log, "#                     Load from backup log\n" );
  fprintf( _load_from_backup_log, "#\n" );
  fprintf( _load_from_backup_log, "# An entry is written whenever the simulation is reloaded from a backup.\n" );
  fprintf( _load_from_backup_log, "######################################################################\n" );
}

void ae_logs::flush( void )
{
  if ( ae_common::logs & REAR )
  {
    fflush( _rear_log );
  }
  if ( ae_common::logs & BARRIER )
  {
    fflush( _barrier_log );
  }
  if ( ae_common::logs & LOADS )
  {
    fflush( _load_from_backup_log );
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
