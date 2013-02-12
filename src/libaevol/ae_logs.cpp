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
#include <ae_exp_setup.h>




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
  _logs = 0;
  
  _transfer_log         = NULL;
  _rear_log             = NULL;
  _barrier_log          = NULL;
  _load_from_backup_log = NULL;
}

// =================================================================
//                             Destructors
// =================================================================
ae_logs::~ae_logs( void )
{
  if ( _logs & LOG_TRANSFER )
  {
    fclose( _transfer_log );
  }
  if ( _logs & LOG_REAR )
  {
    fclose( _rear_log );
  }
  if ( _logs & LOG_BARRIER )
  {
    fclose( _barrier_log );
  }
  if ( _logs & LOG_LOADS )
  {
    fclose( _load_from_backup_log );
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_logs::save( gzFile backup_file ) const
{
  gzwrite( backup_file, &_logs, sizeof(_logs) );
}

void ae_logs::load( gzFile backup_file )
{
  char* line = new char[500];
  char* ret;
  
  // TODO Temporary save, should be gotten rid of to use the value in ae_common::sim
  int32_t num_gener;
  gzread( backup_file, &num_gener, sizeof(num_gener) );
  
  gzread( backup_file, &_logs, sizeof(_logs) );
  
  // Prepare required log files
  if ( _logs & LOG_TRANSFER )
  {
    rename ( "log_transfer.out", "log_transfer.out.old" );
    FILE* old_transfer_log = fopen( "log_transfer.out.old", "r" );
    if ( old_transfer_log == NULL )
    {
      printf( "Error: Failed to open file \"log_transfer.out.old\"\n" );
      exit( EXIT_FAILURE );
    }
    
    _transfer_log = fopen( "log_transfer.out", "w" );
    if ( _transfer_log == NULL )
    {
      printf( "Error: Failed to open file \"log_transfer.out\"\n" );
      exit( EXIT_FAILURE );
    }
    
    // Copy file headers
    ret = fgets( line, 500, old_transfer_log );
    while ( !feof( old_transfer_log ) && line[0] == '#' )
    {
      fputs( line, _transfer_log );
      ret = fgets( line, 500, old_transfer_log );
    }
    // This is the empty line between the header and the values
    fputs( line, _transfer_log );
    
    // Copy log entries until num_gener (included)
    ret = fgets( line, 500, old_transfer_log );
    while ( (int32_t)atol(line) <= num_gener && !feof(old_transfer_log) )
    {
      fputs( line, _transfer_log );
      ret = fgets( line, 500, old_transfer_log );
    }
    
    fclose( old_transfer_log );
  }
  if ( _logs & LOG_REAR )
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
  if ( _logs & LOG_BARRIER )
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
  if ( _logs & LOG_LOADS )
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
    //ret = fgets( line, 500, old_load_from_backup_log );
    //while ( !feof( old_load_from_backup_log ) && line[0] == '#' )
    //{
    //  fputs( line, _load_from_backup_log );
    //  ret = fgets( line, 500, old_load_from_backup_log );
    //}
    // This is the empty line between the header and the values
    //fputs( line, _load_from_backup_log );
    
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

void ae_logs::print_to_file( FILE* file ) const
{
  fprintf( file, "logs        :                %"PRId8"\n", _logs );
}

void ae_logs::set_logs( int8_t logs )
{
  _logs = logs;
    
  // Open required log files
  if ( _logs & LOG_TRANSFER )
  {
    _transfer_log = fopen( "log_transfer.out", "w" );
    if ( _transfer_log == NULL )
    {
      printf( "Error: Failed to open file \"log_transfer.out\"\n" );
      exit( EXIT_FAILURE );
    }
  }
  if ( _logs & LOG_REAR )
  {
    _rear_log = fopen( "log_rear.out", "w" );
    if ( _rear_log == NULL )
    {
      printf( "Error: Failed to open file \"log_rear.out\"\n" );
      exit( EXIT_FAILURE );
    }
  }
  if ( _logs & LOG_BARRIER )
  {
    _barrier_log = fopen( "log_barrier.out", "w" );
    
    if ( _barrier_log == NULL )
    {
      printf( "Error: Failed to open file \"log_barrier.out\"\n" );
      exit( EXIT_FAILURE );
    }
  }
  if ( _logs & LOG_LOADS )
  {
    _load_from_backup_log = fopen( "log_load_from_backup.out", "w" );
    
    if ( _load_from_backup_log == NULL )
    {
      printf( "Error: Failed to open file \"log_load_from_backup.out\"\n" );
      exit( EXIT_FAILURE );
    }
  }
  
  this->write_headers();
}

void ae_logs::flush( void )
{
  if ( _logs & LOG_TRANSFER )
  {
    fflush( _transfer_log );
  }
  if ( _logs & LOG_REAR )
  {
    fflush( _rear_log );
  }
  if ( _logs & LOG_BARRIER )
  {
    fflush( _barrier_log );
  }
  if ( _logs & LOG_LOADS )
  {
    fflush( _load_from_backup_log );
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
void ae_logs::write_headers( void ) const
{
  // ========== TRANSFER LOG ==========
  if ( _logs & LOG_TRANSFER )
  {
    fprintf( _transfer_log, "######################################################################\n" );
    fprintf( _transfer_log, "#                 Horizontal transfer log\n" );
    fprintf( _transfer_log, "#\n" );
    fprintf( _transfer_log, "# Log of every horizontal transfer that occured during the simulation\n" );
    fprintf( _transfer_log, "#\n" );
    fprintf( _transfer_log, "# 1.  Generation\n" );
    fprintf( _transfer_log, "# 2.  Index of the recepient\n" );
    fprintf( _transfer_log, "# 3.  Index of the donor (generation n-1)\n" );
    fprintf( _transfer_log, "# 4.  Type of transfer\n" );
    fprintf( _transfer_log, "# 5.  Length of the transferred segment\n" );
    fprintf( _transfer_log, "# 6.  Length of the replaced segment (if any)\n" );
    fprintf( _transfer_log, "# 7.  Size of the genome before the transfer\n" );
    fprintf( _transfer_log, "# 8.  Size of the genome afetr the transfer\n" );
    fprintf( _transfer_log, "# 9.  Alignment 1 point 1\n" );
    fprintf( _transfer_log, "# 10. Alignment 1 point 2\n" );
    fprintf( _transfer_log, "# 11. Alignment 1 score\n" );
    fprintf( _transfer_log, "# 12. Alignment 2 point 1\n" );
    fprintf( _transfer_log, "# 13. Alignment 2 point 2\n" );
    fprintf( _transfer_log, "# 14. Alignment 2 score\n" );
    fprintf( _transfer_log, "#\n" );
    fprintf( _transfer_log, "######################################################################\n" );
    fprintf( _transfer_log, "#\n" );
    fprintf( _transfer_log, "# Header for R\n" );
    fprintf( _transfer_log, "gener recepient donor t_type seg_len replaced_len size_before size_after align1_pt1 align1_pt2 score1 align2_pt1 align2_pt2 score2\n" );
    fprintf( _transfer_log, "#\n" );
  }
  
  // ========== REAR LOG ==========
  if ( _logs & LOG_REAR )
  {
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
  }
  
  // ========== BARRIER LOG ==========
  if ( _logs & LOG_BARRIER )
  {
    fprintf( _barrier_log, "######################################################################\n" );
    fprintf( _barrier_log, "#                     Genome size limits log\n" );
    fprintf( _barrier_log, "#\n" );
    fprintf( _barrier_log, "# An entry is written whenever a mutation would have produced a\n" );
    fprintf( _barrier_log, "# genome whose size wouldn't lie in [min, max].\n" );
    fprintf( _barrier_log, "# The corresponding mutation is \"cancelled\"\n" );
    fprintf( _barrier_log, "#\n" );
    fprintf( _barrier_log, "# 1. Generation\n" );
    fprintf( _barrier_log, "# 2. Index of the individual\n" );
    fprintf( _barrier_log, "# 3. Type of event\n" );
    fprintf( _barrier_log, "# 4. Segment length\n" );
    fprintf( _barrier_log, "# 5. Replaced segment length\n" );
    fprintf( _barrier_log, "# 6. Genome size (before the event)\n" );
    fprintf( _barrier_log, "#\n" );
    fprintf( _barrier_log, "######################################################################\n" );
  
  }
  
  // ========== LOADS LOG ==========     
  //No header because this file is used by post-treatment like ancstat, it's not meant to be human readable but has to be program readable.
  /*if ( _logs & LOG_LOADS )
  {
    fprintf( _load_from_backup_log, "######################################################################\n" );
    fprintf( _load_from_backup_log, "#                     Load from backup log\n" );
    fprintf( _load_from_backup_log, "#\n" );
    fprintf( _load_from_backup_log, "# An entry is written whenever the simulation is reloaded from a backup.\n" );
    fprintf( _load_from_backup_log, "######################################################################\n" );
  }*/
}
