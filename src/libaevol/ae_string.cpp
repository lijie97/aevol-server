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
#include <assert.h>
#include <string.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_string.h>
#include <ae_exp_setup.h>




//##############################################################################
//                                                                             #
//                               Class ae_string                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_string::ae_string( void )
{
  _nb_blocks = 1;
  _length = 0;
  _data = new char[_nb_blocks * BLOCK_SIZE * sizeof(char)];
  _data[_length] = '\0';
}

ae_string::ae_string( const ae_string &model )
{
  _nb_blocks  = model._nb_blocks;
  _length     = model._length;
  _data = new char[_nb_blocks * BLOCK_SIZE * sizeof(char)];

  memcpy( _data, model._data, (_length+1) * sizeof(char) );
}

/*!
  Creates a new ae_string with enough space to store a string of the given length.
*/
ae_string::ae_string( int32_t length )
{
  _nb_blocks = nb_blocks( length );
  _length = length;
  _data = new char[_nb_blocks * BLOCK_SIZE];

  // Generate a random genome
  /*for ( int32_t i = 0 ; i < _length ; i++ )
  {
    _data[i] = '0' + ae_common::sim->alea->random( NB_BASE );
  }
  _data[_length] = '\0';*/
}

/**
 * Creates a new ae_string with sequence <seq> (having length <length>)
 */
ae_string::ae_string( const char* seq, int32_t length )
{
  _length = length;
  _nb_blocks = nb_blocks( length );
  _data = new char[_nb_blocks * BLOCK_SIZE];
  memcpy( _data, seq, (length+1) * sizeof(char) );
}

/**
 * Creates a new ae_string with sequence <seq> (having length <length>).
 * WARNING : <seq> is used directly which means the caller must not delete it.
 */
ae_string::ae_string( char* seq, int32_t length, bool use_seq )
{
  assert( use_seq );
  
  _length = length;
  _nb_blocks = nb_blocks( length );
  _data = seq;
}

ae_string::ae_string( gzFile* backup_file )
{
  gzread( backup_file, &_nb_blocks,  sizeof(_nb_blocks) );
  //~ printf( "read %d bytes (_nb_blocks : %ld)\n", sizeof(_nb_blocks), _nb_blocks );
  _data = new char[_nb_blocks * BLOCK_SIZE];
  gzread( backup_file, &_length,     sizeof(_length) );
  //~ printf( "read %d bytes (_length : %ld)\n", sizeof(_length), _length );
  gzread( backup_file, _data,        (_length + 1) * sizeof(*_data) );
  //~ printf( "read %ld bytes (_data : %s)\n", (_length + 1) * sizeof(*_data), _data );
}

ae_string::ae_string( char* organism_file_name )
{
  FILE* org_file = fopen( organism_file_name, "r" );
  int length; 
  
  if (org_file==NULL) printf ("Error opening organism file\n");
  else
  {
    fseek (org_file , 0 , SEEK_END);
    length = ftell (org_file);
    rewind (org_file);
    
    _nb_blocks = nb_blocks( length );
    _length = length;
    _data = new char[_nb_blocks * BLOCK_SIZE];
    for ( int32_t i = 0 ; i < _length -1 ; i++ )
    {
      _data[i] = fgetc (org_file);

    }
    
    _data[_length] = '\0';
  }
  
  fclose ( org_file ); 
  
}
// =================================================================
//                             Destructors
// =================================================================
ae_string::~ae_string( void )
{
  delete [] _data;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_string::remove( int32_t pos_1, int32_t pos_2 )
// Remove the sequence between positions 'first' and 'last'
{
  assert( pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= _length );

  // Compute size of new genome
  int32_t new_length    = _length - (pos_2 - pos_1);
  int32_t new_nb_blocks = nb_blocks( new_length );
  char*   new_genome    = new char[new_nb_blocks * BLOCK_SIZE];

  // Copy the remaining of the genome in tmp (preceeding and following parts)
  memcpy( new_genome, _data, pos_1 * sizeof(char) );
  memcpy( &new_genome[pos_1], &_data[pos_2], (new_length - pos_1) * sizeof(char) );
  new_genome[new_length] = '\0';

  // Replace previous genome with the new one
  //~ printf( "genome before - after delete_block( %ld, %ld ) : \n  %s\n  %s\n\n", pos_1, pos_2, _data, new_genome );
  //~ getchar();
  delete [] _data;
  _data = new_genome;

  // Update length data
  _length     = new_length ;
  _nb_blocks  = new_nb_blocks;
}

void ae_string::insert( int32_t pos, const char* seq, int32_t seq_length )
// Insert sequence 'seq' at position 'pos'
{
  assert( pos >= 0 && pos < _length );

  // If the sequence's length was not provided, compute it
  if ( seq_length == -1 )
  {
    seq_length = strlen( seq );
  }

  // Compute size of new genome
  int32_t new_length    = _length + seq_length;
  int32_t new_nb_blocks = nb_blocks( new_length );
  char*   new_genome    = new char[new_nb_blocks * BLOCK_SIZE];

  // Build new genome from previous genome and sequence to insert
  memcpy( new_genome,                   _data,        pos             * sizeof(char) );
  memcpy( &new_genome[pos],             seq,          seq_length      * sizeof(char) );
  memcpy( &new_genome[pos+seq_length],  &_data[pos],  (_length - pos) * sizeof(char) );
  new_genome[new_length] = '\0';

  // Replace the previous genome with the new one
  delete [] _data;
  _data = new_genome;

  // Update length-related data
  _length     = new_length;
  _nb_blocks  = new_nb_blocks;
}

void ae_string::replace( int32_t pos, char* seq, int32_t seq_length )
// Invert the sequence between positions 'first' and 'last'
{
  // Check pos value
  assert( pos >= 0 && pos < _length );

  // If the sequence's length was not provided, compute it
  if ( seq_length == -1 )
  {
    seq_length = strlen( seq );
  }
  
  // Check that the sequence is contiguous
  assert( pos + seq_length <= _length );

  // Perform the replacement
  memcpy( &_data[pos], seq, seq_length * sizeof(char) );
}

void ae_string::write_to_backup( gzFile* backup_file )
{
  gzwrite( backup_file, &_nb_blocks,  sizeof(_nb_blocks) );
  //~ printf( "write %d bytes (_nb_blocks : %ld)\n", sizeof(_nb_blocks), _nb_blocks );
  gzwrite( backup_file, &_length,     sizeof(_length) );
  //~ printf( "write %d bytes (_length : %ld)\n", sizeof(_length), _length );
  gzwrite( backup_file, &_data[0],    (_length + 1) * sizeof(_data[0]) );
  //~ printf( "write %ld bytes (_data : %s)\n", (_length + 1) * sizeof(*_data), _data );
}


// =================================================================
//                           Protected Methods
// =================================================================
