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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_mutation.h>




//##############################################################################
//                                                                             #
//                              Class ae_mutation                              #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_mutation::ae_mutation( void )
{
}




ae_mutation::ae_mutation( const ae_mutation &model )
{
  _mut_type = model._mut_type;
  
  switch( _mut_type )
  {
    case SWITCH :
    {
      _pos = new int32_t( *(model._pos) );
      break;
    }
    case S_INS :
    {
      _pos = new int32_t( *(model._pos) );
      _length = model._length;
      _seq = new char[_length + 1];
      strncpy(_seq, model._seq, _length);
      _seq[_length] = '\0';
      break;
    }
    case S_DEL:
    {
      _pos = new int32_t( *(model._pos) );
      _length = model._length;
      break;
    }
    case DUPL:
    {
      _pos = new int32_t[3];
      _pos[0] = model._pos[0];
      _pos[1] = model._pos[1];
      _pos[2] = model._pos[2];
      _align_score = new int16_t( *(model._align_score) );
      break;
    }
    case DEL:
    {
      _pos = new int32_t[2];
      _pos[0] = model._pos[0];
      _pos[1] = model._pos[1];
      _align_score = new int16_t( *(model._align_score) );
      break;
    }
    case TRANS:
    {
      _pos = new int32_t[4];
      _pos[0] = model._pos[0];
      _pos[1] = model._pos[1];
      _pos[2] = model._pos[2];
      _pos[3] = model._pos[3];
      _invert = model._invert;
      _align_score = new int16_t[2];
      _align_score[0] = model._align_score[0];
      _align_score[1] = model._align_score[1];
      break;
    }
    case INV:
    {
      _pos = new int32_t[2];
      _pos[0] = model._pos[0];
      _pos[1] = model._pos[1];
      _align_score = new int16_t( *(model._align_score) );
      break;
    }
    default :
    {
      fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", model._mut_type, __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      break;
    }
  }
}




ae_mutation::ae_mutation( gzFile* backup_file )
{
  int8_t tmp_mut_type;
  gzread( backup_file, &tmp_mut_type,  sizeof(tmp_mut_type) );
  _mut_type = (ae_mutation_type) tmp_mut_type;
  //~ printf("mut type %d\n", _mut_type);

  switch ( _mut_type )
  {
    case SWITCH :
    {
      _pos = new int32_t;
      gzread( backup_file, _pos,  sizeof(*_pos) );
      break;
    }
    case S_INS :
    {
      _pos = new int32_t;
      gzread( backup_file, _pos,      sizeof(*_pos) );
      gzread( backup_file, &_length,  sizeof(_length) );
      
      _seq = new char[_length + 1];
      gzread( backup_file, _seq,  _length * sizeof(_seq[0]) );
      _seq[_length] = '\0';
      break;
    }
    case S_DEL :
    {
      _pos = new int32_t;
      gzread( backup_file, _pos,      sizeof(*_pos) );
      gzread( backup_file, &_length,  sizeof(_length) );
      break;
    }
    case DUPL :
    {
      _pos = new int32_t[3];
      gzread( backup_file, _pos,  3 * sizeof(_pos[0]) );
      
      _align_score = new int16_t;
      gzread( backup_file, _align_score, sizeof(*_align_score) );
      
      break;
    }
    case DEL :
    {
      _pos = new int32_t[2];
      gzread( backup_file, _pos,  2 * sizeof(_pos[0]) );
      
      _align_score = new int16_t;
      gzread( backup_file, _align_score, sizeof(*_align_score) );
      
      break;
    }
    case TRANS :
    {
      _pos = new int32_t[4];
      gzread( backup_file, _pos,  4 * sizeof(_pos[0]) );
      
      int8_t tmp_invert;
      gzread( backup_file, &tmp_invert,  sizeof(tmp_invert) );
      _invert = (tmp_invert != 0);
      
      _align_score = new int16_t[2];
      gzread( backup_file, _align_score, 2 * sizeof(_align_score[0]) );
      
      break;
    }
    case INV :
    {
      _pos = new int32_t[2];
      gzread( backup_file, _pos,  2 * sizeof(_pos[0]) );
      
      _align_score = new int16_t;
      gzread( backup_file, _align_score, sizeof(*_align_score) );
      
      break;
    }
    default :
    {
      fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", _mut_type, __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      break;
    }
  }


}





// =================================================================
//                             Destructors
// =================================================================
ae_mutation::~ae_mutation( void )
{
  switch ( _mut_type )
  {
    case SWITCH :
      delete _pos;
      break;
    case S_INS :
      delete _pos;
      delete [] _seq;
      break;
    case S_DEL :
      delete _pos;
      break;
    case DUPL :
      delete [] _pos;
      delete _align_score;
      break;
    case DEL :
      delete [] _pos;
      delete _align_score;
      break;
    case TRANS :
      delete [] _pos;
      delete [] _align_score;
      break;
    case INV :
      delete [] _pos;
      delete _align_score;
      break;
    case INSERT :
      delete _pos;
      delete [] _seq;
      break;
    default :
      fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", _mut_type, __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      break;
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_mutation::report_point_mutation( int32_t pos )
{
  //~ printf( "report point mutation at %ld\n", pos );
  
  _mut_type = SWITCH;
  _pos      = new int32_t( pos );
  _length   = 1;
}

void ae_mutation::report_small_insertion( int32_t pos, int32_t length, const char* seq )
{
  //~ printf( "report small insertion of %s (%ld) at %ld\n", seq, length, pos );
  
  _mut_type = S_INS;
  _pos = new int32_t( pos );
  _length = length;
  _seq = new char[length + 1];
  strncpy(_seq, seq, length);
  _seq[length] = '\0';
}

void ae_mutation::report_small_deletion( int32_t pos, int32_t length )
{
  //~ printf( "report small deletion of %ld at %ld\n", length, pos );
  
  _mut_type = S_DEL;
  _pos = new int32_t( pos );
  _length = length;
}

void ae_mutation::report_duplication( int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t length, int16_t align_score )
{
  //~ printf( "report duplication of [%ld, %ld] to %ld\n", pos_1, pos_2, pos_3 );
  
  _mut_type = DUPL;
  
  _pos = new int32_t[3];
  _pos[0] = pos_1;
  _pos[1] = pos_2;
  _pos[2] = pos_3;
  
  _length = length;
  
  _align_score = new int16_t( align_score );
}

void ae_mutation::report_deletion( int32_t pos_1, int32_t pos_2, int32_t length, int16_t align_score )
{
  //~ printf( "report deletion of [%ld, %ld]\n", pos_1, pos_2 );
  
  _mut_type = DEL;
  
  _pos = new int32_t[2];
  _pos[0] = pos_1;
  _pos[1] = pos_2;
  
  _length = length;
  
  _align_score = new int16_t( align_score );
}

void ae_mutation::report_translocation( int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, int32_t length, bool invert, int16_t align_score_1, int16_t align_score_2 )
{
  //~ printf( "report translocation of [%ld, %ld] to %ld through %ld (%s)\n", pos_1, pos_2, pos_3, pos_4, invert?"invert":"plain" );
  
  _mut_type = TRANS;
  
  _pos = new int32_t[4];
  _pos[0] = pos_1;
  _pos[1] = pos_2;
  _pos[2] = pos_3;
  _pos[3] = pos_4;
  _invert = invert;
  
  _length = length;
  
  _align_score = new int16_t[2];
  _align_score[0] = align_score_1;
  _align_score[1] = align_score_2;
}

void ae_mutation::report_inversion( int32_t pos_1, int32_t pos_2, int32_t length, int16_t align_score )
{
  //~ printf( "report inversion of [%ld, %ld]\n", pos_1, pos_2 );
  
  _mut_type = INV;
  
  _pos = new int32_t[2];
  _pos[0] = pos_1;
  _pos[1] = pos_2;
  
  _length = length;
  
  _align_score = new int16_t( align_score );
}

void ae_mutation::report_insertion( int32_t pos, int32_t length, const char* seq )
{
  //~ printf( "report insertion of %s (%ld) at %ld\n", seq, length, pos );
  
  _mut_type = INSERT;
  _pos = new int32_t( pos );
  _length = length;
  _seq = new char[length + 1];
  strncpy(_seq, seq, length);
  _seq[length] = '\0';
}

void ae_mutation::save( gzFile* backup_file ) // Usually <backup_file> is the tree file (?)
{
  int8_t tmp_mut_type = _mut_type;
  gzwrite( backup_file, &tmp_mut_type,  sizeof(tmp_mut_type) );

  switch ( _mut_type )
  {
    case SWITCH :
    {
      gzwrite( backup_file, _pos,  sizeof(*_pos) );
      break;
    }
    case S_INS :
    {
      gzwrite( backup_file, _pos,     sizeof(*_pos) );
      gzwrite( backup_file, &_length, sizeof(_length) );
      gzwrite( backup_file, _seq,     _length * sizeof(_seq[0]) );
      break;
    }
    case S_DEL :
    {
      gzwrite( backup_file, _pos,     sizeof(*_pos) );
      gzwrite( backup_file, &_length, sizeof(_length) );
      break;
    }
    case DUPL :
    {
      gzwrite( backup_file, _pos,         3 * sizeof(_pos[0]) );
      gzwrite( backup_file, _align_score, sizeof(*_align_score) );
      break;
    }
    case DEL :
    {
      gzwrite( backup_file, _pos,         2 * sizeof(_pos[0]) );
      gzwrite( backup_file, _align_score, sizeof(*_align_score) );
      break;
    }
    case TRANS :
    {
      gzwrite( backup_file, _pos,         4 * sizeof(_pos[0]) );
      int8_t tmp_invert = _invert? 1 : 0;
      gzwrite( backup_file, &tmp_invert,  sizeof(tmp_invert) );
      gzwrite( backup_file, _align_score, 2 * sizeof(_align_score[0]) );
      break;
    }
    case INV :
    {
      gzwrite( backup_file, _pos,         2 * sizeof(_pos[0]) );
      gzwrite( backup_file, _align_score, sizeof(*_align_score) );
      break;
    }
    default :
    {
      fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", _mut_type, __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      break;
    }
  }
}




void ae_mutation::get_generic_description_string( char * str )
{
  switch ( _mut_type )
  {
    case SWITCH :
    {
      sprintf( str, "%"PRId8" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId8" %"PRId16" %"PRId16" ", \
              (int8_t) _mut_type, _pos[0], -1, -1, -1, (int8_t) -1, (int16_t) -1, (int16_t) -1 );
      break;
    }
    case S_INS :
    {
      sprintf( str, "%"PRId8" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId8" %"PRId16" %"PRId16" ", \
              (int8_t) _mut_type, _pos[0], -1, -1, -1, (int8_t) -1, (int16_t) -1, (int16_t) -1 );
      break;
    }
    case S_DEL :
    {
      sprintf( str, "%"PRId8" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId8" %"PRId16" %"PRId16" ", \
              (int8_t) _mut_type, _pos[0], -1, -1, -1, (int8_t) -1, (int16_t) -1, (int16_t) -1 );
      break;
    }
    case DUPL :
    {
      sprintf( str, "%"PRId8" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId8" %"PRId16" %"PRId16" ", \
              (int8_t) _mut_type, _pos[0], _pos[1], _pos[2], -1, (int8_t) -1, *_align_score, (int16_t) -1 );
      break;
    }
    case DEL :
    {
      sprintf( str, "%"PRId8" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId8" %"PRId16" %"PRId16" ", \
              (int8_t) _mut_type, _pos[0], _pos[1], -1, -1, (int8_t) -1, *_align_score, (int16_t) -1 );
      break;
    }
    case TRANS :
    {
      int8_t tmp_invert = _invert? 1 : 0;
      sprintf( str, "%"PRId8" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId8" %"PRId16" %"PRId16" ", \
              (int8_t) _mut_type, _pos[0], _pos[1], _pos[2], _pos[3], tmp_invert, _align_score[0], _align_score[1] );
      break;
    }
    case INV :
    {
      sprintf( str, "%"PRId8" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId8" %"PRId16" %"PRId16" ", \
              (int8_t) _mut_type, _pos[0], _pos[1], -1, -1, (int8_t) -1, *_align_score, (int16_t) -1 );
      break;
    }
    default :
    {
      fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", _mut_type, __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      break;
    }
  }

}



/* DEPRECATED, use get_length instead */
int32_t ae_mutation::segment_length( int32_t gen_unit_len )
{
  int32_t res;
  switch ( _mut_type )
  {
    case SWITCH :
    {
      res = 1;
      break;
    }
    case S_INS :
    {
      res = _length;
      break;
    }
    case S_DEL :
    {
      res = _length;
      break;
    }
    case DUPL :
    {
      // pos 0 = begin_seg;
      // pos 1 = end_seg;
      // pos 2 = insertion_point;

      if (_pos[0] <= _pos[1]) res = _pos[1] - _pos[0];
      else res = gen_unit_len - _pos[0] + _pos[1];

      break;
    }
    case DEL :
    {
      // pos 0 = begin_seg;
      // pos 1 = end_seg;

      if (_pos[0] <= _pos[1]) res = _pos[1] - _pos[0];
      else res = gen_unit_len - _pos[0] + _pos[1];
      
      break;
    }
    case TRANS :
    {
      // pos 0 = begin_seg;
      // pos 1 = end_seg;
      // pos 2 = cutting_point_in_translocated_segment (1st reinserted bp)
      // pos 3 = reinsertion_point_in_chromosome;

      if (_pos[0] <= _pos[1]) res = _pos[1] - _pos[0];
      else res = gen_unit_len - _pos[0] + _pos[1];

      break;
    }
    case INV :
    {
      // pos 0 = begin_seg;
      // pos 1 = end_seg;

      if (_pos[0] <= _pos[1]) res = _pos[1] - _pos[0];
      else res = gen_unit_len - _pos[0] + _pos[1];

      break;
    }
    default :
    {
      fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", _mut_type, __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      break;
    }
  }

  return res;
}





void ae_mutation::get_infos_point_mutation( int32_t* pos )
{
  assert( _mut_type == SWITCH );
  *pos = _pos[0];
}



void ae_mutation::get_infos_small_insertion( int32_t* pos, int32_t* length )
{
  assert( _mut_type == S_INS );
  *pos    = *_pos;
  *length = _length;
}


// seq must be a char array, large enough to contain _length+1 characters
void ae_mutation::get_sequence_small_insertion( char* seq )
{
  assert( _mut_type == S_INS );
  strncpy( seq, _seq, _length);
  seq[_length] = '\0';
}

void ae_mutation::get_infos_small_deletion( int32_t* pos, int32_t* length )
{
  assert( _mut_type == S_DEL );
  *pos    = _pos[0];
  *length = _length;
}


void ae_mutation::get_infos_duplication( int32_t* pos1, int32_t* pos2, int32_t* pos3, int16_t* align_score )
{
  assert( _mut_type == DUPL );
  *pos1    = _pos[0];
  *pos2    = _pos[1];
  *pos3    = _pos[2];

  if ( align_score != NULL )
  {
    *align_score = *_align_score;
  }
}


void ae_mutation::get_infos_deletion( int32_t* pos1, int32_t* pos2, int16_t* align_score )
{
  assert( _mut_type == DEL );
  *pos1    = _pos[0];
  *pos2    = _pos[1];

  if ( align_score != NULL )
  {
    *align_score = *_align_score;
  }
}


void ae_mutation::get_infos_translocation( int32_t* pos1, int32_t* pos2, int32_t* pos3, int32_t* pos4, bool* invert,
                                           int16_t* align_score_1, int16_t* align_score_2 )
{
  assert( _mut_type == TRANS );
  *pos1    = _pos[0];
  *pos2    = _pos[1];
  *pos3    = _pos[2];
  *pos4    = _pos[3];
  *invert  = _invert;

  if ( align_score_1 != NULL && align_score_2 != NULL )
  {
    *align_score_1 = _align_score[0];
    *align_score_2 = _align_score[1];
  }
}


void ae_mutation::get_infos_inversion( int32_t* pos1, int32_t* pos2, int16_t* align_score )
{
  assert( _mut_type == INV );
  *pos1    = _pos[0];
  *pos2    = _pos[1];

  if ( align_score != NULL )
  {
    *align_score = *_align_score;
  }
}

void ae_mutation::get_infos_insertion( int32_t* pos, int32_t* length )
{
  assert( _mut_type == INSERT );
  *pos    = *_pos;
  *length = _length;
}



// =================================================================
//                           Protected Methods
// =================================================================
