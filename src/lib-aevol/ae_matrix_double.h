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

#ifndef __AE_MATRIX_DOUBLE_H__
#define __AE_MATRIX_DOUBLE_H__

// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <zlib.h>

// =================================================================
//                            Project Files
// =================================================================

// =================================================================
//                          Class declarations
// =================================================================

class ae_matrix_double
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_matrix_double( const int16_t columns, const int16_t rows )
    {
      _row_size     = columns;
      _column_size  = rows;
      _data         = (double*)new double[ _column_size * _row_size ];
      _data_size    = _row_size * _column_size * sizeof(double);
    }

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_matrix_double( void )
    {
      delete [] _data;
    }

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    // assignement and read
    inline void   set_value( const int16_t row, const int16_t column, double value );
    inline double get_value( const int16_t row, const int16_t column ) const;

    // backup functions
    inline void write_to_backup( gzFile* backup_file );
    inline void read_from_backup( gzFile* backup_file );
    inline void print_to_file( FILE* file );

    // =================================================================
    //                           Public Attributes
    // =================================================================

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_matrix_double( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_matrix_double( const ae_matrix_double &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    int16_t _row_size;
    int16_t _column_size;
    double* _data;
    int32_t _data_size;

};

// =================================================================
//                             Constructors
// =================================================================

// =================================================================
//                             Destructors
// =================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
void ae_matrix_double::set_value( const int16_t row, const int16_t column, double value )
{
  _data[ ( row * _column_size ) + column ] = value;
}

double ae_matrix_double::get_value( const int16_t row, const int16_t column ) const
{
  return _data[ ( row * _column_size ) + column ];
}

void ae_matrix_double::write_to_backup( gzFile* backup_file )
{
  gzwrite( backup_file, _data, _data_size );
  //printf( "write %d bytes (matrix_double) :\n", _data_size );
  //print_to_file( stdout );
}

void ae_matrix_double::read_from_backup( gzFile* backup_file )
{
  gzread( backup_file, _data, _data_size );
  //print_to_file( stdout );
}

void ae_matrix_double::print_to_file( FILE* file )
{
  for( int16_t row = 0 ; row < _column_size ; row++ )
  {
    for( int16_t column = 0 ; column < _row_size ; column++ )
    {
      fprintf( file, "\t%e", get_value( row, column ) );
    }
    fprintf( file, "\n");
  }
}

#endif // __AE_MATRIX_DOUBLE_H__
