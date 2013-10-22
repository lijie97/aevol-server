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
// ****************************************************************************


#ifndef __AE_ARRAY_SHORT_H__
#define __AE_ARRAY_SHORT_H__

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




class ae_array_short
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_array_short( const int32_t size )
    {
      _size       = size;
      _data       = new int16_t[ size ];
    }

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_array_short( void )
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
    inline void set_value( const int32_t index, int16_t value );
    inline int16_t get_value( const int32_t index );

    // backup functions
    inline void save( gzFile backup_file );
    inline void load( gzFile backup_file );
    inline void print_to_file( FILE* file );

    // search
    int32_t search( int16_t value )
    // returns -1 when not found
    {
      for( int32_t index = 0 ; index < _size ; index++ )
      {
        if( get_value( index ) == value )
        {
          return index;
        }
      }
      return -1;
    };

    // sort
    void sort( void )
    {
      // we search the maximum in [0, last_sorted [ and then swap it with the last not sorted = last_sorted-1
      for( int32_t last_sorted = _size ; last_sorted > 0 ; last_sorted-- )
      {
        int32_t index_of_max = 0;
        for( int32_t index_of_current = 1; index_of_current < last_sorted; index_of_current++ )
        {
          index_of_max = ( get_value( index_of_max ) < get_value( index_of_current ) ) ? index_of_current : index_of_max;
        }
        int16_t temp = get_value( last_sorted - 1 );
        set_value( last_sorted - 1, get_value( index_of_max ) );
        set_value( index_of_max, temp );
      }
    }

    // =================================================================
    //                           Public Attributes
    // =================================================================

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_array_short( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_array_short( const ae_array_short &model )
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
    int32_t  _size;
    int16_t* _data;
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
void ae_array_short::set_value( const int32_t index, int16_t value )
{
  _data[ index ] = value;
}

int16_t ae_array_short::get_value( const int32_t index )
{
  return _data[ index ];
}

void ae_array_short::save( gzFile backup_file )
{
  gzwrite( backup_file, &_data[0], _size * sizeof(_data[0]) );
  //print_to_file( stdout );
}

void ae_array_short::load( gzFile backup_file )
{
  gzread( backup_file, _data, _size * sizeof(_data[0]) );
  //print_to_file( stdout );
}

void ae_array_short::print_to_file( FILE* file )
{
  for( int32_t index = 0 ; index < _size ; index++ )
  {
    fprintf( file, "\t%d", get_value( index ) );
  }
  fprintf( file, "\n" );
}

#endif // __AE_ARRAY_SHORT_H__
