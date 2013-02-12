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
 
 
 #ifndef __AE_STRING_H__
#define  __AE_STRING_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <zlib.h>
#include <string.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_utils.h>




// =================================================================
//                          Class declarations
// =================================================================







#define BLOCK_SIZE INT32_C(1024)
 
class ae_string : public ae_object
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_string( void );
    ae_string( const ae_string &model );
    ae_string( int32_t length );
    ae_string( const char* seq, int32_t length );
    ae_string( char* seq, int32_t length, bool use_seq );
    ae_string( gzFile backup_file );
    ae_string( char* organism_file_name );
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_string( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
    inline const char*   get_data( void ) const;
    inline       void    set_data( char* data, int32_t length = -1 );
    inline       int32_t get_length( void ) const;
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    void remove( int32_t first, int32_t last );
    void insert( int32_t pos, const char* seq, int32_t seq_length = -1 );
    void replace( int32_t pos, char* seq, int32_t seq_length = -1 );
    
    void save( gzFile backup_file );
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_string( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    /*    ae_string( const ae_string &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      };*/

  
    // =================================================================
    //                           Protected Methods
    // =================================================================
    static inline int32_t nb_blocks( int32_t length );
  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    char*   _data;
    int32_t _length;
    int32_t _nb_blocks;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline const char* ae_string::get_data( void ) const
{
  return _data;
}

inline void ae_string::set_data( char* data, int32_t length /* = -1 */ )
{
  if ( _data != NULL )
  {
    delete [] _data;
    _data = NULL;
  }
  
  _data       = data;
  _length     = ( length != -1 ) ? length : strlen( _data );
  _nb_blocks  = nb_blocks( _length );
}

inline int32_t ae_string::get_length( void ) const
{
  return _length;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
int32_t ae_string::nb_blocks( int32_t length )
{
  return length/BLOCK_SIZE + 1;
}


#endif // __AE_STRING_H__
