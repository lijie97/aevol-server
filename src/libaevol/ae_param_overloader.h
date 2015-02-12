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


#ifndef __AE_PARAM_OVERLOADER_H__
#define __AE_PARAM_OVERLOADER_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>


// =================================================================
//                            Project Files
// =================================================================
#include <stdlib.h>
#include <string.h>
#include "ae_param_loader.h"

#define NB_OVERLOADABLE INT32_C(256)
#define LEN_TOKEN INT32_C(256)

namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================

class ae_param_overloader
{
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_param_overloader( void );
    
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_param_overloader( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
    inline int16_t get_nb_overloaded( void ) const;
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    inline void store_overload( char* token );
    void overload_params( void ) const;
    void write_log( FILE* log_file ) const;
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
    
    
    
    
 protected :
    
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_param_loader( void )
    //~ {
    //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
    //~ exit( EXIT_FAILURE );
    //~ };
    ae_param_overloader( const ae_param_overloader &model )
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
    int16_t _nb_overloaded;
    char*  _tokens_overloaded[NB_OVERLOADABLE];
};



// =====================================================================
//                          Accessors definitions
// =====================================================================
inline int16_t ae_param_overloader::get_nb_overloaded( void ) const
{
  return _nb_overloaded;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_param_overloader::store_overload( char* token )
{
  _tokens_overloaded[_nb_overloaded] = new char[strlen(token)+1];
  strcpy( _tokens_overloaded[_nb_overloaded++], token );
}

} // namespace aevol

#endif // __AE_PARAM_OVERLOADER_H__
