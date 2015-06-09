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


 #ifndef __AE_CODON_H__
#define  __AE_CODON_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>





// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_dna.h"
#include "macros.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================






class Codon
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Codon( void );
    Codon( const Codon &model );
    Codon( int8_t value );
    Codon( ae_dna* genome, ae_strand strand, int32_t index );
    Codon( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Codon( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline int8_t get_value( void );

    // =================================================================
    //                            Public Methods
    // =================================================================
    inline bool is_start( void );
    inline bool is_stop( void );
    inline Codon * copy( void );
    void   save( gzFile backup_file );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ Codon( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    //~ Codon( const Codon &model )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };

    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    int8_t _value;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
int8_t Codon::get_value( void )
{
  return _value;
}

bool Codon::is_start( void )
{
  return _value == CODON_START;
}

bool Codon::is_stop( void )
{
  return _value == CODON_STOP;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
Codon *Codon::copy( void )
{
  return new Codon( _value );
}

} // namespace aevol
#endif // __AE_CODON_H__
