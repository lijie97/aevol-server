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
//*****************************************************************************




// =================================================================
//                              Libraries
// =================================================================
#include <math.h>



// =================================================================
//                            Project Files
// =================================================================
#include "Codon.h"
#include "Utils.h"

namespace aevol {



//##############################################################################
//                                                                             #
//                                Class Codon                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Codon::Codon( void )
{
  _value = -1;
}

Codon::Codon( const Codon &model )
{
  _value = model._value;
}

Codon::Codon( int8_t value )
{
  _value = value;
}

Codon::Codon( Dna * dna, Strand strand, int32_t index )
{
  const char* gen = dna->get_data();
  int32_t     len = dna->get_length();

  _value = 0;

  if ( strand == LEADING )
  {
    for ( int8_t i = 0 ; i < CODON_SIZE ; i++ )
    {
      if ( gen[Utils::mod((index+i),len)] == '1' )
      {
        //_value += pow( 2, CODON_SIZE - i - 1 );
        _value += 1 << ( CODON_SIZE - i - 1 );
      }
    }
  }
  else // ( strand == LAGGING )
  {
    for ( int8_t i = 0 ; i < CODON_SIZE ; i++ )
    {
      if ( gen[Utils::mod((index-i),len)] != '1' ) // == and not != because we are on the complementary strand...
      {
        //_value += pow( 2, CODON_SIZE - i - 1 );
        _value += 1 << ( CODON_SIZE - i - 1 );
      }
    }
  }
}

Codon::Codon( gzFile backup_file )
{
  gzread( backup_file, &_value,   			sizeof(_value) );
}

// =================================================================
//                             Destructors
// =================================================================
Codon::~Codon( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================
void Codon::save( gzFile backup_file )
{
  gzwrite( backup_file, &_value,   			sizeof(_value) );
}
// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
