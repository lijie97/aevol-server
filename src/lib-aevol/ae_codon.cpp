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
#include <math.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_codon.h>
#include <ae_utils.h>




//##############################################################################
//                                                                             #
//                                Class ae_codon                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_codon::ae_codon( void )
{
  _value = -1;
}

ae_codon::ae_codon( const ae_codon &model )
{
  _value = model._value;
}

ae_codon::ae_codon( int8_t value )
{
  _value = value;
}

ae_codon::ae_codon( ae_dna* dna, ae_strand strand, int32_t index )
{
  const char* gen = dna->get_data();
  int32_t     len = dna->get_length();

  _value = 0;

  if ( strand == LEADING )
  {
    for ( int8_t i = 0 ; i < CODON_SIZE ; i++ )
    {
      if ( gen[ae_utils::mod((index+i),len)] == '1' )
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
      if ( gen[ae_utils::mod((index-i),len)] != '1' ) // == and not != because we are on the complementary strand...
      {
        //_value += pow( 2, CODON_SIZE - i - 1 );
        _value += 1 << ( CODON_SIZE - i - 1 );
      }
    }
  }
}

ae_codon::ae_codon( gzFile backup_file )
{
  gzread( backup_file, &_value,   			sizeof(_value) );
}

// =================================================================
//                             Destructors
// =================================================================
ae_codon::~ae_codon( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_codon::write_to_backup( gzFile backup_file )
{
  gzwrite( backup_file, &_value,   			sizeof(_value) );
}
// =================================================================
//                           Protected Methods
// =================================================================
