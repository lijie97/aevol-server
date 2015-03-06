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

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "ae_rna.h"
#include "ae_genetic_unit.h"
#include "ae_individual.h"

namespace aevol {



//##############################################################################
//                                                                             #
//                                 Class ae_rna                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_rna::ae_rna( ae_genetic_unit* gen_unit, const ae_rna &model )
{
  // Copy "trivial" attributes
  _gen_unit  = gen_unit;
  
  _strand             = model._strand;
  _pos                = model._pos;
  _transcript_length  = model._transcript_length;
  _basal_level        = model._basal_level;
  
  // Copy transcribed proteins
  // WARNING : Since this list do not "own" the proteins (they will not be deleted)
  //            proteins must NOT be CREATED here.

  // TODO : Not needed for the moment...
  // for (const auto& protein: model.transcribed_proteins)
  //   transcribed_proteins.push_back(protein);
}

ae_rna::ae_rna( ae_genetic_unit* gen_unit, ae_strand strand, int32_t pos, int8_t diff )
{
  _gen_unit  = gen_unit;
  _strand = strand;
  _pos    = pos;

  _transcript_length  = -1;
  _basal_level        = 1 - (double)diff / (PROM_MAX_DIFF + 1);
}

/*
ae_rna::ae_rna( ae_rna* parent )
{
  _gen_unit           = parent->_gen_unit;
  _strand             = parent->_strand;
  _pos                = parent->_pos;
  _transcript_length  = parent->_transcript_length;
  _basal_level        = parent->_basal_level;
}
*/

// =================================================================
//                             Destructors
// =================================================================
ae_rna::~ae_rna( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================
int32_t ae_rna::get_first_transcribed_pos( void ) const
{
  if ( _strand == LEADING )
  {
    return ae_utils::mod( _pos + PROM_SIZE, _gen_unit->get_dna()->get_length() );
  }
  else
  {
    return ae_utils::mod( _pos - PROM_SIZE, _gen_unit->get_dna()->get_length() );
  }
}

int32_t ae_rna::get_last_transcribed_pos( void ) const
{
  if ( _strand == LEADING )
  {
    return ae_utils::mod( _pos +  PROM_SIZE + _transcript_length - 1, _gen_unit->get_dna()->get_length() );
  }
  else
  {
    return ae_utils::mod( _pos - (PROM_SIZE + _transcript_length - 1), _gen_unit->get_dna()->get_length() );
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
