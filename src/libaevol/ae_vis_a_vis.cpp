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



// =================================================================
//                            Project Files
// =================================================================
#include <ae_vis_a_vis.h>




//##############################################################################
//                                                                             #
//                              Class ae_vis_a_vis                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_vis_a_vis::ae_vis_a_vis( const ae_dna* chrom_1, const ae_dna* chrom_2, int32_t i_1, int32_t i_2, ae_sense sense /*= DIRECT*/ )
{
  _chrom_1  = chrom_1;
  _chrom_2  = chrom_2;
  _i_1      = i_1;
  _i_2      = i_2;
  _sense    = sense;
}

ae_vis_a_vis::ae_vis_a_vis( const ae_vis_a_vis & orig )
{
  _chrom_1  = orig._chrom_1;
  _chrom_2  = orig._chrom_2;
  _i_1      = orig._i_1;
  _i_2      = orig._i_2;
  _sense    = orig._sense;
}

//~ ae_vis_a_vis::ae_vis_a_vis( const ae_vis_a_vis * orig )
//~ {
  //~ _chrom_1  = orig->_chrom_1;
  //~ _chrom_2  = orig->_chrom_2;
  //~ _i_1      = orig->_i_1;
  //~ _i_2      = orig->_i_2;
  //~ _sense    = orig->_sense;
//~ }

// =================================================================
//                             Destructors
// =================================================================
ae_vis_a_vis::~ae_vis_a_vis( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
