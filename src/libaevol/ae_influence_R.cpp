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

// =================================================================
//                            Project Files
// =================================================================
#include <ae_influence_R.h>

//##############################################################################
//                                                                             #
//                           Class ae_influence_R                              #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_influence_R::ae_influence_R( ae_rna* rna, ae_protein *protein, double enhancing_coef, double operating_coef )
{
  _enhancing_coef = enhancing_coef;
  _operating_coef = operating_coef;
  _protein        = protein;
  _rna            = rna;
}

// =================================================================
//                             Destructors
// =================================================================
ae_influence_R::~ae_influence_R( void )
{
  //printf("delete influence %p\n", this);
}

// =================================================================
//                            Public Methods
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
