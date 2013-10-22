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
#include <ae_phenotype.h>
#include <ae_individual.h>




//##############################################################################
//                                                                             #
//                              Class ae_phenotype                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_phenotype::ae_phenotype( void ) :
#ifdef __NO_X
  ae_fuzzy_set()
#elif defined __X11
  ae_fuzzy_set_X11()
#else
#error You must specify a graphic option
#endif
{
  _indiv = NULL;
}

ae_phenotype::ae_phenotype( ae_individual* indiv, const ae_phenotype &model ) :
#ifdef __NO_X
  ae_fuzzy_set( model )
#elif defined __X11
  ae_fuzzy_set_X11( model )
#else
#error You must specify a graphic option
#endif
{
  _indiv = indiv;
}

ae_phenotype::ae_phenotype( ae_individual* indiv ) :
#ifdef __NO_X
  ae_fuzzy_set()
#elif defined __X11
  ae_fuzzy_set_X11()
#else
#error You must specify a graphic option
#endif
{
  _indiv = indiv;
}

// =================================================================
//                             Destructors
// =================================================================
ae_phenotype::~ae_phenotype( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
