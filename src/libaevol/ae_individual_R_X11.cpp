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
#include <ae_individual_R_X11.h>

//##############################################################################
//                                                                             #
//                           Class ae_individual_R_X11                         #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_individual_R_X11::ae_individual_R_X11( const ae_individual_R_X11 &model ) :
  ae_individual( model ), ae_individual_R( model ), ae_individual_X11( model )
{
  //printf("ae_individual_R_X11( model )");
}

ae_individual_R_X11::ae_individual_R_X11( void )  :
ae_individual(), ae_individual_R(), ae_individual_X11()
{
  //printf("ae_individual_R_X11( void )");
}

ae_individual_R_X11::ae_individual_R_X11( ae_individual_R_X11* parent, int32_t id,
                                          ae_jumping_mt* mut_prng, ae_jumping_mt* stoch_prng ) :
        ae_individual( parent, id, mut_prng, stoch_prng ),
        ae_individual_R( parent, id, mut_prng, stoch_prng  ),
        ae_individual_X11( parent, id, mut_prng, stoch_prng  )
{
  //printf("ae_individual_R_X11( parent )");
}

ae_individual_R_X11::ae_individual_R_X11( gzFile backup_file ) :
ae_individual( backup_file ), ae_individual_R( backup_file ), ae_individual_X11( backup_file )
{
}

// =================================================================
//                             Destructors
// =================================================================
ae_individual_R_X11::~ae_individual_R_X11( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
