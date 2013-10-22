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
#include "Test_ae_jumping_mt.h"



// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================




//############################################################################
//                                                                           #
//                         Class Test_ae_jumping_mt                          #
//                                                                           #
//############################################################################
CPPUNIT_TEST_SUITE_REGISTRATION( Test_ae_jumping_mt );

// ===========================================================================
//                               Static attributes
// ===========================================================================

// ===========================================================================
//                                  Constructors
// ===========================================================================
Test_ae_jumping_mt::Test_ae_jumping_mt( void )
{
}

// ===========================================================================
//                                  Destructors
// ===========================================================================
Test_ae_jumping_mt::~Test_ae_jumping_mt( void )
{
}

// ===========================================================================
//                                   Operators
// ===========================================================================

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void Test_ae_jumping_mt::setUp( void )
{
  obj1 = new ae_jumping_mt( 56873 );
  obj2 = new ae_jumping_mt( *obj1 );
}

void Test_ae_jumping_mt::tearDown( void )
{
  delete obj1;
  delete obj2;
}

void Test_ae_jumping_mt::test1( void )
{
  for ( int i = 0 ; i < 10 ; i++ )
  {
    obj1->random();
    obj2->random();
  }
  
  CPPUNIT_ASSERT( obj1->random() == obj2->random() );
  CPPUNIT_ASSERT( obj1->random( (int8_t) INT8_MAX ) == obj2->random( (int8_t) INT8_MAX ) );
  CPPUNIT_ASSERT( obj1->random( (int16_t) INT16_MAX ) == obj2->random( (int16_t) INT16_MAX ) );
  CPPUNIT_ASSERT( obj1->random( (int32_t) INT32_MAX ) == obj2->random( (int32_t) INT32_MAX ) );
  CPPUNIT_ASSERT( obj1->random( (int64_t) INT64_MAX ) == obj2->random( (int64_t) INT64_MAX ) );
  CPPUNIT_ASSERT( obj1->binomial_random( 435, 0.5 ) == obj2->binomial_random( 435, 0.5 ) );
  CPPUNIT_ASSERT( obj1->gaussian_random() == obj2->gaussian_random() );
}

void Test_ae_jumping_mt::test2( void )
{
  CPPUNIT_ASSERT( true );
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================

// ===========================================================================
//                              Non inline accessors
// ===========================================================================
