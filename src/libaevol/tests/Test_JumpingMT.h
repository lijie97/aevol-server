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


#ifndef TEST_AE_JUMPING_MT
#define TEST_AE_JUMPING_MT


// =================================================================
//                              Libraries
// =================================================================
#include <cstdio>
#include <cstdlib>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>



// =================================================================
//                            Project Files
// =================================================================
#include "JumpingMT.h"

namespace aevol {

// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================
using namespace CppUnit;
using namespace std;






class Test_JumpingMT : public TestFixture
{
  CPPUNIT_TEST_SUITE( Test_JumpingMT );
  CPPUNIT_TEST( test1 );
  CPPUNIT_TEST( test2 );
  CPPUNIT_TEST_SUITE_END();


  public :

    // =======================================================================
    //                                 Enums
    // =======================================================================

    // =======================================================================
    //                               Constructors
    // =======================================================================
    Test_JumpingMT( void );

    // =======================================================================
    //                               Destructors
    // =======================================================================
    virtual ~Test_JumpingMT( void );

    // =======================================================================
    //                            Accessors: getters
    // =======================================================================

    // =======================================================================
    //                            Accessors: setters
    // =======================================================================

    // =======================================================================
    //                                Operators
    // =======================================================================

    // =======================================================================
    //                              Public Methods
    // =======================================================================
    void setUp( void );
    void tearDown( void );
    void test1( void );
    void test2( void );

    // =======================================================================
    //                             Public Attributes
    // =======================================================================



  protected :

    // =======================================================================
    //                            Forbidden Constructors
    // =======================================================================
    /*Test_JumpingMT( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    Test_JumpingMT( const Test_JumpingMT &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =======================================================================
    //                              Protected Methods
    // =======================================================================

    // =======================================================================
    //                             Protected Attributes
    // =======================================================================
    JumpingMT* obj1;
    JumpingMT* obj2;
};


// ===========================================================================
//                              Getters' definitions
// ===========================================================================

// ===========================================================================
//                              Setters' definitions
// ===========================================================================

// ===========================================================================
//                          Inline Operators' definitions
// ===========================================================================

// ===========================================================================
//                          Inline functions' definition
// ===========================================================================

} // namespace aevol

#endif // TEST_AE_JUMPING_MT
