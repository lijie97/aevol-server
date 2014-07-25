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


#ifndef TEST_AE_LIST
#define TEST_AE_LIST


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
#include <ae_list.h>




// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================
using namespace CppUnit;
using namespace std;






class Test_ae_list : public TestFixture
{
  CPPUNIT_TEST_SUITE( Test_ae_list );
  CPPUNIT_TEST( basic_tests1 );
  CPPUNIT_TEST( basic_tests2 );
  CPPUNIT_TEST( test_extract_sublist );
  CPPUNIT_TEST_SUITE_END();


  public :

    // =======================================================================
    //                                 Enums
    // =======================================================================

    // =======================================================================
    //                               Constructors
    // =======================================================================
    Test_ae_list( void );

    // =======================================================================
    //                               Destructors
    // =======================================================================
    virtual ~Test_ae_list( void );

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

    void basic_tests1( void );
    void basic_tests2( void );
    void test_extract_sublist( void );

    // =======================================================================
    //                             Public Attributes
    // =======================================================================



  protected :

    // =======================================================================
    //                            Forbidden Constructors
    // =======================================================================
    /*Test_ae_list( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    Test_ae_list( const Test_ae_list &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =======================================================================
    //                              Protected Methods
    // =======================================================================
    template <typename T> void assert_equal(ae_list<T>* expected,
                                            ae_list<T>* actual,
                                            SourceLine SL);
    //void assert_equal(ae_list<int*>* int_list1, ae_list<int*>* int_list2);

    // =======================================================================
    //                             Protected Attributes
    // =======================================================================
    ae_list<int*>* int_list;
    static const int INT_LIST_SIZE;
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


#endif // TEST_AE_LIST
