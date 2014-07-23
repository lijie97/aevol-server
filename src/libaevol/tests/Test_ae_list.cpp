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
#include "Test_ae_list.h"



// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================




//############################################################################
//                                                                           #
//                            Class Test_ae_list                             #
//                                                                           #
//############################################################################
CPPUNIT_TEST_SUITE_REGISTRATION( Test_ae_list );

// ===========================================================================
//                               Static attributes
// ===========================================================================

// ===========================================================================
//                                  Constructors
// ===========================================================================
Test_ae_list::Test_ae_list( void )
{
}

// ===========================================================================
//                                  Destructors
// ===========================================================================
Test_ae_list::~Test_ae_list( void )
{
}

// ===========================================================================
//                                   Operators
// ===========================================================================

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void Test_ae_list::setUp( void )
{
  int_list = new ae_list<int*>();
  for (int i = 0 ; i < INT_LIST_SIZE ; i++)
  {
    int_list->add(new int(i + 1));
  }
}

void Test_ae_list::tearDown( void )
{
  delete int_list;
}

void Test_ae_list::test_extract_sublist( void )
{
  // Initial check
  int int_list_size = INT_LIST_SIZE;
  CPPUNIT_ASSERT_EQUAL(int_list->get_nb_elts(), int_list_size);
  for (int i = 0 ; i < INT_LIST_SIZE ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(*(int_list->get_object(i)), i + 1);
  }


  // **************************************************************************
  // Extract the first element
  ae_list<int*>* int_list2 = int_list->extract_sublist(0, 1);

  // Check size
  CPPUNIT_ASSERT_EQUAL(int_list->get_nb_elts(), INT_LIST_SIZE - 1);
  CPPUNIT_ASSERT_EQUAL(int_list2->get_nb_elts(), 1);

  // Check content
  CPPUNIT_ASSERT_EQUAL(*(int_list2->get_object(0)), 1);
  for (int i = 0 ; i < INT_LIST_SIZE-1 ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(*(int_list->get_object(i)), i + 2);
  }
  CPPUNIT_ASSERT_EQUAL((void*)int_list->get_last()->get_next(), (void*)NULL);
  CPPUNIT_ASSERT_EQUAL((void*)int_list2->get_last()->get_next(), (void*)NULL);
  CPPUNIT_ASSERT_EQUAL((void*)int_list->get_first()->get_prev(), (void*)NULL);
  CPPUNIT_ASSERT_EQUAL((void*)int_list2->get_first()->get_prev(), (void*)NULL);
  delete int_list2;

  // **************************************************************************
  // Extract elements 10 through 13
  int_list2 = int_list->extract_sublist(10, 4);

  // Check size
  CPPUNIT_ASSERT_EQUAL(int_list->get_nb_elts(), INT_LIST_SIZE - 5);
  CPPUNIT_ASSERT_EQUAL(int_list2->get_nb_elts(), 4);

  // Check content
  for (int i = 0 ; i < 4 ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(*(int_list2->get_object(i)), i + 12);
  }
  for (int i = 0 ; i < 9 ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(*(int_list->get_object(i)), i + 2);
  }
  for (int i = 10 ; i < INT_LIST_SIZE-5 ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(*(int_list->get_object(i)), i + 6);
  }
  CPPUNIT_ASSERT_EQUAL((void*)int_list->get_last()->get_next(), (void*)NULL);
  CPPUNIT_ASSERT_EQUAL((void*)int_list2->get_last()->get_next(), (void*)NULL);
  CPPUNIT_ASSERT_EQUAL((void*)int_list->get_first()->get_prev(), (void*)NULL);
  CPPUNIT_ASSERT_EQUAL((void*)int_list2->get_first()->get_prev(), (void*)NULL);
  delete int_list2;

  // **************************************************************************
  // Extract last 4 elements
  int_list2 = int_list->extract_ending_sublist(4);

  // Check size
  CPPUNIT_ASSERT_EQUAL(int_list->get_nb_elts(), INT_LIST_SIZE - 9);
  CPPUNIT_ASSERT_EQUAL(int_list2->get_nb_elts(), 4);

  // Check content
  for (int i = 0 ; i < 4 ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(*(int_list2->get_object(i)), i + 17);
  }
  for (int i = 0 ; i < 9 ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(*(int_list->get_object(i)), i + 2);
  }
  for (int i = 10 ; i < INT_LIST_SIZE - 9 ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(*(int_list->get_object(i)), i + 6);
  }
  CPPUNIT_ASSERT_EQUAL((void*)int_list->get_last()->get_next(), (void*)NULL);
  CPPUNIT_ASSERT_EQUAL((void*)int_list2->get_last()->get_next(), (void*)NULL);
  CPPUNIT_ASSERT_EQUAL((void*)int_list->get_first()->get_prev(), (void*)NULL);
  CPPUNIT_ASSERT_EQUAL((void*)int_list2->get_first()->get_prev(), (void*)NULL);
  // delete int_list2;
}

void Test_ae_list::test2( void )
{
  CPPUNIT_ASSERT( true );
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================

// ===========================================================================
//                              Non inline accessors
// ===========================================================================
