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
#include <string>



// =================================================================
//                            Project Files
// =================================================================
#include "Test_ae_list.h"

namespace aevol {

// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================




//############################################################################
//                                                                           #
//                            Class Test_ae_list                             #
//                                                                           #
//############################################################################
CPPUNIT_TEST_SUITE_REGISTRATION(Test_ae_list);

// ===========================================================================
//                               Static attributes
// ===========================================================================
// Don't change value (hard use in some tests)
const int Test_ae_list::INT_LIST_SIZE = 20;

// ===========================================================================
//                                  Constructors
// ===========================================================================
Test_ae_list::Test_ae_list(void)
{
}

// ===========================================================================
//                                  Destructors
// ===========================================================================
Test_ae_list::~Test_ae_list(void)
{
}

// ===========================================================================
//                                   Operators
// ===========================================================================

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void Test_ae_list::setUp(void)
{
  printf("Test_ae_list setUp\n");
  int_list = new ae_list<int*>();
  for (int i = 0 ; i < INT_LIST_SIZE ; i++)
  {
    int_list->add(new int(i + 1));
  }
}

void Test_ae_list::tearDown(void)
{
  int_list->erase(true);
  delete int_list;
}

void Test_ae_list::basic_tests1(void)
{
  printf("Test_ae_list basic_tests1\n");
  // Manually check the content of int_list (checks add(T*))
  ae_list_node<int*>* node = int_list->get_first();
  int32_t nb_elts = 0;
  while (node != NULL)
  {
    CPPUNIT_ASSERT_EQUAL(++nb_elts, *node->get_obj());
    node = node->get_next();
  }
  CPPUNIT_ASSERT_EQUAL(INT_LIST_SIZE, nb_elts);


  // Construct the same list using add_front
  ae_list<int*>* expected = new ae_list<int*>();
  for (int i = INT_LIST_SIZE ; i > 0 ; i--)
  {
    expected->add_front(new int(i));
  }

  assert_equal(expected, int_list, CPPUNIT_SOURCELINE());

  expected->erase(true);
  delete expected;


  // Check erase and is_empty
  int_list->erase(true);
  CPPUNIT_ASSERT(int_list->is_empty());
  CPPUNIT_ASSERT_EQUAL(0, int_list->get_nb_elts());
  CPPUNIT_ASSERT_EQUAL((void*)NULL, (void*)int_list->get_first());
  CPPUNIT_ASSERT_EQUAL((void*)NULL, (void*)int_list->get_last());
}

void Test_ae_list::basic_tests2(void)
{
  // Check get_object(int32_t pos)
  CPPUNIT_ASSERT_EQUAL((void*)NULL, (void*)int_list->get_object(-1));
  CPPUNIT_ASSERT_EQUAL((void*)NULL, (void*)int_list->get_object(INT_LIST_SIZE));
  for (int i = 0 ; i < INT_LIST_SIZE ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(i + 1, *(int_list->get_object(i)));
  }


  // Check get_node(int32_t pos)
  CPPUNIT_ASSERT_EQUAL((void*)NULL, (void*)int_list->get_node(-1));
  CPPUNIT_ASSERT_EQUAL((void*)NULL, (void*)int_list->get_node(INT_LIST_SIZE));
  for (int i = 0 ; i < INT_LIST_SIZE-1 ; i++)
  {
    CPPUNIT_ASSERT_EQUAL(i + 1, *(int_list->get_node(i)->get_obj()));
  }


  // Check remove(node) (remove first, last and arbitrary elt)
  int_list->remove(int_list->get_node(0), true, true);
  int_list->remove(int_list->get_node(INT_LIST_SIZE-2), true, true);
  int_list->remove(int_list->get_node(INT_LIST_SIZE/2), true, true);
  ae_list<int*>* expected = new ae_list<int*>();
  for (int i = 0 ; i < INT_LIST_SIZE/2 ; i++)
  {
    expected->add(new int(i + 2));
  }
  for (int i = INT_LIST_SIZE/2 ; i < INT_LIST_SIZE-3 ; i++)
  {
    expected->add(new int(i + 3));
  }

  assert_equal(expected, int_list, CPPUNIT_SOURCELINE());
  expected->erase(true);
  delete expected;
}

void Test_ae_list::test_extract_sublist(void)
{
  // Construct the same list as int_list
  ae_list<int*>* expected = new ae_list<int*>();
  for (int i = 0 ; i < INT_LIST_SIZE ; i++)
  {
    expected->add(new int(i + 1));
  }

  // Initial check
  assert_equal(expected, int_list, CPPUNIT_SOURCELINE());
  //printf("HERE %d %d\n", expected->get_nb_elts(), int_list->get_nb_elts());


  // **************************************************************************
  // Extract the first element
  ae_list<int*>* int_list2 = int_list->extract_sublist(0, 1);

  ae_list<int*>* expected2 = new ae_list<int*>();
  expected->remove(expected->get_first(), true, true);
  expected2->add(new int(1));

  assert_equal(expected, int_list, CPPUNIT_SOURCELINE());
  assert_equal(expected2, int_list2, CPPUNIT_SOURCELINE());

  expected2->erase(true);
  int_list2->erase(true);
  delete int_list2;

  // **************************************************************************
  // Extract elements 10 through 13
  int_list2 = int_list->extract_sublist(10, 4);

  for (int i = 0 ; i < 4 ; i++)
  {
    expected->remove(expected->get_node(10), true, true);
    expected2->add(new int(12 + i));
  }

  assert_equal(expected, int_list, CPPUNIT_SOURCELINE());
  assert_equal(expected2, int_list2, CPPUNIT_SOURCELINE());

  expected2->erase(true);
  int_list2->erase(true);
  delete int_list2;

  // **************************************************************************
  // Extract last 4 elements
  int_list2 = int_list->extract_ending_sublist(4);

  for (int i = 0 ; i < 4 ; i++)
  {
    expected->remove(expected->get_last(), true, true);
    expected2->add(new int(17 + i));
  }

  assert_equal(expected, int_list, CPPUNIT_SOURCELINE());
  assert_equal(expected2, int_list2, CPPUNIT_SOURCELINE());

  expected2->erase(true);
  int_list2->erase(true);
  delete int_list2;

  // **************************************************************************
  // Extract first 3 elements

  int_list2 = int_list->extract_starting_sublist(3);

  for (int i = 0 ; i < 3 ; i++)
  {
    expected2->add(new int(*expected->get_first()->get_obj()));
    expected->remove(expected->get_first(), true, true);
  }

  assert_equal(expected, int_list, CPPUNIT_SOURCELINE());
  assert_equal(expected2, int_list2, CPPUNIT_SOURCELINE());

  expected2->erase(true);
  int_list2->erase(true);
  delete int_list2;
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================
template <typename T>
void Test_ae_list::assert_equal(const ae_list<T>* expected,
                                const ae_list<T>* actual,
                                SourceLine SL)
{
  // Build message string
  char* msg = new char[256];
  sprintf(msg, "From %s:%d", SL.fileName().c_str(), SL.lineNumber());

  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg,
                               expected->get_nb_elts(),
                               actual->get_nb_elts());

  ae_list_node<T>* node1 = expected->get_first();
  ae_list_node<T>* node2 = actual->get_first();

  int32_t nb_elts = 0;
  while (node1 != NULL && node2 != NULL)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg, *node1->get_obj(), *node2->get_obj());
    nb_elts++;

    node1 = node1->get_next();
    node2 = node2->get_next();
  }

  CPPUNIT_ASSERT_EQUAL_MESSAGE(msg, expected->get_nb_elts(), nb_elts);
  delete msg;
}

void Test_ae_list::testfalse() {
  CPPUNIT_ASSERT(false);
}

// ===========================================================================
//                              Non inline accessors
// ===========================================================================

} // namespace aevol
