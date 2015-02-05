// Aevol - An in silico experimental evolution platform
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

#ifndef TEST_FUZZY
#define TEST_FUZZY

#include <cstdio>
#include <cstdlib>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <fuzzy.h>

namespace aevol {

using namespace CppUnit;
using namespace std;

class Test_fuzzy : public TestFixture
{
  CPPUNIT_TEST_SUITE(Test_fuzzy);
  CPPUNIT_TEST(truetest);
  CPPUNIT_TEST(falsetest);
  CPPUNIT_TEST_SUITE_END();
 public:
  void falsetest();
  void truetest();
};

} // namespace aevol

#endif // TEST_FUZZY
