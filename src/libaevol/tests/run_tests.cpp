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
 
 
 
 
 
// ===========================================================================
//                               Include Libraries
// ===========================================================================
#include <cppunit/TestCase.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>



// ===========================================================================
//                             Include Project Files
// ===========================================================================
#include "Test_ae_jumping_mt.h"


// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================
using namespace CppUnit;



// ===========================================================================
//                         Declare Miscellaneous Functions
// ===========================================================================



int main( int argc, char* argv[] )
{
  // Print message
  cout << "Running regression tests";

  // Get the top level suite from the registry
  Test *suite = TestFactoryRegistry::getRegistry().makeTest();

  // Adds the test to the list of test to run
  TextUi::TestRunner runner;
  runner.addTest( suite );

  // Change the default outputter to a compiler error format outputter
  runner.setOutputter( new CompilerOutputter( &runner.result(), cerr ) );
  
  // Run the tests.
  bool wasSucessful = runner.run();

  // Return error code 1 if the one of test failed.
  return wasSucessful ? 0 : 1;
}



// ===========================================================================
//                         Define Miscellaneous Functions
// ===========================================================================
