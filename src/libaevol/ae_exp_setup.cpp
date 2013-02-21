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


/*! \class ae_exp_setup
    \brief This class contains the experiment setup among which the environment
*/


// =================================================================
//                              Libraries
// =================================================================
#include <stdio.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>




// =================================================================
//                            Project Files
// =================================================================
#include <ae_exp_setup.h>
#include <ae_population.h>

#ifdef __X11
  #include <ae_population_X11.h>
#endif




//##############################################################################
//                                                                             #
//                              Class ae_exp_setup                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_exp_setup::ae_exp_setup( ae_exp_manager* exp_m )
{
  _exp_m  = exp_m;
  _env    = new ae_environment();
  _sel    = new ae_selection( exp_m );
}
  

// =================================================================
//                             Destructors
// =================================================================
ae_exp_setup::~ae_exp_setup( void )
{
  delete _env;
  delete _sel;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_exp_setup::save( gzFile env_file, gzFile sp_struct_file ) const
{
  // Save environment state
  _env->save( env_file );

  // Save spatial structure state
  if ( is_spatially_structured() )
  {
    _sel->get_spatial_structure()->save( sp_struct_file );
  }
}

void ae_exp_setup::load( gzFile exp_setup_file, gzFile env_file, gzFile sp_struct_file, bool verbose )
{
  // Retrieve environmental data
  printf( "  Loading environment..." );
  fflush( stdout );
  //~ delete _env;
  //~ _env = new ae_environment();
  _env->load( env_file );
  printf( "OK\n" );

  // Retrieve selection data
  printf( "  Loading selection context..." );
  fflush( stdout );
  _sel->load( exp_setup_file, sp_struct_file );
  printf( "OK\n" );
}

void ae_exp_setup::load( FILE* exp_setup_file, gzFile env_file, gzFile sp_struct_file, bool verbose )
{
  // Retrieve environmental data
  printf( "  Loading environment..." );
  fflush( stdout );
  //~ delete _env;
  //~ _env = new ae_environment();
  _env->load( env_file );
  printf( "OK\n" );

  // Retrieve selection data
  printf( "  Loading selection context..." );
  fflush( stdout );
  //~ delete _sel;
  //~ _sel = new ae_selection();
  _sel->load( exp_setup_file, sp_struct_file );
  printf( "OK\n" );
}




// =================================================================
//                           Protected Methods
// =================================================================
