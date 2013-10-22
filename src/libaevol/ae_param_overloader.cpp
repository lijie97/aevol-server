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
#include <ae_common.h>
#include <ae_param_overloader.h>
#include <ae_param_loader.h>





//##############################################################################
//                                                                             #
//                          Class ae_param_overloader                          #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_param_overloader::ae_param_overloader( void )
{
  _nb_overloaded = 0;
}

// =================================================================
//                             Destructors
// =================================================================
ae_param_overloader::~ae_param_overloader( void )
{
  for( int16_t i = 0 ; i < _nb_overloaded ; i++ )
  {
    delete [] _tokens_overloaded[i];
  }
}

// =================================================================
//                            Public Methods
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
void ae_param_overloader::overload_params( void ) const
{
  for( int16_t i = 0 ; i < _nb_overloaded ; i++ )
  {
    bool interpretable;
    f_line* line = new f_line();
    
    printf( "Overloading : %s\n", _tokens_overloaded[i] );
    ae_param_loader::format_line( line, _tokens_overloaded[i], &interpretable );
    if ( ! interpretable )
    {
      printf( "ERROR in overloaded parameter : \"%s\"\n", _tokens_overloaded[ i ] );
      exit( EXIT_FAILURE );
    }
    
    ae_param_loader::interpret_line( line, 0 );
    delete line;
  }
}

void ae_param_overloader::write_log( FILE* log_file ) const
{
  fprintf( log_file, "NB_PARAM_OVERLOADED %i\n", _nb_overloaded );
  for( int16_t i = 0 ; i < _nb_overloaded ; i++ )
  {
    fprintf( log_file, "    %s\n", _tokens_overloaded[i] );
  }
  fprintf( log_file, "\n");
}
