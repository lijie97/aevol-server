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


#ifndef __PARAM_LOADER_H__
#define __PARAM_LOADER_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>


// =================================================================
//                            Project Files
// =================================================================
#include <params.h>
#include <f_line.h>
#include <ae_params_mut.h>
#include <ae_jumping_mt.h>


// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;
class ae_environment;
class ae_individual;


class param_loader
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    param_loader( const char* file_name );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~param_loader( void );

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void read_file( void );
    void load( ae_exp_manager* exp_m, bool verbose = false );
    
    f_line* get_line( void ); 
    
    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    param_loader( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    param_loader( const param_loader &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    static void format_line( f_line*, char*, bool* );
    void interpret_line( f_line* line, int32_t cur_line );
    ae_individual* create_random_individual( ae_exp_manager* exp_m, ae_params_mut* param_mut, int32_t id ) const;
    ae_individual* create_random_individual_with_good_gene( ae_exp_manager* exp_m, ae_params_mut* param_mut, int32_t id ) const;
    ae_individual* create_clone( ae_individual* dolly, int32_t id ) const;
    
    

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    //~ ae_exp_manager* _exp_m;
    
    ae_jumping_mt* _prng;
    ae_jumping_mt* _mut_prng;
    ae_jumping_mt* _stoch_prng;
    
    char*   _param_file_name;
    FILE*   _param_file;
    
    params* _param_values;
    
    int32_t _cur_line;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

#endif // __param_loader_H__
