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
// ****************************************************************************
 
 
#ifndef __AE_DUMP_H__
#define  __AE_DUMP_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;



#define _AE_DUMP_FILENAME_BUFSIZE INT32_C(256)


class ae_dump : public ae_object
{  
 public :
  
  // =================================================================
  //                             Constructors
  // =================================================================
  ae_dump( ae_exp_manager* exp_m );
  
  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~ae_dump( void )
  {
  }
  
  // =================================================================
  //                              Accessors
  // =================================================================
  
  // =================================================================
  //                            Public Methods
  // =================================================================
  void write_current_generation_dump( void );
  void write_fitness_total( void );
  void write_secretion_present( void );
  void write_fitness_metabolic( void );
  void write_secreted_amount( void );
  void write_individual_probes( void );
  
  // =================================================================
  //                           Public Attributes
  // =================================================================
  
  
  
  protected :
    FILE* current_file;
    char  filename_buffer[_AE_DUMP_FILENAME_BUFSIZE];
    
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_dump( const ae_dump &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
#endif // __AE_DUMP_H__
