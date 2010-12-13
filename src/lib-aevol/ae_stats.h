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


/** \class
 *  \brief
 */
 
 
 #ifndef __AE_STATS_H__
#define  __AE_STATS_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_common.h>
#include <ae_individual.h>




// =================================================================
//                          Class declarations
// =================================================================





 
class ae_stats : public ae_object
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_stats( void );
    ae_stats( int32_t num_gener );
    ae_stats( const char * ancstat_file_name ); // for post-treatments

 
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_stats( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
  
    // =================================================================
    //                            Public Methods
    // =================================================================

    void write_current_generation_statistics( void );
    void write_statistics_of_this_indiv( ae_individual * indiv, int32_t t ); // for post-treatments
    
    void flush( void );

    inline void write_header( FILE* file_name, const char* header );
    inline void write_header( FILE* file_name, const char* header, int8_t key );

    void write_headers( void );



    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_stats( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    ae_stats( const ae_stats &model )
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
    
    FILE** _stat_files_best;
    FILE** _stat_files_glob;
    char** _stat_files_best_names;
    char** _stat_files_glob_names;
    
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


void ae_stats::write_header( FILE* file_name, const char* header )
{
  if ( file_name != NULL) fprintf( file_name, "# %s\n", header );
}

void ae_stats::write_header( FILE* file_name, const char* header, int8_t key )
{
  if ( file_name != NULL) fprintf( file_name, "# %2d. %s\n", key, header );
}
        


#endif // __AE_STATS_H__
