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


/*! \class ae_stats
    \brief Manage statistics files
*/
 
 
#ifndef __AE_STATS_H__
#define __AE_STATS_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_individual.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;





 
class ae_stats : public ae_object
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_stats( ae_exp_manager* exp_m, const char * prefix = "stat", bool best_indiv_only = false );
    ae_stats( ae_exp_manager* exp_m, int32_t num_gener, const char * prefix = "stat", bool best_indiv_only = false, bool addition_old_stats = true, bool delete_old_stats = true );

 
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_stats( void );

    // =================================================================
    //                        Accessors: getters
    // =================================================================

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    inline void set_exp_m( ae_exp_manager* exp_m);
    
    // =================================================================
    //                            Public Methods
    // =================================================================

    void write_current_generation_statistics( void );
    void write_statistics_of_this_indiv( ae_individual* indiv, int32_t num_gener );
    
    void flush( void );

    void write_headers( bool ancstats_stats = false );



    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_stats( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_stats( const ae_stats &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
  
    // =================================================================
    //                           Protected Methods
    // =================================================================
    void init_data( void );
    void set_file_names( const char * prefix, bool best_indiv_only );
    void open_files( void );
    
    inline void write_header( FILE* file_name, const char* header );
    inline void write_header( FILE* file_name, const char* header, int8_t key );
  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;
    
    // 3D tables of stat files (FILE*) and their names (char*)
    // Dimensions are given by:
    //    * genetic unit (ALL_GU, CHROM or PLASMIDS)
    //    * BEST or GLOB
    //    * stat type (FITNESS_STATS, MUTATION_STATS, GENES_STATS, BP_STATS or REAR_STATS)
    // Files that are not wanted MUST have their name set to NULL.
    // The files themselves are also NULL because we don't fopen() them.
    FILE**** _stat_files;
    char**** _stat_files_names;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_stats::set_exp_m( ae_exp_manager* exp_m)
{
  _exp_m = exp_m;
}

inline void ae_stats::write_header( FILE* file_name, const char* header )
{
  if ( file_name != NULL) fprintf( file_name, "# %s\n", header );
}

inline void ae_stats::write_header( FILE* file_name, const char* header, int8_t key )
{
  if ( file_name != NULL) fprintf( file_name, "# %2d. %s\n", key, header );
}
        


#endif // __AE_STATS_H__
