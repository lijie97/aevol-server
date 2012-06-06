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


/*! \class ae_simulation
    \brief This class contains all the data regarding a simulation (either main run or post-treatment)
*/
 
 
#ifndef  __AE_SIMULATION_H__
#define  __AE_SIMULATION_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_population.h>
#include <ae_environment.h>
#include <ae_stats.h>
#include <ae_logs.h>
#include <ae_tree.h>
#include <ae_dump.h>
#include <ae_rand_mt.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_param_loader;
class ae_param_overloader;




 
class ae_simulation : public ae_object
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_simulation( void );
    //~ ae_simulation( ae_param_overloader* param_overloader /*= NULL*/ );
    //~ ae_simulation( char* backup_file_name, bool to_be_run = true, ae_param_overloader* param_overloader = NULL );
    //~ ae_simulation( char* organism_file_name, ae_param_overloader* param_overloader = NULL );
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_simulation( void );
  
    // =================================================================
    //                         Accessors: getters
    // =================================================================
    inline ae_environment*  get_env( void )                       const;
    inline ae_population*   get_pop( void )                       const;
    inline ae_tree*         get_tree( void )                      const;
    inline int32_t          get_first_gener( void )               const;
    inline int32_t          get_num_gener( void )                 const;
  
    // =================================================================
    //                         Accessors: setters
    // =================================================================
    inline void set_env( ae_environment* env );
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    // Initialization
    void load_params( ae_param_loader* param_loader, ae_param_overloader* param_overloader = NULL );
    void load_backup( char* backup_file_name, bool to_be_run = true, ae_param_overloader* param_overloader = NULL );
    
    void run( void );
    
    void set_total_generations( int32_t nb_generations );
    
    void write_backup( void );
    void write_tree( void );
    void write_envir( void );
    
    inline bool fitness_is_composite( void );
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
    ae_rand_mt* alea;
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_simulation( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    ae_simulation( const ae_simulation &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    void read_from_backup( gzFile* backup_file );
    virtual void display( void ){};
  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    int32_t _num_gener;
    int32_t _first_gener;
    int32_t _last_gener;
  
    ae_environment*   _env;
    ae_population*    _pop;
    
    ae_stats* _stats;
    ae_tree*  _tree;
    ae_dump*  _dump;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline ae_environment* ae_simulation::get_env( void ) const
{
  return _env;
}

inline ae_population* ae_simulation::get_pop( void ) const
{
  return _pop;
}

inline ae_tree* ae_simulation::get_tree( void ) const
{
  return _tree;
}

inline int32_t ae_simulation::get_first_gener( void ) const
{
  return _first_gener;
}

inline int32_t ae_simulation::get_num_gener( void ) const
{
  return _num_gener;
}

inline void ae_simulation::set_env( ae_environment* env )
{
  _env = env;
}

inline bool ae_simulation::fitness_is_composite( void )
{
  return _env->fitness_is_composite();
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_SIMULATION_H__
