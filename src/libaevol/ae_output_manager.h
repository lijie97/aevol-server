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


/*! \class ae_output_manager
    \brief
*/


#ifndef __AE_OUPUT_MANAGER_H__
#define __AE_OUPUT_MANAGER_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_enums.h>
#include <ae_stats.h>
#include <ae_tree.h>
#include <ae_dump.h>
#include <ae_logs.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;






class ae_output_manager : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_output_manager( ae_exp_manager* exp_m );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_output_manager( void );

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    
    // Tree
    inline bool         get_record_tree( void ) const;
    inline int32_t      get_tree_step( void ) const;
    inline ae_tree_mode get_tree_mode( void ) const;
    inline ae_tree*     get_tree( void ) const;
  
    // Logs
    inline FILE* get_log( ae_log_type log_type )   const;
    inline bool  is_logged( ae_log_type log_type ) const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    inline void set_backup_step( int32_t backup_step );
    inline void set_big_backup_step( int32_t big_backup_step );
    inline void init_tree( ae_tree_mode _tree_mode, int32_t _tree_step );
    inline void set_dump_step( int32_t dump_step );

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void write_setup_file( gzFile* setup_file ) const;
    void write_setup_file( FILE* setup_file ) const;
    void load( gzFile* file, bool verbose );
    void load( FILE* file, bool verbose );
    void write_current_generation_outputs( void ) const;
    inline void flush( void );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_output_manager( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_output_manager( const ae_output_manager &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };


    // =================================================================
    //                           Protected Methods
    // =================================================================
    void write_tree( void ) const;

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;
    
    // Backups
    int32_t _backup_step;
    int32_t _big_backup_step;
    
    // Stats
    ae_stats* _stats;
    
    // Tree
    bool      _record_tree;
    ae_tree*  _tree;
    
    // Dumps
    bool      _make_dumps;
    int32_t   _dump_step;
    ae_dump*  _dump;
    
    // Logs
    ae_logs*  _logs;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// Tree
inline bool ae_output_manager::get_record_tree( void ) const
{
  return _record_tree;
}

inline int32_t ae_output_manager::get_tree_step( void ) const
{
  return _tree->get_tree_step();
}

inline ae_tree_mode ae_output_manager::get_tree_mode( void ) const
{
  return _tree->get_tree_mode();
}

inline ae_tree* ae_output_manager::get_tree( void ) const
{
  return _tree;
}

// Logs
inline FILE* ae_output_manager::get_log( ae_log_type log_type )   const
{
  return _logs->get_log( log_type );
}

inline bool  ae_output_manager::is_logged( ae_log_type log_type ) const
{
  return _logs->is_logged( log_type );
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void ae_output_manager::set_backup_step( int32_t backup_step )
{
  _backup_step = backup_step;
}

inline void ae_output_manager::set_big_backup_step( int32_t big_backup_step )
{
  _big_backup_step = big_backup_step;
}

inline void ae_output_manager::init_tree( ae_tree_mode _tree_mode, int32_t _tree_step )
{
  _record_tree = true;
  _tree = new ae_tree( _tree_mode, _tree_step );
}

inline void ae_output_manager::set_dump_step( int32_t dump_step )
{
  _make_dumps = true;
  _dump_step  = dump_step;
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_output_manager::flush( void )
{
  _stats->flush();
}

#endif // __AE_OUPUT_MANAGER_H__
