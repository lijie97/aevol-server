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






class ae_output_manager
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
    
    // Backup
    inline int32_t	get_backup_step(void) const;
    inline int32_t	get_big_backup_step(void) const;
    
    // Tree
    inline bool         get_record_tree( void ) const;
    inline int32_t      get_tree_step( void ) const;
    inline ae_tree_mode get_tree_mode( void ) const;
    inline ae_tree*     get_tree( void ) const;
  
    // Logs
    inline FILE* get_log( ae_log_type log_type )   const;
    inline bool  is_logged( ae_log_type log_type ) const;
  
    // Stats
    inline bool get_compute_phen_contrib_by_GU( void ) const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    inline void set_backup_step( int32_t backup_step );
    inline void set_big_backup_step( int32_t big_backup_step );
    inline void init_tree( ae_exp_manager* exp_m, ae_tree_mode _tree_mode, int32_t _tree_step );
    inline void set_dump_step( int32_t dump_step );
    inline void set_compute_phen_contrib_by_GU( bool compute_phen_contrib_by_GU );
    inline void set_logs (int8_t logs);
  
    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void write_setup_file( gzFile setup_file ) const;
    void write_setup_file( FILE* setup_file ) const;
    void load( gzFile file, bool verbose, bool to_be_run  );
    void load( FILE* file, bool verbose, bool to_be_run );
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
    bool      _compute_phen_contrib_by_GU;
  
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

// Backup
inline int32_t ae_output_manager::get_backup_step(void) const
{
  return _backup_step;
}

inline int32_t ae_output_manager::get_big_backup_step(void) const
{
  return _big_backup_step;
}

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

// Stats
inline bool ae_output_manager::get_compute_phen_contrib_by_GU( void ) const
{
  return _compute_phen_contrib_by_GU;
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

inline void ae_output_manager::init_tree( ae_exp_manager* exp_m, ae_tree_mode _tree_mode, int32_t _tree_step )
{
  _record_tree = true;
  _tree = new ae_tree( exp_m, _tree_mode, _tree_step );
}

inline void ae_output_manager::set_dump_step( int32_t dump_step )
{
  _make_dumps = true;
  _dump_step  = dump_step;
}

inline void ae_output_manager::set_compute_phen_contrib_by_GU( bool compute_phen_contrib_by_GU )
{
  _compute_phen_contrib_by_GU = compute_phen_contrib_by_GU;
}

inline void ae_output_manager::set_logs (int8_t logs)
{
  _logs->set_logs(logs);
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
