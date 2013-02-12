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


/*! \class ae_params_record
    \brief Singleton. Contains all the parameters regarding outputs (backups, stats, logs, ...).
*/


#ifndef __AE_PARAMS_RECORD_H__
#define __AE_PARAMS_RECORD_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_enums.h>
#include <ae_logs.h>




// =================================================================
//                          Class declarations
// =================================================================






class ae_params_record : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_params_record( void );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_params_record( void );

    // =================================================================
    //                       Accessors: getters
    // =================================================================
    // Backups
    inline int32_t get_backup_step( void ) const;
    inline int32_t get_big_backup_step( void ) const;
    
    // Tree
    inline int32_t      get_tree_step( void ) const;
    inline bool         get_record_tree( void ) const;
    inline ae_tree_mode get_tree_mode( void ) const;
    
    // Logs
    inline bool is_logged( ae_log_type log_type ) const;
    inline FILE* get_log( ae_log_type log_type ) const;
    
    // Other
    inline bool    get_more_stats( void ) const;
    inline int32_t get_dump_period( void ) const;

    // =================================================================
    //                       Accessors: setters
    // =================================================================
    // Backups
    inline void set_backup_step( int32_t backup_step );
    inline void set_big_backup_step( int32_t big_backup_step );
    
    // Tree
    inline void set_tree_step( int32_t tree_step );
    inline void set_record_tree( bool record_tree );
    inline void set_tree_mode( ae_tree_mode tree_mode );
    
    // Logs
    inline void init_logs( int8_t to_be_logged );
    
    // Other
    inline void set_more_stats( bool more_stats );
    inline void set_dump_period( int32_t dump_period );

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void save( gzFile backup_file ) const;
    void load( gzFile backup_file, bool verbose = false );
    void print_to_file( FILE* file ) const;

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_params_record( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_params_record( const ae_params_record &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    // Backups
    int32_t _backup_step;
    int32_t _big_backup_step;
    
    // Tree
    int32_t      _tree_step;
    bool         _record_tree;
    ae_tree_mode _tree_mode;
    
    // Logs
    ae_logs* _logs;
    
    // Other
    bool    _more_stats;   // TODO : explain
    int32_t _dump_period;  // TODO : explain
};


// =====================================================================
//                          Getter definitions
// =====================================================================
// Backups
inline int32_t ae_params_record::get_backup_step( void ) const
{
  return _backup_step;
}

inline int32_t ae_params_record::get_big_backup_step( void ) const
{
  return _big_backup_step;
}

// Tree
inline int32_t ae_params_record::get_tree_step( void ) const
{
  return _tree_step;
}

inline bool ae_params_record::get_record_tree( void ) const
{
  return _record_tree;
}

inline ae_tree_mode ae_params_record::get_tree_mode( void ) const
{
  return _tree_mode;
}
    
// Logs
inline bool ae_params_record::is_logged( ae_log_type log_type ) const
{
  return _logs->is_logged( log_type );
}

inline FILE* ae_params_record::get_log( ae_log_type log_type ) const
{
  return _logs->get_log( log_type );
}

// Other
inline bool ae_params_record::get_more_stats( void ) const
{
  return _more_stats;
}

inline int32_t ae_params_record::get_dump_period( void ) const
{
  return _dump_period;
}


// =================================================================
//                          Setter definitions
// =================================================================
// Backups
inline void ae_params_record::set_backup_step( int32_t backup_step )
{
  _backup_step = backup_step;
}

inline void ae_params_record::set_big_backup_step( int32_t big_backup_step )
{
  _big_backup_step = big_backup_step;
}

// Tree
inline void ae_params_record::set_tree_step( int32_t tree_step )
{
  _tree_step = tree_step;
}

inline void ae_params_record::set_record_tree( bool record_tree )
{
  _record_tree = record_tree;
}

inline void ae_params_record::set_tree_mode( ae_tree_mode tree_mode )
{
  _tree_mode = tree_mode;
}
    
// Logs
inline void ae_params_record::init_logs( int8_t to_be_logged )
{
  _logs->set_logs( to_be_logged );
}

// Other
inline void ae_params_record::set_more_stats( bool more_stats )
{
  _more_stats = more_stats;
}

inline void ae_params_record::set_dump_period( int32_t dump_period )
{
  _dump_period = dump_period;
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_PARAMS_RECORD_H__
