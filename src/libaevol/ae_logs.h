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


#ifndef __AE_LOGS_H__
#define __AE_LOGS_H__


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




// =================================================================
//                          Class declarations
// =================================================================



class ae_logs : public ae_object
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_logs( void );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_logs( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline FILE* get_log( ae_log_type log_type )   const;
    inline int8_t get_logs() const;
    inline bool  is_logged( ae_log_type log_type ) const;

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    //void save( gzFile backup_file ) const;
    void load( int8_t logs, int32_t num_gener  );
    void print_to_file( FILE* file ) const;
    
    void set_logs( int8_t logs );
    void flush( void );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_logs( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_logs( const ae_logs &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================
    void write_headers( void ) const;

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    int8_t  _logs; // Which logs are "turned on" (bitmap)
    FILE*   _transfer_log;
    FILE*   _rear_log;
    FILE*   _barrier_log;
    //FILE*   _param_modification_log;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline FILE* ae_logs::get_log( ae_log_type log_type ) const
{
  switch ( log_type )
  {
    case LOG_TRANSFER :
    {
      return _transfer_log;
    }
    case LOG_REAR :
    {
      return _rear_log;
    }
    case LOG_BARRIER :
    {
      return _barrier_log;
    }
    /*case LOG_LOADS :
    {
      return _param_modification_log;
    }*/
    default:
    {
      printf( "ERROR: unknown log_type in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    }
  }
}

inline int8_t ae_logs::get_logs() const
{
  return _logs;
}

inline bool ae_logs::is_logged( ae_log_type log_type ) const
{
  switch ( log_type )
  {
    case LOG_TRANSFER :
    {
      return ( _logs & LOG_TRANSFER );
    }
    case LOG_REAR :
    {
      return ( _logs & LOG_REAR );
    }
    case LOG_BARRIER :
    {
      return ( _logs & LOG_BARRIER );
    }
    /*case LOG_LOADS :
    {
      return ( _logs & LOG_LOADS );
    }*/
    default:
    {
      printf( "ERROR: unknown log_type in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    }
  }
}

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_LOGS_H__
