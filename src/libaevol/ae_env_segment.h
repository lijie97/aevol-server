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


/*! \class ae_env_segment
  \brief
*/
 
 
#ifndef __AE_ENV_SEGMENT_H__
#define __AE_ENV_SEGMENT_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_macros.h>
#include <ae_enums.h>




// =================================================================
//                          Class declarations
// =================================================================





 
class ae_env_segment : public ae_object
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    inline ae_env_segment( double start, double stop, ae_env_axis_feature feature );
    inline ae_env_segment( const ae_env_segment& source );
    inline ae_env_segment( gzFile* backup_file );
  
    // =================================================================
    //                             Destructors
    // =================================================================
    inline virtual ~ae_env_segment( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    inline void save( gzFile* backup_file ) const;
    inline void load( gzFile* backup_file );
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
    double start;
    double stop;
    ae_env_axis_feature feature;
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_env_segment( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    //~ ae_env_segment( const ae_env_segment &model )
    //~ {
      //~ printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };

  
    // =================================================================
    //                           Protected Methods
    // =================================================================
  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
};




//##############################################################################
//                                                                             #
//                             Class ae_env_segment                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
//~ inline ae_env_segment::ae_env_segment( void )
//~ {
  //~ start   = X_MIN;
  //~ stop    = X_MAX;
  //~ feature = NEUTRAL;
//~ }

inline ae_env_segment::ae_env_segment( double start, double stop, ae_env_axis_feature feature )
{
  this->start   = start;
  this->stop    = stop;
  this->feature = feature;
}

inline ae_env_segment::ae_env_segment( const ae_env_segment& source )
{
  this->start   = source.start;
  this->stop    = source.stop;
  this->feature = source.feature;
}

inline ae_env_segment::ae_env_segment( gzFile* backup_file )
{
  load( backup_file );
}

// =================================================================
//                             Destructors
// =================================================================
inline ae_env_segment::~ae_env_segment( void )
{
}

// =====================================================================
//                          Accessors definitions
// =====================================================================

// =================================================================
//                            Public Methods
// =================================================================
inline void ae_env_segment::save( gzFile* backup_file ) const
{
  gzwrite( backup_file, &start, sizeof(start) );
  gzwrite( backup_file, &stop,  sizeof(stop) );
  int8_t tmp_feature = feature;
  gzwrite( backup_file, &tmp_feature, sizeof(tmp_feature) );
}

inline void ae_env_segment::load( gzFile* backup_file )
{
  gzread( backup_file, &start,  sizeof(start) );
  gzread( backup_file, &stop,   sizeof(stop) );
  int8_t tmp_feature;
  gzread( backup_file, &tmp_feature, sizeof(tmp_feature) );
  feature = (ae_env_axis_feature) tmp_feature;
}

// =================================================================
//                           Protected Methods
// =================================================================



#endif // __AE_ENV_SEGMENT_H__
