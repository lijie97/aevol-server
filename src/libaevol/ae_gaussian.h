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
 
 
 #ifndef __AE_GAUSSIAN_H__
#define  __AE_GAUSSIAN_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <stdlib.h>
#include <math.h>



// =================================================================
//                            Project Files
// =================================================================




// =================================================================
//                          Class declarations
// =================================================================





 
class ae_gaussian : public ae_object
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    inline ae_gaussian( double heigth, double mean, double width );
    inline ae_gaussian( const ae_gaussian &model );
    inline ae_gaussian( gzFile* backup_file );
  
    // =================================================================
    //                             Destructor
    // =================================================================
    virtual inline ~ae_gaussian( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
    inline double get_height( void ) const;
    inline double get_mean( void ) const;
    inline double get_width( void ) const;
    inline void   set_height( double height );
    inline void   set_mean( double mean );
    inline void   set_width( double width );
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    inline double compute_y( double x );
    
    inline void save( gzFile* backup_file );
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_gaussian( void )
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
    double _height;
    double _mean;
    double _width; // In fact half-width to the inflexion points
};


// =====================================================================
//                               Constructors
// =====================================================================
inline ae_gaussian::ae_gaussian( double height, double mean, double width )
{
  _height = height;
  _mean   = mean;
  _width  = width;
}

inline ae_gaussian::ae_gaussian( const ae_gaussian &model )
{
  _height = model._height;
  _mean   = model._mean;
  _width  = model._width;
}

inline ae_gaussian::ae_gaussian( gzFile* backup_file )
{
  gzread( backup_file, &_height,  sizeof(_height) );
  gzread( backup_file, &_mean,    sizeof(_mean) );
  gzread( backup_file, &_width,   sizeof(_width) );
}

// =====================================================================
//                               Destructor
// =====================================================================
inline ae_gaussian::~ae_gaussian( void )
{
}

// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline double ae_gaussian::get_height( void ) const
{
  return _height;
}

inline double ae_gaussian::get_mean( void ) const
{
  return _mean;
}

inline double ae_gaussian::get_width( void ) const
{
  return _width;
}

inline void ae_gaussian::set_height( double height )
{
  _height = height;
}

inline void ae_gaussian::set_mean( double mean )
{
  _mean = mean;
}

inline void ae_gaussian::set_width( double width )
{
  _width = width;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
double ae_gaussian::compute_y( double x )
{
  return _height * exp( -(x-_mean)*(x-_mean) / (2*_width*_width) );
}

void ae_gaussian::save( gzFile* backup_file )
{
  gzwrite( backup_file, &_height, sizeof(_height) );
  gzwrite( backup_file, &_mean, sizeof(_mean) );
  gzwrite( backup_file, &_width, sizeof(_width) );
}


#endif // __AE_GAUSSIAN_H__
