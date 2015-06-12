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
 
 
#ifndef AEVOL_GAUSSIAN_H__
#define AEVOL_GAUSSIAN_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <cstdlib>
#include <cmath>

#include <zlib.h>

namespace aevol {


// =================================================================
//                            Project Files
// =================================================================




// =================================================================
//                          Class declarations
// =================================================================





 
class Gaussian
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    inline Gaussian(double heigth, double mean, double width);
    inline Gaussian(const Gaussian &model);
    inline Gaussian(gzFile backup_file);
  
    // =================================================================
    //                             Destructor
    // =================================================================
    virtual inline ~Gaussian(void);
  
    // =================================================================
    //                              Accessors
    // =================================================================
    inline double get_height(void) const;
    inline double get_mean(void) const;
    inline double get_width(void) const;
    inline void   set_height(double height);
    inline void   set_mean(double mean);
    inline void   set_width(double width);
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    inline double compute_y(double x) const;
    
    inline void save(gzFile backup_file) const;
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    Gaussian(void)
    {
      printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
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
inline Gaussian::Gaussian(double height, double mean, double width)
{
  _height = height;
  _mean   = mean;
  _width  = width;
}

inline Gaussian::Gaussian(const Gaussian &model)
{
  _height = model._height;
  _mean   = model._mean;
  _width  = model._width;
}

inline Gaussian::Gaussian(gzFile backup_file)
{
  gzread(backup_file, &_height,  sizeof(_height));
  gzread(backup_file, &_mean,    sizeof(_mean));
  gzread(backup_file, &_width,   sizeof(_width));
}

// =====================================================================
//                               Destructor
// =====================================================================
inline Gaussian::~Gaussian(void)
{
}

// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline double Gaussian::get_height(void) const
{
  return _height;
}

inline double Gaussian::get_mean(void) const
{
  return _mean;
}

inline double Gaussian::get_width(void) const
{
  return _width;
}

inline void Gaussian::set_height(double height)
{
  _height = height;
}

inline void Gaussian::set_mean(double mean)
{
  _mean = mean;
}

inline void Gaussian::set_width(double width)
{
  _width = width;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
double Gaussian::compute_y(double x) const
{
  return _height * exp(-(x-_mean)*(x-_mean) / (2*_width*_width));
}

void Gaussian::save(gzFile backup_file) const
{
  gzwrite(backup_file, &_height, sizeof(_height));
  gzwrite(backup_file, &_mean, sizeof(_mean));
  gzwrite(backup_file, &_width, sizeof(_width));
}

} // namespace aevol

#endif // AEVOL_GAUSSIAN_H__
