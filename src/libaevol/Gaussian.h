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

#ifndef AEVOL_GAUSSIAN_H_
#define AEVOL_GAUSSIAN_H_

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

class Gaussian {  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    Gaussian(double heigth, double mean, double width);
    Gaussian(const Gaussian& model);
    Gaussian(gzFile backup_file);
  
    // =================================================================
    //                             Destructor
    // =================================================================
    virtual ~Gaussian() {}
  
    // =================================================================
    //                              Accessors
    // =================================================================
    double get_height() const { return _height; }
    double get_mean() const { return _mean; }
    double get_width() const { return _width; }
    void   set_height(double height) { _height = height; }
    void   set_mean(double mean) { _mean = mean; }
    void   set_width(double width) { _width = width; }
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    double compute_y(double x) const { return _height * exp(-(x-_mean)*(x-_mean) / (2*_width*_width)); }
    void save(gzFile backup_file) const;
  
    // =================================================================
    //                           Public Attributes
    // =================================================================

  protected :
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    Gaussian() = delete;

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
inline Gaussian::Gaussian(double height, double mean, double width) {
  _height = height;
  _mean   = mean;
  _width  = width;
}

inline Gaussian::Gaussian(const Gaussian &model) {
  _height = model._height;
  _mean   = model._mean;
  _width  = model._width;
}

inline Gaussian::Gaussian(gzFile backup_file) {
  gzread(backup_file, &_height,  sizeof(_height));
  gzread(backup_file, &_mean,    sizeof(_mean));
  gzread(backup_file, &_width,   sizeof(_width));
}

// =====================================================================
//                               Destructor
// =====================================================================

// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       functions' definition
// =====================================================================
inline void Gaussian::save(gzFile backup_file) const {
  gzwrite(backup_file, &_height, sizeof(_height));
  gzwrite(backup_file, &_mean, sizeof(_mean));
  gzwrite(backup_file, &_width, sizeof(_width));
}

} // namespace aevol

#endif // AEVOL_GAUSSIAN_H_
