//
// Created by arrouan on 31/07/15.
//


#include <math.h>

#ifdef __BLAS__
#include <cblas.h>
#endif

#include <iostream>
#include "HybridFuzzy.h"

namespace aevol {


HybridFuzzy::HybridFuzzy( )
{
  _pheno_size = PHENO_SIZE;
  _points = new double[_pheno_size];
  for (int i = 0; i < _pheno_size; i++)
    _points[i] = 0.0;
}

HybridFuzzy::HybridFuzzy( const HybridFuzzy& model )
{
  _pheno_size = PHENO_SIZE;
  _points = new double[_pheno_size];

#ifdef __BLAS__
  cblas_dcopy(_pheno_size,model._points,1,_points,1);
#else
	for (int i=0; i < _pheno_size; i++)
		_points[i] = model._points[i];
#endif
}

HybridFuzzy::HybridFuzzy( gzFile backup_file )
{
  _pheno_size = PHENO_SIZE;
  _points = new double[_pheno_size];
  for (int i = 0; i < _pheno_size; i++)
    _points[i] = 0.0;

  load( backup_file );
}

// =================================================================
//                             Destructors
// =================================================================

HybridFuzzy::~HybridFuzzy( void )
{
  if (_points != NULL) delete [] _points;
  _points = NULL;
}

// =================================================================
//                            Public Methods
// =================================================================
void HybridFuzzy::simplify( )
{

}

void HybridFuzzy::reset( )
{
  for (int i = 0; i < _pheno_size; i++)
    _points[i] = 0.0;
}

void HybridFuzzy::add_triangle( double mean, double width, double height )
{
  if ( fabs(width) < 1e-15 || fabs(height) < 1e-15 ) return;

  // Compute triangle points' coordinates
  double x0 = mean - width;
  double x1 = mean;
  double x2 = mean + width;

  int ix0 = (int) (x0 * _pheno_size);
  int ix1 = (int) (x1 * _pheno_size);
  int ix2 = (int) (x2 * _pheno_size);

  if (ix0 < 0) ix0 = 0; else if (ix0 > (_pheno_size-1)) ix0 = _pheno_size-1;
  if (ix1 < 0) ix1 = 0; else if (ix1 > (_pheno_size-1)) ix1 = _pheno_size-1;
  if (ix2 < 0) ix2 = 0; else if (ix2 > (_pheno_size-1)) ix2 = _pheno_size-1;

  // Compute the first equation of the triangle
  double incY = height / (ix1 - ix0);
  int count = 1;
  // Updating value between x0 and x1
//  _points[ix0] = 0.0;
  for (int i = ix0+1; i < ix1; i++) {
    _points[i]+=incY*(count++);
  }

  _points[ix1]+= height;

  // Compute the second equation of the triangle
  incY = height / (ix2 - ix1);
  count = 1;

  // Updating value between x1 and x2
  for (int i = ix1+1; i < ix2; i++) {
    _points[i]+=(height-(incY*(count++)));
  }
}

void HybridFuzzy::add( const AbstractFuzzy& f )
{
  const HybridFuzzy to_add = (HybridFuzzy&)(f);
#ifdef __BLAS__
  cblas_daxpy(_pheno_size, 1.0, to_add.points(), 1, _points, 1);
#else
		for (int i = 0; i < _pheno_size; i++) {
			if (to_add._points[i] != 0) _points[i] = _points[i] + to_add._points[i];
		}
#endif
}

void HybridFuzzy::sub( const AbstractFuzzy& f )
{
  const HybridFuzzy to_sub = (HybridFuzzy&)(f);
#ifdef __BLAS__
  cblas_daxpy(_pheno_size, -1.0, to_sub.points(), 1, _points, 1);
#else
		for (int i = 0; i < _pheno_size; i++) {
			if (to_sub._points[i] !=0 ) _points[i] = _points[i] - to_sub._points[i];
		}
#endif
}

double HybridFuzzy::get_geometric_area( ) const
{
  return get_geometric_area(X_MIN,X_MAX);
}

double HybridFuzzy::get_geometric_area( double start_segment, double end_segment ) const
{
  double area = 0;

  int istart_segment = (int) (start_segment  * _pheno_size);
  int iend_segment = (int) (end_segment  * _pheno_size);

  if (istart_segment < 0) istart_segment = 0; else if (istart_segment > (_pheno_size-1)) istart_segment = _pheno_size-1;
  if (iend_segment < 0) iend_segment = 0; else if (iend_segment > (_pheno_size-1)) iend_segment = _pheno_size-1;
  for (int i = istart_segment; i < iend_segment; i++) {
    area+=((fabs(_points[i]) + fabs(_points[i+1])) / (2.0*_pheno_size));
  }

  return area;
}



bool HybridFuzzy::is_identical_to( const AbstractFuzzy& f, double tolerance  ) const
{
  const HybridFuzzy fs = (HybridFuzzy&)(f);
  // Since list::size() has constant complexity since C++ 11, checking
  // size is an inexpensive first step.
  if (get_pheno_size() != fs.get_pheno_size())
    return false;

  for (int i = 0; i < _pheno_size; i++)
    if (fabs(_points[i] - fs.points()[i]) > tolerance * (fabs(_points[i]) + fabs(fs.points()[i])) or
        fabs(_points[i] - fs.points()[i]) > tolerance * (fabs(_points[i]) + fabs(fs.points()[i])))
      return false;
  return true;
}


void HybridFuzzy::save( gzFile backup_file ) const
{
  gzwrite(backup_file, &_pheno_size, sizeof(_pheno_size));
  std::cout << __FILE__ << ":" << __LINE__ << ":" << gztell(backup_file) << std::endl;
  std::cout << __FILE__ << ":" << __LINE__ << ":" << _pheno_size << std::endl;

  for (int i = 0; i < _pheno_size; i++)
    gzwrite(backup_file, &_points[i], sizeof(_points[i]));
}


void HybridFuzzy::load( gzFile backup_file ) {
  gzread(backup_file, &_pheno_size, sizeof(_pheno_size));
  std::cout << __FILE__ << ":" << __LINE__ << ":" << gztell(backup_file) <<
                                                     std::endl;
  std::cout << __FILE__ << ":" << __LINE__ << ":" << _pheno_size << std::endl;
  for (int i = 0; i < _pheno_size; i++) {
    gzread(backup_file, &_points[i], sizeof(_points[i]));
  }
}


void HybridFuzzy::clip(clipping_direction direction, double bound) {

  if (direction == clipping_direction::min)
    for (int i = 0; i < _pheno_size; i++)
      _points[i] = _points[i] < bound ? bound : _points[i];
  else if (direction == clipping_direction::max)
    for (int i = 0; i < _pheno_size; i++)
      _points[i] = _points[i] > bound ? bound : _points[i];
}

void HybridFuzzy::add_point(double x, double y) {
  int ix = (int) ( x * _pheno_size);
  _points[ix] = y;
}

// =================================================================
//                           Protected Methods
// =================================================================
double HybridFuzzy::get_y( double x ) const
{
  int ix = (int) ( x * _pheno_size);

  double retValue = _points[ix];

  return retValue;
}

void HybridFuzzy::print() const
{
  for (int i = 0; i < _pheno_size; i++)
    if (_points[i]!=0) printf("[%d : %f] ",i,_points[i]);
  printf("\n");
}
}
