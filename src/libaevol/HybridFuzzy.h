//
// Created by arrouan on 31/07/15.
//

#ifndef AEVOL_HYBRIDFUZZY_H
#define AEVOL_HYBRIDFUZZY_H

#define PHENO_SIZE 300

#include <vector>

#include "macros.h"
#include "Point.h"
#include "AbstractFuzzy.h"

namespace aevol {

class HybridFuzzy : public AbstractFuzzy
{
 public:
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  HybridFuzzy();
  HybridFuzzy(const HybridFuzzy& f);
  HybridFuzzy(const gzFile backup);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~HybridFuzzy();

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void save(gzFile backup) const;
  void load(gzFile backup);
  void reset();
  void simplify();
  void add_triangle(double mean, double width, double height);
  void add(const AbstractFuzzy& f);
  void sub(const AbstractFuzzy& f);
  void add_point(double x, double y);

  void clip(clipping_direction direction, double bound);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  double* points() const { return _points; };

  double get_geometric_area() const;

  double get_geometric_area(double start_segment, double end_segment) const;

  double get_y(double x) const;
  // get_x should be moved out of fuzzy class as it really applies to pair of points

  bool is_identical_to(const AbstractFuzzy& fs, double tolerance) const;

  int get_pheno_size() const { return _pheno_size; };

  void print() const;
  inline void clear()  {reset();};
  // ==========================================================================
  //                                 Setters
  // ==========================================================================

  // ==========================================================================
  //                                Operators
  // ==========================================================================

 protected:
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================



  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  double* _points = NULL;
  int _pheno_size = PHENO_SIZE;
};

double trapezoid_area(const Point& p1, const Point& p2);
} // namespace aevol


#endif //AEVOL_HYBRIDFUZZY_H
