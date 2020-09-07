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

/// TODO: add unit tests
/// Why should there always be points (X_MIN, 0),(X_MAX, 0) ?
/// Many tests for ProteinConcentration-type equality ==. Should't we check mod ε?

#include "Vector_Fuzzy.h"

#include <cstdint>
#include <cassert>
#include <cmath>
#include <set>

#include <iterator>
#include <iostream>
#include <algorithm>

#include "Point.h"
#include "macros.h"
#include "Fuzzy.h"

using std::set;
using std::prev;
using std::next;
#ifndef __OPENMP_GPU
using std::find_if;
#endif
using std::fabs;
using std::cout;
using std::endl;

namespace aevol {

/// Get probability of x belonging to fuzzy set.
///
/// If there is an actual point in `points_` set with abscissa `x`
/// return it ordinate. Otherwise interpolate surrounding points.
///
/// TODO: use it! (vld, 2014-12-19)
///
ProteinConcentration Vector_Fuzzy::y(ProteinConcentration x, set<Point>::const_iterator begin) const {
  assert(x >= X_MIN and x <= X_MAX);
  assert(points_.size() >= 2);

  // Get the first point having abscissa ≥ x

  set<Point>::const_iterator p2 = points_.lower_bound(Point(x,0.0));
  assert(p2 != points_.end());

  if (p2->x == x) // If p2 has abscissa x, he's the guy
    return p2->y;
  else { // Otherwise interpolate
    set<Point>::const_iterator p1 = prev(p2);
    return p1->y +
      (x - p1->x) * ((p2->y - p1->y) /
                         (p2->x - p1->x));
  }
}

ProteinConcentration Vector_Fuzzy::y(ProteinConcentration x) const {
  return y(x, points_.begin());
}

/// Get abscissa of point interpolated between `p1` and `p2`, at
/// ordinate `y`.
///
/// \pre{`p1` and `p2` should not be on the same horizontal line}
/// otherwise any x would fit.
///
/// \pre{`y` should be between `p1` and `p2` ordinates} Despite the
/// fact that the reverse is mathematically sound, it's not supposed
/// to happend here.
ProteinConcentration Vector_Fuzzy::x(const Point& p1, const Point& p2, ProteinConcentration y) const {
  assert((p2.y <= y and y <= p1.y) or
         (p1.y <= y and y <= p2.y));
  assert(p1.y != p2.y);
  ProteinConcentration x = p1.x + (y - p1.y) * (p2.x - p1.x) /
                                 (p2.y - p1.y);
  assert((p2.x <= x and x <= p1.x) or
         (p1.x <= x and x <= p2.x));
  return x;
}

/// Remove (some) superfluous points_.
///
/// When several (≥3) consecutive points_ in the fuzzy set are on the
/// same segment, the inner ones don't add information as they could
/// be interpolated from the outer ones. This function trims down the
/// `points_` set from such points_, only when they are on the same
/// ordinate or on the same abscissa.
///
/// For instance, on the following probability function, A and B would
/// be removed:
/// \verbatim
///   ^      B
///   |   +--+-----+
///   |   |         \     ---
///   |  A+          \   /
///   |   |           \ /
///   |---+------------+-------+----->
/// \endverbatim
///
/// TODO: ProteinConcentration check if using this function is beneficial. Removed
/// points_ could then be recreated.
///
/// TODO: test with points_
/// {(X_MIN,0), ((X_MIN+X_MAX)/2,0), (X_MAX,0)} → {(X_MIN,0), (X_MAX,0)}
/// d:= (X_MIN+X_MAX)/n {(X_MIN,0), ..., (k * (X_MIN + X_MAX) / n,0), ..., (X_MAX,0)} → {(X_MIN,0), (X_MAX,0)}
/// idem on non-null ordinate
/// idem on //y-axis
/// test with points_ starting/ending with constant
void Vector_Fuzzy::simplify() {
  // assert(invariant());

  for (set<Point>::iterator p = points_.begin();
       p != points_.end() and p != prev(points_.end()) and p != prev(points_.end(), 2);
       ++p)
    if (p->x == next(p)->x and p->x == next(p,2)->x) {
      auto tmp = prev(
          find_if(p, points_.end(), [p](const Point& q){return q.x != p->x;})
      );
      points_.erase(next(p),tmp);
    }
    else if (p->y == next(p)->y and p->y == next(p,2)->y) {
      auto tmp = prev(
          find_if(p, points_.end(),
                       [p](const Point& q) {
                           return q.y != p->y;
                       })
      );
      points_.erase(next(p), tmp);
    }
  // postcondition:
  // there are no 3 points_ that all share the same abscissas or that all share the same ordinates
  // all the points_ come from previous `points_` set
  // assert(invariant());
}

/// Add a triangle to the fuzzy set.
/// \param mean abscissa of its apex
/// \param width of the side opposite to the apex
/// \param height ordinate of the apex
void Vector_Fuzzy::add_triangle(ProteinConcentration mean, ProteinConcentration width, ProteinConcentration height, bool verbose) {
  // assert(invariant());

  assert(width > 0.0);
  assert(X_MIN <= mean and mean <= X_MAX);
  assert(W_MIN <= width); // the maximum width depends on each individual
  // assert(MIN_H <= height and height <= MAX_H); Not necessarily because the concentration can be > 1

  const ProteinConcentration threshold = 1e-15; // TODO: should it not be the machine epsilon?
                                  // if not, it should at least be a class constant

  if (fabs(width) < threshold or fabs(height) < threshold)
    return;

  set<Point>::iterator p0, p1, p2;
  p0 = p1 = points_.begin();
  p2 = prev(points_.end());

  ProteinConcentration x0 = mean - width;
  ProteinConcentration x1 = mean;
  ProteinConcentration x2 = mean + width;

  // TODO: bugfix? if points on borders X_MIN,MAX, should not the ordinate be appropriately set?
  // TODO: create_interpolated_point should return an ITERATOR to point set
  if (x0 >= X_MIN)  p0 = create_interpolated_point(x0);
  p1 = create_interpolated_point(mean, p0);
  if (x2 <= X_MAX)  p2 = create_interpolated_point(x2, p1);

  // Update points with abscissas in (x0;x1)
  for (set<Point>::iterator p = p0 ; p != std::next(p1) ; ++p) {
      p->y += (p->x - x0) / (x1 - x0) * height;
  }

  // Update points with abscissas in (x0;x1)


            for (set<Point>::iterator p = std::next(p1); p != std::next(p2); ++p) {

                p->y += height * (x2 - p->x) / (x2 - x1);
            }
  return;
}

/// Add a fuzzy set to the current one.
///
/// Should actually be called `operator+=()`.
///
/// Semantically speaking, we deal with fuzzy sets over the same
/// range. So adding two fuzzy sets sums up to adding the probability
/// functions.
void Vector_Fuzzy::add(const AbstractFuzzy& f) {
  const Vector_Fuzzy fs = (Vector_Fuzzy&)(f);
  // assert(invariant());

  // Add interpolated points_ to current fuzzy set so that
  // `fs.points_` ⊂ `points_`
  for (const Point& q: fs.points_)
    create_interpolated_point(q.x);

  // each point in `points_` gets `fs`'s counterpart ordinate added to
  // it.
  for (Point p: points_)
    p.y += fs.y(p.x);

  // assert(invariant());
}

    void Vector_Fuzzy::add(Vector_Fuzzy* f) {
        for (const Point& q: f->points_) {
            create_interpolated_point(q.x);
        }

        // each point in `points_` gets `fs`'s counterpart ordinate added to
        // it.
        for (std::set<Point>::iterator p = points_.begin(); p != points_.end(); ++p)
            (*p).y += f->y((*p).x);

        // assert(invariant());
    }

/// Substract to the current fuzzy set.
///
/// TODO: Dumb version (?), to be completed.
void Vector_Fuzzy::sub(const AbstractFuzzy& f) {
  const Vector_Fuzzy fs = (Vector_Fuzzy&)(f);
  // assert(invariant());

          fs.print();

  for (const Point q: fs.points_) {
      create_interpolated_point(q.x);
  }

  for (Point p: points_)
    p.y -= fs.y(p.x);

  // assert(invariant());
}

    void Vector_Fuzzy::sub(Vector_Fuzzy* f) {
        for (const Point q: f->points_) {
            create_interpolated_point(q.x);
        }
        for (std::set<Point>::iterator p = points_.begin(); p != points_.end(); ++p)
            (*p).y -= f->y((*p).x);

        // assert(invariant());
    }

ProteinConcentration Vector_Fuzzy::get_geometric_area() const {
  return get_geometric_area(points_.begin(), points_.end());
}

/// Get integral of the absolute of probability function.
///
ProteinConcentration Vector_Fuzzy::get_geometric_area(set<Point>::const_iterator begin,
                             set<Point>::const_iterator end) const {
  // Precondition would be along the lines of:
  // assert(points_.begin() <= begin < end < points_.end());
  ProteinConcentration area = 0;
  for (set<Point>::const_iterator p = begin ; next(p) != end ; ++p)
    area += trapezoid_area(*p, *next(p));
  return area;
}

ProteinConcentration Vector_Fuzzy::get_geometric_area(ProteinConcentration x_start, ProteinConcentration x_stop) const {
  // assert(invariant());
  // Precondition: X_MIN ≤ x_start < x_stop ≤ X_MAX
  assert(X_MIN <= x_start and x_start < x_stop and x_stop <= X_MAX);

  // first point with abscissa ≥ x_start
  set<Point>::const_iterator begin =
                                    find_if(points_.begin(), points_.end(),
                                              [x_start](const Point& p){return p.x >= x_start;});
  // point following the last one with abscissa ≤ x_stop
  set<Point>::const_iterator end =
                                    find_if(begin, points_.end(),
                                            [x_stop](const Point& p){return p.x > x_stop;});
  // area before begin
  ProteinConcentration first_part = trapezoid_area(Point(x_start, y(x_start)), *begin);
  // area after prev(end)
  ProteinConcentration last_part = trapezoid_area(*prev(end), Point(x_stop, y(x_stop)));

  return first_part + get_geometric_area(begin, end) + last_part;
}

/// Probability function gets clipped either upwise ou downwise.
///
/// `pf` := min(`pf`, `upper_bound`)
///
///            X    above: removed               |
///           / \                                |
///          /   \               X      bound    |
/// --------o-----o-------------o-o--------      |
///        /       \   X       /   \             |
///       X         \ / \     /     \            |
///                  X   \   /       X           |
///                       \ /                    |
///      underneath: kept  X                     |

/// `pf` := max(`pf`, `lower_bound`)
void Vector_Fuzzy::clip(clipping_direction direction, ProteinConcentration bound) {
  // assert(invariant());

  for (set<Point>::iterator p = points_.begin() ; p != points_.end() ; ++p) {
    if (next(p) != points_.end() and
        ((p->y < bound and bound < next(p)->y) or
         (p->y > bound and bound > next(p)->y))) { // ie if p and next(p) are across bound
      // insert interpolated point
      //           *after* p

      points_.insert(next(p), Point(x(*p, *next(p), bound),
                                   bound));
      // could now fast forward over created point... TODO?
    }
    if ((direction == clipping_direction::min and p->y < bound) or
        (direction == clipping_direction::max and p->y > bound))
      p->y = bound;
  }

  // assert(invariant());
}


bool Vector_Fuzzy::is_identical_to(const AbstractFuzzy& f, ProteinConcentration tolerance ) const {
  const Vector_Fuzzy fs = (Vector_Fuzzy&)(f);
  // Since set::size() has constant complexity since C++ 11, checking
  // size is an inexpensive first step.
  if (points_.size() != fs.points_.size())
    return false;

  for (set<Point>::const_iterator p = points_.begin(), q = fs.points_.begin() ;
       p != points_.end() ; // no need to check q because both sets have same size
       ++p, ++q)
    if (fabs(p->x - q->x) > tolerance * (fabs(p->x) + fabs(q->x)) or
        fabs(p->y - q->y) > tolerance * (fabs(p->y) + fabs(q->y))) {
        return false;
    }
  return true;
}


void Vector_Fuzzy::save(gzFile backup_file) const {
  int16_t nb_points = points_.size();
  gzwrite(backup_file, &nb_points, sizeof(nb_points));
  cout << __FILE__ << ":" << __LINE__ << ":" << gztell(backup_file) << endl;
  cout << __FILE__ << ":" << __LINE__ << ":" << nb_points << endl;

  for (const Point& p : points_)
    writepoint(p, backup_file);
}


void Vector_Fuzzy::load(gzFile backup_file) {
  // assert(invariant());

  int16_t nb_points;
  gzread(backup_file, &nb_points, sizeof(nb_points));
  cout << __FILE__ << ":" << __LINE__ << ":" << gztell(backup_file) << endl;
  cout << __FILE__ << ":" << __LINE__ << ":" << nb_points << endl;
  for (int16_t i = 0 ; i < nb_points ; i++)
    points_.insert(Point(readpoint(backup_file)));

  // assert(invariant());
}

set<Point>::iterator Vector_Fuzzy::create_interpolated_point(ProteinConcentration x) {
  return create_interpolated_point(x, points_.begin());
}

/// Find first point before abscissa `x`, starting from `start`.
///
/// `start_point` must refer to a point before abscissa `x`
///
/// idempotent: creating existing point returns existing point
set<Point>::iterator Vector_Fuzzy::create_interpolated_point(ProteinConcentration x, std::set<Point>::iterator start) {
  // assert(invariant());
  assert(x >= X_MIN and x <= X_MAX);

  auto result = points_.insert(Point(x, y(x)));
  return result.first;
}

/// Check that set of `points_`' abscissas is (strictly) increasing.
bool Vector_Fuzzy::is_increasing() const {
  for (set<Point>::const_iterator p = points_.begin() ; p != prev(points_.end(), 2) ; ++p)
    if (p->x > next(p)->x)
      return false;
  return true;
}

/// Set all points ordinate to 0
///
// TODO <david.parsons@inria.fr> Not sure if it's useful.
void Vector_Fuzzy::reset() {
  for (Point p: points_)
    p.y = 0;

  // assert(invariant());
}

void Vector_Fuzzy::clear() {
  points_.clear();
}

void Vector_Fuzzy::add_point(ProteinConcentration x, ProteinConcentration y) {
  auto result = points_.insert(Point(x, y));

  if (!result.second)
    (*result.first).y += y;
}


void Vector_Fuzzy::print() const
{
  for (const Point& p : points_)
    printf("[%f : %f] ",p.x,p.y);
  printf("\n");
}
} // namespace aevol
