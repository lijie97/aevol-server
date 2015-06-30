/// TODO: add unit tests
/// Why should there always be points (X_MIN, 0),(X_MAX, 0) ?
/// Many tests for double-type equality ==. Should't we check mod ε?

// Aevol - An in silico experimental evolution platform
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

#include "Fuzzy.h"

#include <cstdint>
#include <cassert>
#include <cmath>

#include <iterator>
#include <iostream>
#include <algorithm>

#include "Point.h"
#include "macros.h"

using std::list;
using std::prev;
using std::next;
using std::find_if;
using std::fabs;
using std::cout;
using std::endl;

namespace aevol {

/// Get probability of x belonging to fuzzy set.
///
/// If there is an actual point in `points` list with abscissa `x`
/// return it ordinate. Otherwise interpolate surrounding points.
///
/// TODO: use it! (vld, 2014-12-19)
/// 
double Fuzzy::get_y(double x, list<Point>::const_iterator begin) const {
  assert(x >= X_MIN and x <= X_MAX);
  assert(points.size() >= 2);

  // Get the first point having abscissa ≥ x
  list<Point>::const_iterator p2 =
    find_if(begin, points.end(),
            [x](const Point& m){return m.x >= x;});
  assert(p2 != points.end());
  
  if (p2->x == x) // If p2 has abscissa x, he's the guy
    return p2->y;
  else { // Otherwise interpolate
    list<Point>::const_iterator p1 = prev(p2);
    return p1->y +
      (x - p1->x) * ((p2->y - p1->y) /
                         (p2->x - p1->x));
  }
}

double Fuzzy::get_y(double x) const {
  return get_y(x, points.begin());
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
double Fuzzy::get_x(const Point& p1, const Point& p2, double y) const {
  assert((p2.y <= y and y <= p1.y) or
         (p1.y <= y and y <= p2.y));
  assert(p1.y != p2.y);
  double x = p1.x + (y - p1.y) * (p2.x - p1.x) /
                                 (p2.y - p1.y);
  assert((p2.x <= x and x <= p1.x) or
         (p1.x <= x and x <= p2.x));
  return x;
}

/// Remove (some) superfluous points.
///
/// When several (≥3) consecutive points in the fuzzy set are on the
/// same segment, the inner ones don't add information as they could
/// be interpolated from the outer ones. This function trims down the
/// `points` list from such points, only when they are on the same
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
/// TODO: double check if using this function is beneficial. Removed
/// points could then be recreated.
///
/// TODO: test with points
/// {(X_MIN,0), ((X_MIN+X_MAX)/2,0), (X_MAX,0)} → {(X_MIN,0), (X_MAX,0)}
/// d:= (X_MIN+X_MAX)/n {(X_MIN,0), ..., (k * (X_MIN + X_MAX) / n,0), ..., (X_MAX,0)} → {(X_MIN,0), (X_MAX,0)}
/// idem on non-null ordinate
/// idem on //y-axis
/// test with points starting/ending with constant
void Fuzzy::simplify() {
  // assert(invariant());

  for (list<Point>::iterator p = points.begin();
       p != points.end() and p != prev(points.end()) and p != prev(points.end(), 2);
       ++p)
    if (p->x == next(p)->x and p->x == next(p,2)->x)
      points.erase(next(p), prev(find_if(p, points.end(), [p](const Point& q){return q.x != p->x;})));
    else if (p->y == next(p)->y and p->y == next(p,2)->y)
      points.erase(next(p), prev(find_if(p, points.end(), [p](const Point& q){return q.y != p->y;})));

  // postcondition:
  // there are no 3 points that all share the same abscissas or that all share the same ordinates
  // all the points come from previous `points` list
  // assert(invariant());
}

/// Add a triangle to the fuzzy set.
/// \param mean abscissa of its apex
/// \param width of the side opposite to the apex
/// \param height ordinate of the apex
void Fuzzy::add_triangle(double mean, double width, double height) {
  // assert(invariant());

  assert(width > 0.0);
  assert(X_MIN <= mean and mean <= X_MAX);
  assert(W_MIN <= width); // the maximum width depends on each individual
  // assert(MIN_H <= height and height <= MAX_H); Not necessarily because the concentration can be > 1

  const double threshold = 1e-15; // TODO: should it not be the machine epsilon?
                                  // if not, it should at least be a class constant

  if (fabs(width) < threshold or fabs(height) < threshold)
    return;

  list<Point>::iterator p0, p1, p2;
  p0 = p1 = points.begin();
  p2 = prev(points.end());
  
  double x0 = mean - width;
  double x1 = mean;
  double x2 = mean + width;

  // TODO: bugfix? if points on borders X_MIN,MAX, should not the ordinate be appropriately set?
  // TODO: create_interpolated_point should return an ITERATOR to point list
  if (x0 >= X_MIN)  p0 = create_interpolated_point(x0);
  p1 = create_interpolated_point(mean, p0);
  if (x2 <= X_MAX)  p2 = create_interpolated_point(x2, p1);

  // Update points with abscissas in (x0;x1)
  for (list<Point>::iterator p = p0 ; p != p1 ; ++p)
    p->y += (p->x - x0) / (x1 - x0) * height;
  
  // Update points with abscissas in (x0;x1)
  for (list<Point>::iterator p = p1 ; p != p2 ; ++p)
    p->y += height * (x2 - p->x) / (x2 - x1);

  // assert(invariant());
  return;
}

/// Add a fuzzy set to the current one.
///
/// Should actually be called `operator+=()`.
///
/// Semantically speaking, we deal with fuzzy sets over the same
/// range. So adding two fuzzy sets sums up to adding the probability
/// functions.
void Fuzzy::add(const Fuzzy& fs) {
  // assert(invariant());

  // Add interpolated points to current fuzzy set so that
  // `fs.points` ⊂ `points`
  for (const Point& q: fs.points)
    create_interpolated_point(q.x);

  // each point in `points` gets `fs`'s counterpart ordinate added to
  // it.
  for (Point& p: points)
    p.y += fs.get_y(p.x);

  // assert(invariant());
}

/// Substract to the current fuzzy set.
///
/// TODO: Dumb version (?), to be completed.
void Fuzzy::sub(const Fuzzy& fs) {
  // assert(invariant());

  for (const Point& q: fs.points)
    create_interpolated_point(q.x);

  for (Point& p: points)
    p.y -= fs.get_y(p.x);

  // assert(invariant());
}

/// Absolute area between x-axis and segment [p1,p2].
///
/// The area of a crossed trapezoid can be computed just the same as a
/// normal one if the bases are counted algebrically (±).
double trapezoid_area(const Point& p1, const Point& p2) {
  return fabs((p1.y + p2.y) / 2.0 *
              (p2.x - p1.x));
}

double Fuzzy::get_geometric_area() const {
  return get_geometric_area(points.begin(), points.end());
}

/// Get integral of the absolute of probability function.
///
double Fuzzy::get_geometric_area(list<Point>::const_iterator begin,
                                 list<Point>::const_iterator end) const {
  // Precondition would be along the lines of:
  // assert(points.begin() <= begin < end < points.end());
  double area = 0;
  for (list<Point>::const_iterator p = begin ; next(p) != end ; ++p)
    area += trapezoid_area(*p, *next(p));
  return area;
}

double Fuzzy::get_geometric_area(double x_start, double x_stop) const {
  // assert(invariant());
  // Precondition: X_MIN ≤ x_start < x_stop ≤ X_MAX
  assert(X_MIN <= x_start and x_start < x_stop and x_stop <= X_MAX);

  // first point with abscissa ≥ x_start
  list<Point>::const_iterator begin = find_if(points.begin(), points.end(),
                                              [x_start](const Point& p){return p.x >= x_start;});
  // point following the last one with abscissa ≤ x_stop
  list<Point>::const_iterator end = find_if(begin, points.end(),
                                            [x_stop](const Point& p){return p.x > x_stop;});

  // area before begin
  double first_part = trapezoid_area(Point(x_start, get_y(x_start)), *begin);
  // area after prev(end)
  double last_part = trapezoid_area(*prev(end), Point(x_stop, get_y(x_stop)));

  return first_part + get_geometric_area(begin, end) + last_part;
}

// double Fuzzy::get_geometric_area(double start_segment, double end_segment) const {
//   // Precondition: X_MIN ≤ start_segment < end_segment ≤ X_MAX
//   assert(X_MIN <= start_segment and start_segment < end_segment and end_segment <= X_MAX);

//   Fuzzy copy(*this);

//   return copy.get_geometric_area(copy.create_interpolated_point(start_segment), next(copy.create_interpolated_point(end_segment)));
// }

double area_test() {
  Fuzzy f;
  f.add_triangle(0.5, 1.0, 0.5);
  double a = f.get_geometric_area(0.0,1.0);
  return a;
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
void Fuzzy::clip(clipping_direction direction, double bound) {
  // assert(invariant());

  for (list<Point>::iterator p = points.begin() ; p != points.end() ; ++p) {
    if (next(p) != points.end() and
        ((p->y < bound and bound < next(p)->y) or
         (p->y > bound and bound > next(p)->y))) { // ie if p and next(p) are across bound
      // insert interpolated point
      //           *after* p
      points.insert(next(p), Point(get_x(*p, *next(p), bound),
                                   bound));
      // could now fast forward over created point... TODO?
    }
    if ((direction == clipping_direction::min and p->y < bound) or
        (direction == clipping_direction::max and p->y > bound))
      p->y = bound;
  }

  // assert(invariant());
}


bool Fuzzy::is_identical_to(const Fuzzy& fs, double tolerance ) const {
  // Since list::size() has constant complexity since C++ 11, checking
  // size is an inexpensive first step.
  if (points.size() != fs.points.size()) 
    return false;

  for (list<Point>::const_iterator p = points.begin(), q = fs.points.begin() ;
       p != points.end() ; // no need to check q because both lists have same size
       ++p, ++q)
    if (fabs(p->x - q->x) > tolerance * (fabs(p->x) + fabs(q->x)) or
        fabs(p->y - q->y) > tolerance * (fabs(p->y) + fabs(q->y)))
      return false;
  return true;
}


void Fuzzy::save(gzFile backup_file) const {
  int16_t nb_points = points.size();
  gzwrite(backup_file, &nb_points, sizeof(nb_points));
  cout << __FILE__ << ":" << __LINE__ << ":" << gztell(backup_file) << endl;
  cout << __FILE__ << ":" << __LINE__ << ":" << nb_points << endl;

  for (const Point& p : points)
    writepoint(p, backup_file);
}


void Fuzzy::load(gzFile backup_file) {
  // assert(invariant());

  int16_t nb_points;
  gzread(backup_file, &nb_points, sizeof(nb_points));
  cout << __FILE__ << ":" << __LINE__ << ":" << gztell(backup_file) << endl;
  cout << __FILE__ << ":" << __LINE__ << ":" << nb_points << endl;
  for (int16_t i = 0 ; i < nb_points ; i++)
    points.push_back(Point(readpoint(backup_file)));

  // assert(invariant());
}

list<Point>::iterator Fuzzy::create_interpolated_point(double x) {
  return create_interpolated_point(x, points.begin());
}

/// Find first point before abscissa `x`, starting from `start`.
///
/// `start_point` must refer to a point before abscissa `x`
///
/// idempotent: creating existing point returns existing point
list<Point>::iterator Fuzzy::create_interpolated_point(double x, std::list<Point>::iterator start) {
  // assert(invariant());
  assert(x >= X_MIN and x <= X_MAX );

  // TODO: probably denotes a logic error
  if (start->x <= x )
    start = points.begin();
  
  // get first point with abscissa stricly greater than x
  list<Point>::iterator p = find_if(start, points.end(), [x](Point& q){return q.x > x;});
  if (prev(p)->x == x) {
    // point already in points
    // assert(invariant());
    return prev(p);
  }
  // insert point before p
  // assert(invariant());
  return points.insert(p, Point(x, get_y(x)));
}

/// Check that list of `points`' abscissas is (strictly) increasing.
bool Fuzzy::is_increasing() const {
  for (list<Point>::const_iterator p = points.begin() ; p != prev(points.end(), 2) ; ++p)
    if (p->x > next(p)->x)
      return false;
  return true;
}

/// Set all points ordinate to 0
///
// TODO <david.parsons@inria.fr> Not sure if it's useful.
void Fuzzy::reset() {
  // assert(invariant());

  for (Point& p: points)
    p.y = 0;

  // assert(invariant());
}
} // namespace aevol
