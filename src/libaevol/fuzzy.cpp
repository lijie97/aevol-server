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

#include <cstdint>
#include <cassert>
#include <cmath>
#include <iterator>
#include <algorithm>

#include "fuzzy.h"
#include "point.h"
#include "ae_macros.h"

using std::list;
using std::prev;
using std::next;
using std::find_if;
using std::fabs;

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
/// \pre{`y` should be between `p1` and `p2` ordinates} because
/// despite the fast that it might be mathematically sound otherwise,
/// it suggest the user has mixed things up.
double Fuzzy::get_x(const Point& p1, const Point& p2, double y) const {
  assert((p2.y <= y and y <= p1.y) or
  	 (p1.y <= y and y <= p2.y));
  assert(p1.y != p2.y);
  double x = p1.x + (y - p1.y) * ((p2.x - p1.x) /
                                           (p2.y - p1.y));
  assert(p1.x <= x and x <= p2.x);
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
  for (list<Point>::iterator p = points.begin();
       p != points.end() and p != prev(points.end()) and p != prev(points.end(), 2);
       ++p)
    if (p->x == next(p)->x and p->x == next(p,2)->x)
      points.erase(next(p), prev(find_if(p, points.end(), [p](const Point& q){return q.x != p->x;})));
    else if (p->y == next(p)->y and p->y == next(p,2)->y)
      points.erase(next(p), prev(find_if(p, points.end(), [p](const Point& q){return q.y != p->y;})));

  assert(points.begin()->x == X_MIN);
  assert(prev(points.end())->x == X_MAX);
}

/// Add a triangle to the fuzzy set.
/// \param mean abscissa of its apex
/// \param width of the side opposite to the apex
/// \param height ordinate of the apex
void Fuzzy::add_triangle(double mean, double width, double height) {
  assert(points.begin()->x == X_MIN);
  assert(prev(points.end())->x == X_MAX);

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

  assert(points.begin()->x == X_MIN);
  assert(prev(points.end())->x == X_MAX);
  assert(is_increasing());
  return;
}

//************************************** proof reading mark

/// Add a fuzzy set to the current one. Should actually be called `operator+=()`.
///
/// Semantically speaking, we deal with fuzzy sets over the same
/// range. So adding two fuzzy sets sums up to adding the probability
/// functions.
void Fuzzy::add(const Fuzzy& fs) {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  // Add interpolated points to current fuzzy set so that
  // `fs.points` ⊂ `points`
  for (const Point& q: fs.points)
    create_interpolated_point(q.x);

  // each point in `points` gets `fs`'s counterpart ordinate added to
  // it.
  for (Point& p: points)
    p.y += fs.get_y(p.x);

  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
  // assert(is_increasing());
}

/// Substract to the current fuzzy set.
///
/// TODO: Dumb version (?), to be completed.
void Fuzzy::sub(const Fuzzy& fs) {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  for (const Point& q: fs.points)
    create_interpolated_point(q.x);

  for (Point& p: points)
    p.y -= fs.get_y(p.x);

  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
  // assert(is_increasing());
}

double Fuzzy::get_geometric_area() const {
  return get_geometric_area(points.begin(), points.end());
}

/// Get integral of the absolute of probability function.
double Fuzzy::get_geometric_area(list<Point>::const_iterator begin,
				 list<Point>::const_iterator end) const {
  double area = 0;
  for (list<Point>::const_iterator p = begin ; p != points.end() and next(p) != end ; ++p)
    // Note that the area of a crossed trapezoid can be computed just
    // the same as a normal one provided the bases are counted
    // algebrically (+/-).
    area += fabs(((p->y) + next(p)->y) * (next(p)->x - p->x) / 2);
  return area;
}

/// TODO: test case with discontinuity
double Fuzzy::get_geometric_area(double start_segment, double end_segment) const {

  // assert(start_segment < end_segment);

  // Fuzzy set first (resp last) point must be at x = X_MIN (resp x = X_MAX)
  assert(not points.empty());
  assert(points.begin()->x == X_MIN);
  assert(prev(points.end())->x == X_MAX);

  // We must have (X_MIN <= start_segment < end_segment <= X_MAX)
  assert(start_segment >= X_MIN and start_segment < end_segment and end_segment <= X_MAX);

  // ****************************************************************
  // TODO: check!!! (vld, 2014-12-16)
  // ****************************************************************

  list<Point>::const_iterator begin = find_if(points.begin(), points.end(),
                                        [start_segment](const Point& p){return p.x > start_segment;});
  list<Point>::const_iterator end = find_if(begin, points.end(),
                                      [end_segment](const Point& p){return p.x < end_segment;});

  assert(begin != points.end());
  assert(end != points.end());
  
  double first_part = fabs((get_y(start_segment) + begin->y) * (begin->x - start_segment) / 2.0);
  double last_part = fabs((get_y(end_segment) + end->y) * (end->x - end_segment) / 2.0);

  return first_part + get_geometric_area(begin, end) + last_part;
}

double area_test() {
  Fuzzy f;
  f.add_triangle(0.5, 1.0, 0.5);
  double a = f.get_geometric_area(0.0,1.0);
  return a;
}

/// Probability function gets "ceilinged out".
///
/// `pf` := min(`pf`, `upper_bound`)
///
///            X    above: removed		 |
///           / \				 |
///          /   \               X   upper_bound |
/// --------o-----o-------------o-o--------	 |
///        /       \   X       /   \		 |
///       X         \ / \     /     \		 |
///                  X   \   /       X		 |
///                       \ /			 |
///      underneath: kept  X			 |
///
/// TODO: prevent adding superfluous points by jumping over whole above zones
void Fuzzy::add_upper_bound(double upper_bound) {
  add_lower_bound(upper_bound, false);
}

/// Probability function gets floored out.
///
/// `pf` := max(`pf`, `lower_bound`)
///
/// TODO: rename function to clip and make arguments more explicit (vld, 2015-02-10)
void Fuzzy::add_lower_bound(double lower_bound, bool lower) {
  assert(points.begin()->x == X_MIN);
  assert(prev(points.end())->x == X_MAX);
  assert(is_increasing());

  for (list<Point>::iterator p = points.begin() ; p != points.end() ; ++p) {
    if (next(p) != points.end() and
        ((p->y < lower_bound and lower_bound < next(p)->y) or
         (p->y > lower_bound and lower_bound > next(p)->y))) { // ie if p and next(p) are across lower_bound
      // insert interpolated point
      //           *after* p
      points.insert(next(p), Point(get_x(*p, *next(p), lower_bound),
                                   lower_bound));
      // could now fast forward over created point... TODO?
    }
    if (lower     and p->y < lower_bound or
        not lower and p->y > lower_bound)
      p->y = lower_bound;
  }

  assert(points.begin()->x == X_MIN);
  assert(prev(points.end())->x == X_MAX);
  assert(is_increasing());
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

  for (const Point& p : points)
    writepoint(p, backup_file);
}


void Fuzzy::load(gzFile backup_file) {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  int16_t nb_points;
  gzread(backup_file, &nb_points, sizeof(nb_points));

  for (int16_t i = 0 ; i < nb_points ; i++)
    points.push_back(Point(readpoint(backup_file)));

  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
}

list<Point>::iterator Fuzzy::create_interpolated_point(double x) {
  return create_interpolated_point(x, points.begin());
}

/// Find first point before abscissa `x`, starting from `start`.
///
/// `start_point` must refer to a point before abscissa `x`
list<Point>::iterator Fuzzy::create_interpolated_point(double x, std::list<Point>::iterator start) {
  assert(points.begin()->x == X_MIN);
  assert(prev(points.end())->x == X_MAX);

  assert(x >= X_MIN and x <= X_MAX );
  assert(points.size() >= 2 );

  // TODO: probably denotes a logic error
  if (start->x <= x )
    start = points.begin();
  
  static int call_counter;
  ++call_counter;
  // cout << call_counter << std::endl;

  // get first point stricly greater than x and return its predecessor
  list<Point>::iterator p = find_if(start, points.end(), [x](Point& q){return q.x > x;});
  if (prev(p)->x == x)
    // point already in points
    return prev(p);
  else {
    auto newp = points.insert(p, Point(x, get_y(x)));
    // insert point *before* p
    return newp; // points.insert(p, point(x, get_y(x)));
  }
}

bool Fuzzy::is_increasing() const {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  for (list<Point>::const_iterator p = points.begin() ; p != prev(points.end(), 2) ; ++p)
    if (p->x > next(p)->x)
      return false;
  return true;
}

/// Set all points ordinate to 0
///
/// Used in ae_environment::apply_noise(). Not sure if it's useful.
void Fuzzy::reset() {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  for (Point& p: points)
    p.y = 0;

  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
}
} // namespace aevol
