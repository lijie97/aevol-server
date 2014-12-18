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
double fuzzy::get_y(double x, list<point>::const_iterator begin) const {
  assert(x >= X_MIN and x <= X_MAX);
  assert(points.size() >= 2);

  // Get the first point having abscissa ≥ x
  list<point>::const_iterator p2 =
    find_if(begin, points.end(),
            [x](const point& m){return m.first >= x;});
  assert(p2 != points.end());
  
  if (p2->first == x) // If p2 has abscissa x, he's the guy
    return p2->second;
  else { // Otherwise interpolate
    list<point>::const_iterator p1 = prev(p2);
    return p1->second +
      (x - p1->first) * ((p2->second - p1->second) /
                         (p2->first - p1->first));
  }
}

double fuzzy::get_y(double x) const {
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
double fuzzy::get_x(const point& p1, const point& p2, double y) const {
  assert((p2.second <= y and y <= p1.second) or
  	 (p1.second <= y and y <= p2.second));
  assert(p1.second != p2.second);
  double x = p1.first + (y - p1.second) * ((p2.first - p1.first) /
                                           (p2.second - p1.second));
  assert(p1.first <= x and x <= p2.first);
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
void fuzzy::simplify() {
  for (list<point>::iterator p = points.begin();
       p != points.end() and p != prev(points.end()) and p != prev(points.end(), 2);
       ++p)
    if (p->first == next(p)->first and p->first == next(p,2)->first)
      points.erase(next(p), prev(find_if(p, points.end(), [p](const point& q){return q.first != p->first;})));
    else if (p->second == next(p)->second and p->second == next(p,2)->second)
      points.erase(next(p), prev(find_if(p, points.end(), [p](const point& q){return q.second != p->second;})));

  assert(points.begin()->first == X_MIN);
  assert(prev(points.end())->first == X_MAX);
}

/// Add a triangle to the fuzzy set.
/// \param mean abscissa of its apex
/// \param width of the side opposite to the apex
/// \param height ordinate of the apex
void fuzzy::add_triangle(double mean, double width, double height) {
  assert(points.begin()->first == X_MIN);
  assert(prev(points.end())->first == X_MAX);

  assert(width > 0.0);
  assert(X_MIN <= mean and mean <= X_MAX);
  assert(W_MIN <= width); // the maximum width depends on each individual
  // assert(MIN_H <= height and height <= MAX_H); Not necessarily because the concentration can be > 1

  const double threshold = 1e-15; // TODO: should it not be the machine epsilon?
                                  // if not, it should at least be a class constant

  if (fabs(width) < threshold or fabs(height) < threshold)
    return;

  list<point>::iterator p0, p1, p2;
  p0 = p1 = p2 = points.begin();
  
  double x0 = mean - width;
  double x1 = mean;
  double x2 = mean + width;

  // TODO: bugfix? if points on borders X_MIN,MAX, should not the ordinate be appropriately set?
  // TODO: create_interpolated_point should return an ITERATOR to point list
  if (x0 >= X_MIN)  p0 = create_interpolated_point(x0);
  p1 = create_interpolated_point(mean, p0);
  if (x2 <= X_MAX)  p2 = create_interpolated_point(x2, p1);

  // Update points with abscissas in (x0;x1)
  for (list<point>::iterator p = p0 ; p != p1 ; ++p)
    p->second += height * (p->first - x0) / (x1 - x0);
  
  // Update points with abscissas in (x0;x1)
  for (list<point>::iterator p = p1 ; p != p2 ; ++p)
    p->second += height * (x2 - p->first) / (x2 - x1);

  assert(points.begin()->first == X_MIN);
  assert(prev(points.end())->first == X_MAX);
  assert(is_increasing());
  return;
}

//************************************** proof reading mark

/// Add a fuzzy set to the current one. Should actually be called `operator+=()`.
///
/// Semantically speaking, we deal with fuzzy sets over the same
/// range. So adding two fuzzy sets sums up to adding the probability
/// functions.
void fuzzy::add(const fuzzy& fs) {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  // Add interpolated points to current fuzzy set so that
  // `fs.points` ⊂ `points`
  for (const point& q: fs.points)
    create_interpolated_point(q.first);

  // each point in `points` gets `fs`'s counterpart ordinate added to
  // it.
  for (point& p: points)
    p.second += fs.get_y(p.first);

  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
  // assert(is_increasing());
}

/// Substract to the current fuzzy set.
///
/// TODO: Dumb version (?), to be completed.
void fuzzy::sub(const fuzzy& fs) {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  for (const point& q: fs.points)
    create_interpolated_point(q.first);

  for (point& p: points)
    p.second -= fs.get_y(p.first);

  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
  // assert(is_increasing());
}

double fuzzy::get_geometric_area() const {
  return get_geometric_area(points.begin(), points.end());
}

/// Get integral of the absolute of probability function.
double fuzzy::get_geometric_area(list<point>::const_iterator begin,
				 list<point>::const_iterator end) const {
  double area = 0;
  for (list<point>::const_iterator p = begin ; p != points.end() and next(p) != end ; ++p)
    // Note that the area of a crossed trapezoid can be computed just
    // the same as a normal one provided the bases are counted
    // algebrically (+/-).
    area += fabs(((p->second) + next(p)->second) * (next(p)->first - p->first) / 2);
  return area;
}

/// TODO: test case with discontinuity
double fuzzy::get_geometric_area(double start_segment, double end_segment) const {

  // assert(start_segment < end_segment);

  // Fuzzy set first (resp last) point must be at x = X_MIN (resp x = X_MAX)
  assert(not points.empty());
  assert(points.begin()->first == X_MIN);
  assert(prev(points.end())->first == X_MAX);

  // We must have (X_MIN <= start_segment < end_segment <= X_MAX)
  assert(start_segment >= X_MIN and start_segment < end_segment and end_segment <= X_MAX);

  // ****************************************************************
  // TODO: check!!! (vld, 2014-12-16)
  // ****************************************************************

  list<point>::const_iterator begin = find_if(points.begin(), points.end(),
                                        [start_segment](const point& p){return p.first > start_segment;});
  list<point>::const_iterator end = find_if(begin, points.end(),
                                      [end_segment](const point& p){return p.first < end_segment;});

  assert(begin != points.end());
  assert(end != points.end());
  
  double first_part = fabs((get_y(start_segment) + begin->second) * (begin->first - start_segment) / 2.0);
  double last_part = fabs((get_y(end_segment) + end->second) * (end->first - end_segment) / 2.0);

  return first_part + get_geometric_area(begin, end) + last_part;
}

double area_test() {
  fuzzy f;
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
void fuzzy::add_upper_bound(double upper_bound) {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  for (list<point>::iterator p = points.begin() ; p != points.end() ; ++p) {
    if (p->second <= upper_bound)
      continue;
    // Surrounding segments must be cut if they fall down the upper
    // bound. Add cutting points.
    if (prev(p)->second < upper_bound)
      points.insert(p, point(get_x(*prev(p), *p, upper_bound), upper_bound));
    if (next(p)->second < upper_bound)
      points.insert(next(p), point(get_x(*p, *next(p), upper_bound), upper_bound));
    p->second = upper_bound;
  }
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
  // assert(is_increasing());
}

/// Probability function gets floored out.
///
/// `pf` := max(`pf`, `lower_bound`)
///
/// TODO: refactor with add_upper_bound (vld, 2014-12-17)
void fuzzy::add_lower_bound(double lower_bound) {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  for (list<point>::iterator p = points.begin() ; p != points.end() ; ++p) {
    if (p->second >= lower_bound)
      continue;
    // Surrounding segments must be cut if they fall down the upper
    // bound. Add cutting points.
    if (prev(p)->second > lower_bound)
      points.insert(p, point(get_x(*prev(p), *p, lower_bound), lower_bound));
    if (next(p)->second > lower_bound)
      points.insert(next(p), point(get_x(*p, *next(p), lower_bound), lower_bound));
    p->second = lower_bound;
  }

  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
  // assert(is_increasing());
}


bool fuzzy::is_identical_to(const fuzzy& fs, double tolerance ) const {
  // Since list::size() has constant complexity since C++ 11, checking
  // size is an inexpensive first step.
  if (points.size() != fs.points.size()) 
    return false;

  for (list<point>::const_iterator p = points.begin(), q = fs.points.begin() ;
       p != points.end() ; // no need to check q because both lists have same size
       ++p, ++q)
    if (fabs(p->first - q->first) > tolerance * (fabs(p->first) + fabs(q->first)) or
        fabs(p->second - q->second) > tolerance * (fabs(p->second) + fabs(q->second)))
      return false;
  return true;
}


void fuzzy::save(gzFile backup_file) const {
  int16_t nb_points = points.size();
  gzwrite(backup_file, &nb_points, sizeof(nb_points));

  for (const point& p : points)
    writepoint(p, backup_file);
}


void fuzzy::load(gzFile backup_file) {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  int16_t nb_points;
  gzread(backup_file, &nb_points, sizeof(nb_points));

  for (int16_t i = 0 ; i < nb_points ; i++)
    points.push_back(point(readpoint(backup_file)));

  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
}

list<point>::iterator fuzzy::create_interpolated_point(double x) {
  return create_interpolated_point(x, points.begin());
}

/// Find first point before abscissa `x`, starting from `start`.
///
/// `start_point` must refer to a point before abscissa `x`
list<point>::iterator fuzzy::create_interpolated_point(double x, std::list<point>::iterator start) {
  assert(points.begin()->first == X_MIN);
  assert(prev(points.end())->first == X_MAX);

  assert(x >= X_MIN and x <= X_MAX );
  assert(points.size() >= 2 );

  // TODO: probably denotes a logic error
  if (start->first <= x )
    start = points.begin();
  
  static int call_counter;
  ++call_counter;
  // cout << call_counter << std::endl;

  // get first point *greater* than x
  list<point>::iterator p = find_if(start, points.end(), [x](point& q){return q.first > x;});
  if (prev(p)->first == x)
    // point already in points
    return p;
  else {
    auto newp = points.insert(p, point(x, get_y(x)));
    // insert point *before* p
    return newp; // points.insert(p, point(x, get_y(x)));
  }
}

bool fuzzy::is_increasing() const {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  for (list<point>::const_iterator p = points.begin() ; p != prev(points.end(), 2) ; ++p)
    if (p->first > next(p)->first)
      return false;
  return true;
}

/// Set all points ordinate to 0
///
/// Used in ae_environment::apply_noise(). Not sure if it's useful.
void fuzzy::reset() {
  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);

  for (point& p: points)
    p.second = 0;

  // assert(points.begin()->first == X_MIN);
  // assert(prev(points.end())->first == X_MAX);
}
} // namespace aevol
