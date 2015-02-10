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

#ifndef AEVOL_FUZZY_H
#define AEVOL_FUZZY_H

#include <list>

#include "ae_list.h" // TODO: not needed here; add where needed transitively
#include "ae_macros.h"
#include "point.h"

namespace aevol {

/// Triangular fuzzy sets.
///
/// This class provides management tools for "fuzzy sets" abscissa are
/// bound between X_MIN and X_MAX (defined in ae_macros.h) A "fuzzy
/// set" should always have at least two points of abscissa X_MIN and
/// X_MAX.
///
/// A fuzzy set holds elements in a set X together with a probability
/// function X → [0,1] which tells how likely element x ∈ X beholds to
/// the fuzzy set. With these triangular fuzzy sets, the probability
/// function is a finite sum of isosceles triangles.
///
/// The current class models fuzzy sets over range X = [X_MIN; X_MAX]
/// by representing the probability function as the list of singular
/// points on its graph.
///
/// \verbatim
/// \\code{.unparsed}
///           ^
///       y2  +...              X                                            ....
///           |...             / \                                           ....
///       y6  +...            /   \               X                          ....
///           |...           /     \             / \                         ....
///       y4  +...          /       \   X       /   \                        ....
///           |...         /         \ / \     /     \                       ....
///      y3,y9+...        /           X   \   /       \             X        ....
///           |...       /                 \ /         \           / \       ....
///       y5  +...      /                   X           \         /   \      ....
///           |...     /                                 \       /     \     ....
///       0   +--X----X---------|-----|-|---|-----|-------X-----X-------X----X--->
///            X_MIN x1        x2    x3 x4 x5    x6      x7    x8  x9  x10 X_MAX
/// \\endcode
/// \endverbatim
/// fs.points would hold the list {(X_MIN,0),(x1,y1),...,(x10,y10)(X_MAX,0)}
class Fuzzy
{
 public:
  Fuzzy(): points({Point(X_MIN, 0.0), Point(X_MAX, 0.0)}) {};
  // fuzzy();
  Fuzzy(const Fuzzy& f): points(f.points) {};
  Fuzzy(const gzFile backup) {load(backup);};
  virtual ~Fuzzy() {};
  // breaks encapsulation
  // TODO: remove foreign dependency and nuke this function
  std::list<Point>& get_points() {return points;};
  // called from ae_environment::apply_noise(), should otherwise be private
  // TODO: remove call from apply_noise() and make create_interpolated_point protected
  std::list<Point>::iterator create_interpolated_point(double x);
  std::list<Point>::iterator create_interpolated_point(double x, std::list<Point>::iterator start);
  void simplify();
  void add_triangle(double mean, double width, double height);
  void add(const Fuzzy& f);
  void sub(const Fuzzy& f);
  double get_geometric_area() const;
  double get_geometric_area(std::list<Point>::const_iterator begin,
                            std::list<Point>::const_iterator end) const;
  double get_geometric_area(double start_segment, double end_segment) const;
  double get_y(double x, std::list<Point>::const_iterator begin) const;
  double get_y(double x) const;
  // get_x should be moved out of fuzzy class as it really applies to pair of points
  double get_x(const Point& left, const Point& right, double y) const;
  

  void add_upper_bound(double up);
  void add_lower_bound(double low, bool lower = true);

  bool is_identical_to(const Fuzzy& fs, double tolerance) const;

  void save(gzFile backup) const;
  void load(gzFile backup);
  void reset();

 protected:
  bool is_increasing() const;
  std::list<Point> points;
};
} // namespace aevol
#endif // AEVOL_FUZZY_H
