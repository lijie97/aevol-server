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

#ifndef FUZZY_H
#define FUZZY_H

#include <ae_list.h>
#include <ae_point_2d.h>
#include <ae_macros.h>

namespace aevol {

class exp_manager;

/// Triangle type-2 fuzzy sets.
///
/// This class provides management tools for "fuzzy sets" abscissa are
/// bound between MIN_X and MAX_X (defined elsewhere) A "fuzzy set"
/// should always have at least two points of abscissa MIN_X and MAX_X.
class fuzzy
{
  public:
    fuzzy();
    fuzzy(const fuzzy& f);
    fuzzy(const gzFile backup);
    virtual ~fuzzy();
    ae_list<ae_point_2d*>* get_points();
    void initialize();
    void simplify();
    void print_points() const;
    void add_triangle(double mean, double width, double height);
    void add(fuzzy* f);
    void sub(fuzzy* f);
    double get_geometric_area() const;
    double get_geometric_area(double begin, double end) const;
    double get_y(double x) const;
    double get_x(double y, ae_point_2d* left, ae_point_2d* right) const;
    ae_list_node<ae_point_2d*>* create_interpolated_point(double x);

    void add_upper_bound(double up);
    void add_lower_bound(double low);

    bool is_identical_to(const fuzzy* fs, double tolerance) const;

    void save(gzFile backup) const;
    void load(gzFile backup);

  protected:
    void _assert_order();
    exp_manager* _exp_m;
    ae_list<ae_point_2d*>* _points;

  private:
    // WARNING: The 2 following functions are private because they are
    // dangerous. They are optimized versions of their homonyms,
    // starting their search at the given list entry.

    // The following conditions are MANDATORY : 
    //  - list_entry must be a node of the point list (_points) 
    //  - the corresponding point's abscissa must be ≤ x
    double get_y(double x, ae_list_node<ae_point_2d*>* start) const;
    ae_list_node<ae_point_2d*>* create_interpolated_point(double x, ae_list_node<ae_point_2d*>* start);
};
} // namespace
#endif // FUZZY_H
