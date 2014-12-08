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

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "fuzzy.h"
#include "point.h"
#include "ae_macros.h"

namespace aevol {

fuzzy::fuzzy() {
  _points = new ae_list<point*>();

  initialize();
}

fuzzy::fuzzy(const fuzzy &model) {
  _points = new ae_list<point*>();

  ae_list_node<point*>* model_point_parser = model._points->get_first();
  point* model_point;

  while (model_point_parser != NULL) {
    model_point = model_point_parser->get_obj();

    _points->add(new point(*model_point));

    model_point_parser = model_point_parser->get_next();
  }
}

fuzzy::fuzzy(gzFile backup_file) {
  _points = new ae_list<point*>();

  load(backup_file);
}

fuzzy::~fuzzy() {
  if (_points != NULL) {
    _points->erase(true);
    delete _points;
  }
}

ae_list<point*>* fuzzy::get_points() {
  return _points;
}

void fuzzy::initialize() {
  _points->erase(true);
  _points->add(new point(X_MIN, 0.0));
  _points->add(new point(X_MAX, 0.0));
}

double fuzzy::get_y(double x) const {
  return get_y(x, _points->get_first());
}

double fuzzy::get_x(double y, point* left_point, point* right_point) const {
  assert(left_point != NULL);
  assert(right_point != NULL);
  assert(left_point->second != right_point->second);
  assert((y >= right_point->second && y <= left_point->second) || (y >= left_point->second && y <= right_point->second));

  double tmp = left_point->first + (right_point->first - left_point->first) * (y - left_point->second) / (right_point->second - left_point->second);

  if (tmp < left_point->first) return left_point->first;
  else if (tmp > right_point->first) return right_point->first;
  else return tmp;
}

ae_list_node<point*>* fuzzy::create_interpolated_point(double x) {
  return create_interpolated_point(x, _points->get_first());
}

void fuzzy::simplify() {
  // Removes useless points like A or B:
  //
  //  ^      B
  //  |   +--+-----+
  //  |   |         \     ---
  //  |  A+          \   /
  //  |   |           \ /
  //  |---+------------+-------+----->

  //  printf("initial nb of points %d\n",_points.get_count());

  ae_list_node<point*>* node = _points->get_first();
  ae_list_node<point*>* next_node = NULL;
  ae_list_node<point*>* other_node = NULL;
  point* pt = NULL;
  point* next_pt = NULL;
  point* other_pt = NULL;
  double x, y, next_x, next_y, other_x, other_y;

  while (node != NULL) {
    pt = node->get_obj();
    x = pt->first;
    y = pt->second;
    next_node = node->get_next();
    if (next_node != NULL) {
      next_pt = next_node->get_obj();
      next_x = next_pt->first;
      next_y = next_pt->second;
      if (next_x == x) {
        // discontinuity on x
        // is the next point also on the same x ?
        other_node = next_node->get_next();
        while (other_node != NULL) {
          other_pt = other_node->get_obj();
          other_x = other_pt->first;
          other_y = other_pt->second;
          if (other_x == x) {
            _points->remove(next_node, true, true);
            next_node = other_node;
          }
          else {
            break;
          }
          other_node = other_node->get_next();
        }
      }
    }
    node = next_node;
  }


  node = _points->get_first();
  next_node = NULL;
  other_node = NULL;
  pt = NULL;
  next_pt = NULL;
  other_pt = NULL;

  while (node != NULL) {
    pt = node->get_obj();
    x = pt->first;
    y = pt->second;
    next_node = node->get_next();
    if (next_node != NULL) {
      next_pt = next_node->get_obj();
      next_x = next_pt->first;
      next_y = next_pt->second;
      if (next_y == y) {
        // plateau on y
        // is the next point also on the same y ?
        other_node = next_node->get_next();
        while (other_node != NULL) {
          other_pt = other_node->get_obj();
          other_x = other_pt->first;
          other_y = other_pt->second;
          if (other_y == y) {
            _points->remove(next_node, true, true);
            next_node = other_node;
            next_pt = next_node->get_obj();
          }
          else break;
          other_node = other_node->get_next();
        }
      }
    }
    node = next_node;
  }

  assert(is_increasing());
}

void fuzzy::print_points() const {
  ae_list_node<point*>* point_node = _points->get_first();
  point* m;

  while (point_node != NULL) {
    m = point_node->get_obj();
    printf("  (%f, %f)\n", m->first, m->second);
    point_node = point_node->get_next();
  }
}

void fuzzy::add_triangle(double mean, double width, double height) {
  assert(X_MIN <= mean   && mean <= X_MAX);
  assert(W_MIN <= width); // the maximum width depends on each individual
  // assert(MIN_H <= height && height <= MAX_H); Not necessarily because the concentration can be > 1

  if (fabs(width) < 1e-15 || fabs(height) < 1e-15) return;

  // Compute triangle points' coordinates
  double x0 = mean - width;
  double x1 = mean;
  double x2 = mean + width;
  //~ double y0 = 0.0;
  double y1 = height;
  //~ double y2 = 0.0;

  // We will need at least a point at abscissa x0, x1 and x2,
  // Create them if they don't exist or update them if they do.
  ae_list_node<point*>* node_0 = NULL;
  ae_list_node<point*>* node_1 = NULL;
  ae_list_node<point*>* node_2 = NULL;

  if (x0 >= X_MIN)  node_0 = create_interpolated_point(x0);
                      node_1 = create_interpolated_point(x1, node_0);
  if (x2 <= X_MAX)  node_2 = create_interpolated_point(x2, node_1);


  // Update all the points in the list having their abscisse in ]x0;x2[
  ae_list_node<point*>* point_node;
  point*  m;
  if (node_0 != NULL) point_node = node_0->get_next();
  else point_node = _points->get_first();

  while (point_node != node_2) {
    m = point_node->get_obj();

    if (point_node == node_1) {
      m->second += y1;
    }
    else {
      if (m->first < x1) {
        m->second += y1 * (m->first - x0) / (x1 - x0);
      }
      else {
        m->second += y1 * (x2 - m->first) / (x2 - x1);
      }
    }

    point_node = point_node->get_next();
  }

  assert(is_increasing());
}

void fuzzy::add(fuzzy* to_add) {
  // Add interpolated points (one point is needed for each abscissa present in either of the two point lists)
  ae_list_node<point*>* point_node = _points->get_first();
  ae_list_node<point*>* point_to_add_node = to_add->get_points()->get_first();
  point*  m;
  point*  point_to_add;

  while (point_to_add_node != NULL) {
    point_to_add = point_to_add_node->get_obj();

    point_node = create_interpolated_point(point_to_add->first, point_node);

    point_to_add_node = point_to_add_node->get_next();
  }


  // Update each point's ordinate
  point_node = _points->get_first();

  while (point_node != NULL) {
    m = point_node->get_obj();

    m->second += to_add->get_y(m->first);  // It is possible to gain some execution time here
                                            // by parsing  to_add's point list and using the
                                            // get_y(double x, ae_list_node* list_entry) function
                                            // instead of get_y(double x)

    point_node = point_node->get_next();
  }

  assert(is_increasing());
}

void fuzzy::sub(fuzzy* to_sub) {
  // 1. Remove multiple discontinuities.
  simplify();
  to_sub->simplify();


  // 2. Add.
  ae_list_node<point*>* bn = _points->get_first();
  ae_list_node<point*>* rn = (to_sub->_points)->get_first();
  ae_list_node<point*>* next_bn = NULL;
  ae_list_node<point*>* next_rn = NULL;
  point* bpt = NULL;
  point* rpt = NULL;
  point* next_bpt = NULL;
  point* next_rpt = NULL;
  point* new_pt = NULL;
  double xb, yb, xr, yr, next_xb, next_yb, next_xr, next_yr;
  double yb_interpol, yr_interpol;
  bool red_discont = false;
  bool black_discont = false;
  ae_list<point*>* result = new ae_list<point*>();

  while (bn != NULL || rn != NULL) { // there are more points to deal with 
    // get data for the black function
    if (bn != NULL) {
      black_discont = false;
      bpt = bn->get_obj();
      xb = bpt->first;
      yb = bpt->second;
      next_bn = bn->get_next();
      if (next_bn != NULL) {
        next_bpt = next_bn->get_obj();
        next_xb = next_bpt->first;
        next_yb = next_bpt->second;

        if (next_xb == xb) {
          black_discont = true;
        }
      }
    }
    else {
      // no more black point, do as if the black function were 0:
      // simply copy the remaining red points (if any) in the result
      // xb and yb have the same values as those of the last black point
      yr_interpol = to_sub->get_y(xb);

      if (rn != NULL && yb != 0) {
        new_pt = new point(xb, -yr_interpol);
        assert(new_pt);
        result->add(new_pt);
      }

      while (rn != NULL) {
        rpt = rn->get_obj();
        xr = rpt->first;
        yr = rpt->second;
        new_pt = new point(xr, -yr);
        assert(new_pt);
        result->add(new_pt);

        rn = rn->get_next();
      }

      break;
    }

    // get data for the red function
    if (rn != NULL) {
      red_discont = false;
      rpt = rn->get_obj();
      xr = rpt->first;
      yr = rpt->second;
      next_rn = rn->get_next();

      if (next_rn != NULL) {
        next_rpt  = next_rn->get_obj();
        next_xr   = next_rpt->first;
        next_yr   = next_rpt->second;

        if (next_xr == xr) {
          red_discont = true;
        }
      }
    }
    else {
      // no more red point, do as if the red function were 0:
      // simply copy the (eventual) remaining black points in the result
      // xr and yr have the same values as those of the last red point
      // printf("no more red point\n");
      yb_interpol = get_y(xr);

      if (bn != NULL && yr != 0) {
        new_pt = new point(xr, yb_interpol);
        assert(new_pt);
        result->add(new_pt);
      }

      while (bn != NULL) {
        bpt = bn->get_obj();
        xb = bpt->first;
        yb = bpt->second;
        new_pt = new point(xb, yb);
        assert(new_pt);
        result->add(new_pt);

        bn = bn->get_next();
      }
      break;
    }


    // compare xb and xr
    if (xb == xr) {
      // ----------------- case 1 : xb == xr -----------------------
      // there will be at at least one point P on x == xb == x_c

      if (! red_discont && ! black_discont) {
        // case 1a : no discontinuity
        // => P(x, yb + yr)
        new_pt = new point(xb, yb - yr);
        assert(new_pt);
        result->add(new_pt);

        bn = next_bn;
        if (bn != NULL) next_bn = bn->get_next();
        rn = next_rn;
        if (rn != NULL) next_rn = rn->get_next();
      }
      else {
        if (black_discont && ! red_discont) {
          // case 1b : discontinuity in the black function only
          // => P-(x, yb- + yr) and P+(x, yb+ + yr)
          // P+ is computed using the following black point
          new_pt = new point(xb, yb - yr);
          assert(new_pt);
          result->add(new_pt);

          new_pt = new point(xb, next_yb - yr);
          assert(new_pt);
          result->add(new_pt);

          if (next_bn != NULL) bn = next_bn->get_next(); // jump over the next black point
          if (bn != NULL) next_bn = bn->get_next();
          rn = next_rn;
          if (rn != NULL) next_rn = rn->get_next();
        }
        else {
          if (! black_discont && red_discont) {
            // case 1c : discontinuity in the red function only
            // => P-(x, yb + yr-) and P+(x, yb + yr+)
            // P+ is computed using the following red point
            new_pt = new point(xb, yb - yr);
            assert(new_pt);
            result->add(new_pt);

            new_pt = new point(xb, yb - next_yr);
            assert(new_pt);
            result->add(new_pt);

            bn = next_bn;
            if (bn != NULL) next_bn = bn->get_next();
            if (next_rn != NULL) rn = next_rn->get_next(); // jump over the next red point
            if (rn != NULL) next_rn = rn->get_next();
          }
          else {
            // case 1d : discontinuity in both functions
            // => P-(x, yb- + yr-) and P+(x, yb+ + yr+)
            // P+ is computed using both the following red and black points
            new_pt = new point(xb, yb - yr);
            assert(new_pt);
            result->add(new_pt);

            new_pt = new point(xb, next_yb - next_yr);
            assert(new_pt);
            result->add(new_pt);

            if (next_bn != NULL) bn = next_bn->get_next(); // jump over the next black point
            if (bn != NULL) next_bn = bn->get_next();
            if (next_rn != NULL) rn = next_rn->get_next(); // jump over the next red point
            if (rn != NULL) next_rn = rn->get_next();
          }
        }
      }
    }
    else {
      // xb is the === zone:
      //
      // ==============|---|---|================
      //           xr-eps  xr  xr+eps

      if(xb < xr) {
        // ----------------- case 2 : xb < xr -----------------------
        // xb is the === zone:
        //
        // ==============|---|---|-----------------
        //           xr-eps  xr  xr+eps
        //
        // there will be at at least one point P on x == xb

        yr_interpol = to_sub->get_y(xb);

        if (! black_discont) {
          // case 2a : no discontinuity
          // => P(xb, yb + red(xb))
          new_pt = new point(xb, yb - yr_interpol);
          assert(new_pt);
          result->add(new_pt);

          bn = next_bn;
          if (bn != NULL) next_bn = bn->get_next();
        }
        else {
          // case 2b : discontinuity on xb in the black function
          // => P-(xb, yb- + red(xb)) and P+(xb, yb+ + red(xb))
          new_pt = new point(xb, yb - yr_interpol);
          assert(new_pt);
          result->add(new_pt);

          new_pt = new point(xb, next_yb - yr_interpol);
          assert(new_pt);
          result->add(new_pt);

          if (next_bn != NULL) bn = next_bn->get_next(); // jump over the next bn
          if (bn != NULL) next_bn = bn->get_next();
        }
        // can there be a discontinuity on xr ? yes, but we will deal
        // with it later, when min{xb,xr} will be xr

        // bn goes one step further and rn doesn't move
        // if there is a discontinuity on xb, the same thing will
        // happen with the following black point
      }
      else {
        // ----------------- case 3 : xb > xr -----------------------
        // xb is the === zone:
        //
        // --------------|---|---|===============
        //           xr-eps  xr  xr+eps
        //
        // there will be at at least one point P on x == xr

        yb_interpol = get_y(xr);

        if (! red_discont) {
          // case 3a : no discontinuity
          // => P(xr, yr + black(xr))
          new_pt = new point(xr, yb_interpol - yr);
          assert(new_pt);
          result->add(new_pt);

          rn = next_rn;
          if (rn != NULL) next_rn = rn->get_next();
        }
        else {
          // case 2b : discontinuity on xr
          // => P-(xr, yr- + black(xr)) and P+(xr, yr+ + black(xr))
          new_pt = new point(xr, yb_interpol - yr);
          assert(new_pt);
          result->add(new_pt);

          new_pt = new point(xr, yb_interpol - next_yr);
          assert(new_pt);
          result->add(new_pt);

          if (next_rn != NULL) rn = next_rn->get_next(); // jump over next red point
          if (rn != NULL) next_rn = rn->get_next();
        }

        // store xr, rn goes one step further and bn doesn't move
        // if there is a discontinuity on xr, the same thing will
        // happen with the following red point
      }
    }
  }

  _points->erase(true);
  delete _points;
  _points = result;

  assert(is_increasing());
}

double fuzzy::get_geometric_area() const {
  double area = 0;
  double tmp, tmp2;

  point*  p =  _points->get_first()->get_obj();
  ae_list_node<point*>* next_point_node = _points->get_first()->get_next();
  point*  next_point = NULL;

  while (next_point_node != NULL) {
    next_point = next_point_node->get_obj();

    if (p->first != next_point->first) {
      if ((p->second > 0 && next_point->second < 0) || (p->second < 0 && next_point->second > 0)) {
        // The segment [p, next_point] crosses the abscissa axis => We will compute the area in 2 parts
        point* pt_zero  = new point(get_x(0.0, p, next_point), 0.0);

        tmp = (p->second + pt_zero->second) * (pt_zero->first - p->first) / 2.0;              // Negative if y1 is negative
        tmp2 =  (pt_zero->second + next_point->second) * (next_point->first - pt_zero->first) / 2.0;  // Negative if y2 is negative

        delete pt_zero;

        area += (tmp  >= 0)? tmp  : -tmp;
        area += (tmp2 >= 0)? tmp2 : -tmp2;
      }
      else {
        // The segment is entirely above or underneath the abscissa axis
        tmp = (p->second + next_point->second) * (next_point->first - p->first) / 2.0;

        area += (tmp >= 0)? tmp : -tmp;
      }
    }

    p = next_point;
    next_point_node = next_point_node->get_next();
  }

  return area;
}

double fuzzy::get_geometric_area(double start_segment, double end_segment) const {
  // Fuzzy set first (resp last) point must be at x = X_MIN (resp x = X_MAX)
  assert(_points->get_first() != _points->get_last());
  assert(_points->get_first()->get_obj()->first == X_MIN);
  assert(_points->get_last()->get_obj()->first == X_MAX);

  // We must have (X_MIN <= start_segment < end_segment <= X_MAX)
  assert(start_segment >= X_MIN && start_segment < end_segment && end_segment <= X_MAX);

  double area = 0;
  double tmp, tmp2;

  ae_list_node<point*>* point_node      = _points->get_first();
  ae_list_node<point*>* next_point_node = point_node->get_next();

  point*  m       = point_node->get_obj();
  point*  next_point  = NULL;

  while (next_point_node != NULL) {
    next_point = next_point_node->get_obj();

    // If there are no points at x = start_segment and x = end_segment, we must interpolate them
    if (m->first < start_segment && next_point->first > start_segment) {
      next_point      = new point(start_segment, get_y(start_segment, point_node));
      next_point_node = _points->add_after(next_point, point_node);

      continue;
    }
    if (m->first < end_segment && next_point->first > end_segment) {
      next_point      = new point(end_segment, get_y(end_segment, point_node));
      next_point_node = _points->add_after(next_point, point_node);

      continue;
    }

    // If we are within the segment, compute the area
    if (m->first >= start_segment && next_point->first <= end_segment) {
      if (m->first != next_point->first) {
        if ((m->second > 0 && next_point->second < 0) || (m->second < 0 && next_point->second > 0)) {
          // The segment [m, next_point] crosses the abscissa axis => We will compute the area in 2 parts
          point* pt_zero  = new point(get_x(0.0, m, next_point), 0.0);

          tmp = (m->second + pt_zero->second) * (pt_zero->first - m->first) / 2.0;              // Negative if y1 is negative
          tmp2 =  (pt_zero->second + next_point->second) * (next_point->first - pt_zero->first) / 2.0;  // Negative if y2 is negative

          delete pt_zero;

          area += (tmp  >= 0)? tmp  : -tmp;
          area += (tmp2 >= 0)? tmp2 : -tmp2;
        }
        else {
          // The segment is entirely above or underneath the abscissa axis
          tmp = (m->second + next_point->second) * (next_point->first - m->first) / 2.0;

          area += (tmp >= 0)? tmp : -tmp;
        }
      }
    }

    m       = next_point;
    point_node  = next_point_node;
    next_point_node = next_point_node->get_next();
  }

  return area;
}


void fuzzy::add_upper_bound(double upper_bound) {
  ae_list_node<point*>* point_node = _points->get_first();
  point*  m;

  ae_list_node<point*>* prev_node;
  ae_list_node<point*>* next_node;
  point*  prev_point;
  point*  next_point;

  while (point_node != NULL) {
    m = point_node->get_obj();

    if (m->second > upper_bound) {
      prev_node = point_node->get_prev();
      next_node = point_node->get_next();

      if (prev_node != NULL) {
        prev_point = prev_node->get_obj();

        if (prev_point->second < upper_bound) { // In fact it can only be < or == 
          // Create a point at the intersection of the segment and the ceiling defined by the upper bound
          _points->add_after(new point(get_x(upper_bound, prev_point, m), upper_bound), prev_node);
        }
      }

      if (next_node != NULL) {
        next_point = next_node->get_obj();

        if (next_point->second < upper_bound) {
          // Create a point at the intersection of the segment and the ceiling defined by the upper bound
          _points->add_after(new point(get_x(upper_bound, m, next_point), upper_bound), point_node);
        }
      }

      // Lower current point down to upper_bound
      m->second = upper_bound;
    }

    point_node = point_node->get_next();
  }

  assert(is_increasing());
}

void fuzzy::add_lower_bound(double lower_bound) {
  ae_list_node<point*>* point_node = _points->get_first();
  point*  m;

  ae_list_node<point*>* prev_node;
  ae_list_node<point*>* next_node;
  point*  prev_point;
  point*  next_point;

  while (point_node != NULL) {
    m = point_node->get_obj();

    if (m->second < lower_bound) {
      prev_node = point_node->get_prev();
      next_node = point_node->get_next();

      if (prev_node != NULL) {
        prev_point = prev_node->get_obj();

        if (prev_point->second > lower_bound) { // In fact it can only be > or ==
          // Create a point at the intersection of the segment and the floor defined by the lower bound
          _points->add_after(new point(get_x(lower_bound, prev_point, m), lower_bound), prev_node);

	  assert(is_increasing());
        }
        // else nothing to do
      }

      if (next_node != NULL) {
        next_point = next_node->get_obj();

        if (next_point->second > lower_bound) {
          // Create a point at the intersection of the segment and the floor defined by the lower bound
          _points->add_after(new point(get_x(lower_bound, m, next_point), lower_bound), point_node);

	  assert(is_increasing());
        }
        // else nothing to do
      }

      // Raise current point up to lower_bound
      m->second = lower_bound;
    }

    point_node = point_node->get_next();
  }

  assert(is_increasing());
}


bool fuzzy::is_identical_to(const fuzzy * other, double tolerance ) const {
  if (_points->get_nb_elts() != other->_points->get_nb_elts()) {
    return false;
  }

  ae_list_node<point*>* point_node = _points->get_first();
  point*  m = NULL;

  ae_list_node<point*>* other_node = other->_points->get_first();
  point*  other_point = NULL;

  bool ok = true;

  while (ok && (point_node != NULL)) {
    m = point_node->get_obj();
    other_point = other_node->get_obj();

    if (fabs(m->first - other_point->first)> tolerance*(fabs(m->first) + fabs(other_point->first)) || fabs(m->second - other_point->second)>tolerance*(fabs(m->second) + fabs(other_point->second))) {ok = false;}

    point_node = point_node->get_next();
    other_node = other_node->get_next();
  }

  return ok;
}


void fuzzy::save(gzFile backup_file) const {
  int16_t nb_points = (_points == NULL) ? 0 : _points->get_nb_elts();
  gzwrite(backup_file, &nb_points, sizeof(nb_points));

  if (_points != NULL) {
    ae_list_node<point*>* point_node = _points->get_first();
    point*  m;
    for (int16_t i = 0 ; i < nb_points ; i++) {
      m = point_node->get_obj();

      writepoint(*m, backup_file);

      point_node = point_node->get_next();
    }
  }
}


void fuzzy::load(gzFile backup_file) {
  int16_t nb_points;
  gzread(backup_file, &nb_points, sizeof(nb_points));

  if (nb_points > 0) _points = new ae_list<point*>();

  for (int16_t i = 0 ; i < nb_points ; i++) {
    _points->add(new point(readpoint(backup_file)));
  }
}


// =================================================================
//                           Protected Methods
// =================================================================
double fuzzy::get_y(double x, ae_list_node<point*>* list_entry) const {
  if (list_entry == NULL) list_entry = _points->get_first();

  assert(x >= X_MIN && x <= X_MAX);
  assert(list_entry != NULL);
  assert(list_entry->get_obj()->first <= x);
  assert(_points->get_nb_elts() >= 2);

  // Look for the first point of the list having m->first >= x
  ae_list_node<point*>* point_node  = list_entry;
  point*  m       = point_node->get_obj();

  while (m->first < x) {
    point_node = point_node->get_next();
    m      = point_node->get_obj();
  }

  // If a point with abscissa x exists, return it's y, otherwise compute the needed y by interpolation
  // from the 2 flanking points
  if (m->first == x) return m->second;
  else {
    point* prev_point = point_node->get_prev()->get_obj();
    return prev_point->second + (m->second - prev_point->second) * (x - prev_point->first) / (m->first - prev_point->first);
  }
}

ae_list_node<point*>* fuzzy::create_interpolated_point(double x, ae_list_node<point*>* list_entry) {
  if (list_entry == NULL) list_entry = _points->get_first();

  assert(x >= X_MIN && x <= X_MAX);
  assert(list_entry->get_obj()->first <= x);
  assert(_points->get_nb_elts() >= 2);

  // Look for the first point of the list having m->first >= x
  ae_list_node<point*>* point_node  = list_entry;
  point*  m       = point_node->get_obj();

  while (m->first < x) {
    point_node = point_node->get_next();
    m      = point_node->get_obj();
  }

  if (m->first == x) {
    // The point already exists, nothing to do
    return point_node;
  }
  else {
    point_node = _points->add_before(new point(x, get_y(x)), point_node);
    return point_node;
  }

  assert(is_increasing());
}

bool fuzzy::is_increasing() const {
  assert((_points->get_first()->get_obj())->first == X_MIN);
  assert((_points->get_last()->get_obj())->first == X_MAX);

  ae_list_node<point*>* point_node  = _points->get_first();
  ae_list_node<point*>* next_point_node = point_node->get_next();
  point*  m;
  point*  next_point;

  while (next_point_node != NULL) {
    m = point_node->get_obj();
    next_point = next_point_node->get_obj();

    if (m->first > next_point->first)
      return false;

    point_node = point_node->get_next();
    next_point_node = next_point_node->get_next();
  }
  return true;
}  

} // namespace aevol
