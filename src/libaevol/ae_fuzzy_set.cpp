//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, Virginie Lefort
//                    David Parsons
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************


/** \class
 *  \brief
 */


// =================================================================
//                              Libraries
// =================================================================
#include <math.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_fuzzy_set.h>




//##############################################################################
//                                                                             #
//                             Class ae_fuzzy_set                              #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_fuzzy_set::ae_fuzzy_set( void )
{
  _points = new ae_list();
  
  initialize();
}

ae_fuzzy_set::ae_fuzzy_set( const ae_fuzzy_set &model )
{
  _points = new ae_list();

  ae_list_node* model_point_parser = model._points->get_first();
  ae_point_2d*  model_point;

  while ( model_point_parser != NULL )
  {
    model_point = (ae_point_2d*) model_point_parser->get_obj();

    _points->add( new ae_point_2d( *model_point ) );

    model_point_parser = model_point_parser->get_next();
  }
}

ae_fuzzy_set::ae_fuzzy_set( gzFile backup_file )
{
  _points = new ae_list();
  
  load( backup_file );
}

// =================================================================
//                             Destructors
// =================================================================
ae_fuzzy_set::~ae_fuzzy_set( void )
{
  if ( _points != NULL )
  {
    _points->erase( DELETE_OBJ );
    delete _points;
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_fuzzy_set::simplify( void )
{
  // Removes useless points like A or B:
  //
  //  ^      B
  //  |   +--+-----+
  //  |   |         \     ---
  //  |  A+          \   /
  //  |   |           \ /
  //  |---+------------+-------+----->

  //  printf("initial nb of points %d\n",_points.get_count());

  ae_list_node * node = _points->get_first();
  ae_list_node * next_node = NULL;
  ae_list_node * other_node = NULL;
  ae_point_2d * pt = NULL;
  ae_point_2d * next_pt = NULL;
  ae_point_2d * other_pt = NULL;
  double x, y, next_x, next_y, other_x, other_y;

  while ( node != NULL )
  {
    pt = (ae_point_2d *) node->get_obj();
    x = pt->x;
    y = pt->y;
    next_node = node->get_next();
    if ( next_node != NULL )
    {
      next_pt = (ae_point_2d *) next_node->get_obj();
      next_x = next_pt->x;
      next_y = next_pt->y;
      if ( next_x == x )
      {
        // discontinuity on x
        // is the next point also on the same x ?
        other_node = next_node->get_next();
        while ( other_node != NULL )
        {
          other_pt = (ae_point_2d *) other_node->get_obj();
          other_x = other_pt->x;
          other_y = other_pt->y;
          if ( other_x == x )
          {
            _points->remove( next_node, DELETE_OBJ, DELETE_OBJ );
            next_node = other_node;
          }
          else
          {
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

  while ( node != NULL )
  {
    pt = (ae_point_2d *) node->get_obj();
    x = pt->x;
    y = pt->y;
    next_node = node->get_next();
    if ( next_node != NULL )
    {
      next_pt = (ae_point_2d *) next_node->get_obj();
      next_x = next_pt->x;
      next_y = next_pt->y;
      if ( next_y == y )
      {
        // plateau on y
        // is the next point also on the same y ?
        other_node = next_node->get_next();
        while ( other_node != NULL )
        {
          other_pt = (ae_point_2d *) other_node->get_obj();
          other_x = other_pt->x;
          other_y = other_pt->y;
          if ( other_y == y )
          {
            _points->remove( next_node, DELETE_OBJ, DELETE_OBJ );
            next_node = other_node;
            next_pt = (ae_point_2d *) next_node->get_obj();
          }
          else break;
          other_node = other_node->get_next();
        }
      }
    }
    node = next_node;
  }
  
  #ifdef DEBUG
    _assert_order();
  #endif
}

void ae_fuzzy_set::print_points( void ) const
{
  ae_list_node* point_node = _points->get_first();
  ae_point_2d*  point;

  while ( point_node != NULL )
  {
    point = (ae_point_2d*) point_node->get_obj();

    printf( "  (%f, %f)\n", point->x, point->y );

    point_node = point_node->get_next();
  }
}

void ae_fuzzy_set::add_triangle( double mean, double width, double height )
{
  assert( X_MIN <= mean   && mean <= X_MAX );
  assert( W_MIN <= width ); // the maximum width depends on each individual
  // assert( MIN_H <= height && height <= MAX_H ); Not necessarily because the concentration can be > 1
  
  if ( fabs(width) < 1e-15 || fabs(height) < 1e-15 ) return;
  
  // Compute triangle points' coordinates
  double x0 = mean - width;
  double x1 = mean;
  double x2 = mean + width;
  //~ double y0 = 0.0;
  double y1 = height;
  //~ double y2 = 0.0;

  // We will need at least a point at abscissa x0, x1 and x2,
  // Create them if they don't exist or update them if they do.
  ae_list_node* node_0 = NULL;
  ae_list_node* node_1 = NULL;
  ae_list_node* node_2 = NULL;

  if ( x0 >= X_MIN )  node_0 = create_interpolated_point( x0 );
                      node_1 = create_interpolated_point( x1, node_0 );
  if ( x2 <= X_MAX )  node_2 = create_interpolated_point( x2, node_1 );


  // Update all the points in the list having their abscisse in ]x0;x2[
  ae_list_node* point_node;
  ae_point_2d*  point;
  if ( node_0 != NULL ) point_node = node_0->get_next();
  else point_node = _points->get_first();

  while ( point_node != node_2 )
  {
    point = (ae_point_2d*)point_node->get_obj();

    if ( point_node == node_1 )
    {
      point->y += y1;
    }
    else
    {
      if ( point->x < x1 )
      {
        point->y += y1 * (point->x - x0) / (x1 - x0);
      }
      else
      {
        point->y += y1 * (x2 - point->x) / (x2 - x1);
      }
    }

    point_node = point_node->get_next();
  }
  
  #ifdef DEBUG
    _assert_order();
  #endif
}

void ae_fuzzy_set::add( ae_fuzzy_set* to_add )
{
  // Add interpolated points (one point is needed for each abscissa present in either of the two point lists)
  ae_list_node* point_node = _points->get_first();
  ae_list_node* point_to_add_node = to_add->get_points()->get_first();
  ae_point_2d*  point;
  ae_point_2d*  point_to_add;

  while ( point_to_add_node != NULL )
  {
    point_to_add = (ae_point_2d*) point_to_add_node->get_obj();

    point_node = create_interpolated_point( point_to_add->x, point_node );

    point_to_add_node = point_to_add_node->get_next();
  }


  // Update each point's ordinate
  point_node = _points->get_first();

  while ( point_node != NULL )
  {
    point = (ae_point_2d*) point_node->get_obj();

    point->y += to_add->get_y( point->x );  // It is possible to gain some execution time here
                                            // by parsing  to_add's point list and using the
                                            // get_y( double x, ae_list_node* list_entry ) function
                                            // instead of get_y( double x )

    point_node = point_node->get_next();
  }
  
  #ifdef DEBUG
    _assert_order();
  #endif
}

void ae_fuzzy_set::sub( ae_fuzzy_set* to_sub )
{
  // 1. Remove multiple discontinuities.
  simplify();
  to_sub->simplify();


  // 2. Add.
  ae_list_node* bn = _points->get_first();
  ae_list_node* rn = (to_sub->_points)->get_first();
  ae_list_node* next_bn = NULL;
  ae_list_node* next_rn = NULL;
  ae_point_2d* bpt = NULL;
  ae_point_2d* rpt = NULL;
  ae_point_2d* next_bpt = NULL;
  ae_point_2d* next_rpt = NULL;
  ae_point_2d* new_pt = NULL;
  double xb, yb, xr, yr, next_xb, next_yb, next_xr, next_yr;
  double yb_interpol, yr_interpol;
  bool red_discont = false;
  bool black_discont = false;
  ae_list* result = new ae_list();

  while ( bn != NULL || rn != NULL ) // there are more points to deal with
  {
    // get data for the black function
    if ( bn != NULL )
    {
      black_discont = false;
      bpt = (ae_point_2d *) bn->get_obj();
      xb = bpt->x;
      yb = bpt->y;
      next_bn = bn->get_next();
      if ( next_bn != NULL )
      {
        next_bpt = (ae_point_2d *) next_bn->get_obj();
        next_xb = next_bpt->x;
        next_yb = next_bpt->y;

        if ( next_xb == xb )
        {
          black_discont = true;
        }
      }
    }
    else
    {
      // no more black point, do as if the black function were 0:
      // simply copy the remaining red points (if any) in the result
      // xb and yb have the same values as those of the last black point
      yr_interpol = to_sub->get_y( xb );

      if ( rn != NULL && yb != 0 )
      {
        new_pt = new ae_point_2d( xb, -yr_interpol );
        assert( new_pt );
        result->add( new_pt );
      }

      while ( rn != NULL )
      {
        rpt = (ae_point_2d *) rn->get_obj();
        xr = rpt->x;
        yr = rpt->y;
        new_pt = new ae_point_2d( xr, -yr );
        assert( new_pt );
        result->add( new_pt );

        rn = rn->get_next();
      }

      break;
    }

    // get data for the red function
    if ( rn != NULL )
    {
      red_discont = false;
      rpt = (ae_point_2d *) rn->get_obj();
      xr = rpt->x;
      yr = rpt->y;
      next_rn = rn->get_next();

      if ( next_rn != NULL )
      {
        next_rpt = (ae_point_2d *) next_rn->get_obj();
        next_xr = next_rpt->x;
        next_yr = next_rpt->y;

        if ( next_xr == xr )
        {
          red_discont = true;
        }
      }
    }
    else
    {
      // no more red point, do as if the red function were 0:
      // simply copy the (eventual) remaining black points in the result
      // xr and yr have the same values as those of the last red point
      // printf("no more red point\n");
      yb_interpol = get_y( xr );

      if ( bn != NULL && yr != 0 )
      {
        new_pt = new ae_point_2d( xr, yb_interpol );
        assert( new_pt );
        result->add( new_pt );
      }

      while ( bn != NULL )
      {
        bpt = (ae_point_2d *) bn->get_obj();
        xb = bpt->x;
        yb = bpt->y;
        new_pt = new ae_point_2d( xb, yb );
        assert( new_pt );
        result->add( new_pt );

        bn = bn->get_next();
      }
      break;
    }


    // compare xb and xr
    if ( xb == xr )
    {
      // ----------------- case 1 : xb == xr -----------------------
      // there will be at at least one point P on x == xb == x_c

      if ( ! red_discont && ! black_discont )
      {
        // case 1a : no discontinuity
        // => P(x, yb + yr)
        new_pt = new ae_point_2d( xb, yb - yr );
        assert( new_pt );
        result->add( new_pt );

        bn = next_bn;
        if ( bn != NULL ) next_bn = bn->get_next();
        rn = next_rn;
        if ( rn != NULL ) next_rn = rn->get_next();
      }
      else
      {
        if ( black_discont && ! red_discont )
        {
          // case 1b : discontinuity in the black function only
          // => P-(x, yb- + yr) and P+(x, yb+ + yr)
          // P+ is computed using the following black point
          new_pt = new ae_point_2d( xb, yb - yr );
          assert( new_pt );
          result->add( new_pt );

          new_pt = new ae_point_2d( xb, next_yb - yr );
          assert( new_pt );
          result->add( new_pt );

          if ( next_bn != NULL ) bn = next_bn->get_next(); // jump over the next black point
          if ( bn != NULL ) next_bn = bn->get_next();
          rn = next_rn;
          if ( rn != NULL ) next_rn = rn->get_next();
        }
        else
        {
          if ( ! black_discont && red_discont )
          {
            // case 1c : discontinuity in the red function only
            // => P-(x, yb + yr-) and P+(x, yb + yr+)
            // P+ is computed using the following red point
            new_pt = new ae_point_2d( xb, yb - yr );
            assert( new_pt );
            result->add( new_pt );

            new_pt = new ae_point_2d( xb, yb - next_yr );
            assert( new_pt );
            result->add( new_pt );

            bn = next_bn;
            if ( bn != NULL ) next_bn = bn->get_next();
            if ( next_rn != NULL ) rn = next_rn->get_next(); // jump over the next red point
            if ( rn != NULL ) next_rn = rn->get_next();
          }
          else
          {
            // case 1d : discontinuity in both functions
            // => P-(x, yb- + yr-) and P+(x, yb+ + yr+)
            // P+ is computed using both the following red and black points
            new_pt = new ae_point_2d( xb, yb - yr );
            assert( new_pt );
            result->add( new_pt );

            new_pt = new ae_point_2d( xb, next_yb - next_yr );
            assert( new_pt );
            result->add( new_pt );

            if ( next_bn != NULL ) bn = next_bn->get_next(); // jump over the next black point
            if ( bn != NULL ) next_bn = bn->get_next();
            if ( next_rn != NULL ) rn = next_rn->get_next(); // jump over the next red point
            if ( rn != NULL ) next_rn = rn->get_next();
          }
        }
      }
    }
    else
    {
      // xb is the === zone:
      //
      // ==============|---|---|================
      //           xr-eps  xr  xr+eps

      if( xb < xr )
      {
        // ----------------- case 2 : xb < xr -----------------------
        // xb is the === zone:
        //
        // ==============|---|---|-----------------
        //           xr-eps  xr  xr+eps
        //
        // there will be at at least one point P on x == xb

        yr_interpol = to_sub->get_y( xb );

        if ( ! black_discont )
        {
          // case 2a : no discontinuity
          // => P(xb, yb + red(xb))
          new_pt = new ae_point_2d( xb, yb - yr_interpol );
          assert( new_pt );
          result->add( new_pt );

          bn = next_bn;
          if ( bn != NULL ) next_bn = bn->get_next();
        }
        else
        {
          // case 2b : discontinuity on xb in the black function
          // => P-(xb, yb- + red(xb)) and P+(xb, yb+ + red(xb))
          new_pt = new ae_point_2d( xb, yb - yr_interpol );
          assert( new_pt );
          result->add( new_pt );

          new_pt = new ae_point_2d( xb, next_yb - yr_interpol );
          assert( new_pt );
          result->add( new_pt );

          if ( next_bn != NULL ) bn = next_bn->get_next(); // jump over the next bn
          if ( bn != NULL ) next_bn = bn->get_next();
        }
        // can there be a discontinuity on xr ? yes, but we will deal
        // with it later, when min{xb,xr} will be xr

        // bn goes one step further and rn doesn't move
        // if there is a discontinuity on xb, the same thing will
        // happen with the following black point
      }
      else
      {
        // ----------------- case 3 : xb > xr -----------------------
        // xb is the === zone:
        //
        // --------------|---|---|===============
        //           xr-eps  xr  xr+eps
        //
        // there will be at at least one point P on x == xr

        yb_interpol = get_y( xr );

        if ( ! red_discont )
        {
          // case 3a : no discontinuity
          // => P(xr, yr + black(xr))
          new_pt = new ae_point_2d( xr, yb_interpol - yr );
          assert( new_pt );
          result->add( new_pt );

          rn = next_rn;
          if ( rn != NULL ) next_rn = rn->get_next();
        }
        else
        {
          // case 2b : discontinuity on xr
          // => P-(xr, yr- + black(xr)) and P+(xr, yr+ + black(xr))
          new_pt = new ae_point_2d( xr, yb_interpol - yr );
          assert( new_pt );
          result->add( new_pt );

          new_pt = new ae_point_2d( xr, yb_interpol - next_yr );
          assert( new_pt );
          result->add( new_pt );

          if ( next_rn != NULL ) rn = next_rn->get_next(); // jump over next red point
          if ( rn != NULL ) next_rn = rn->get_next();
        }

        // store xr, rn goes one step further and bn doesn't move
        // if there is a discontinuity on xr, the same thing will
        // happen with the following red point
      }
    }
  }

  _points->erase( DELETE_OBJ );
  delete _points;
  _points = result;
  
  #ifdef DEBUG
    _assert_order();
  #endif
}

double ae_fuzzy_set::get_geometric_area( void ) const
{
  double area = 0;
  double tmp, tmp2;
  
  ae_point_2d*  point           = (ae_point_2d*) _points->get_first()->get_obj();
  ae_list_node* next_point_node = _points->get_first()->get_next();
  ae_point_2d*  next_point      = NULL;

  while ( next_point_node != NULL )
  {
    next_point = (ae_point_2d*) next_point_node->get_obj();

    if ( point->x != next_point->x )
    {
      if ( (point->y > 0 && next_point->y < 0) || (point->y < 0 && next_point->y > 0) )
      {
        // The segment [point, next_point] crosses the abscissa axis => We will compute the area in 2 parts
        ae_point_2d* pt_zero  = new ae_point_2d( get_x( 0.0, point, next_point ), 0.0 );

        tmp = (point->y + pt_zero->y) * (pt_zero->x - point->x) / 2.0;              // Negative if y1 is negative
        tmp2 =  (pt_zero->y + next_point->y) * (next_point->x - pt_zero->x) / 2.0;  // Negative if y2 is negative

        delete pt_zero;

        area += ( tmp  >= 0 )? tmp  : -tmp;
        area += ( tmp2 >= 0 )? tmp2 : -tmp2;
      }
      else
      {
        // The segment is entirely above or underneath the abscissa axis
        tmp = (point->y + next_point->y) * (next_point->x - point->x) / 2.0;

        area += ( tmp >= 0 )? tmp : -tmp;
      }
    }

    point = next_point;
    next_point_node = next_point_node->get_next();
  }

  return area;
}

double ae_fuzzy_set::get_geometric_area( double start_segment, double end_segment ) const
{
  // Fuzzy set first (resp last) point must be at x = X_MIN (resp x = X_MAX)
  assert( _points->get_first() != _points->get_last() );
  assert( ((ae_point_2d*)_points->get_first()->get_obj())->x == X_MIN );
  assert( ((ae_point_2d*)_points->get_last()->get_obj())->x == X_MAX );
  
  // We must have ( X_MIN <= start_segment < end_segment <= X_MAX )
  assert( start_segment >= X_MIN && start_segment < end_segment && end_segment <= X_MAX );
  
  double area = 0;
  double tmp, tmp2;
  
  ae_list_node* point_node      = _points->get_first();
  ae_list_node* next_point_node = point_node->get_next();
  
  ae_point_2d*  point       = (ae_point_2d*) point_node->get_obj();
  ae_point_2d*  next_point  = NULL;

  while ( next_point_node != NULL )
  {
    next_point = (ae_point_2d*) next_point_node->get_obj();
    
    // If there are no points at x = start_segment and x = end_segment, we must interpolate them
    if ( point->x < start_segment && next_point->x > start_segment )
    {
      next_point      = new ae_point_2d( start_segment, get_y( start_segment, point_node ) );
      next_point_node = _points->add_after( next_point, point_node );
      
      continue;
    }
    if ( point->x < end_segment && next_point->x > end_segment )
    {
      next_point      = new ae_point_2d( end_segment, get_y( end_segment, point_node ) );
      next_point_node = _points->add_after( next_point, point_node );
      
      continue;
    }

    // If we are within the segment, compute the area
    if ( point->x >= start_segment && next_point->x <= end_segment ) 
    {
      if ( point->x != next_point->x )
      {
        if ( (point->y > 0 && next_point->y < 0) || (point->y < 0 && next_point->y > 0) )
        {
          // The segment [point, next_point] crosses the abscissa axis => We will compute the area in 2 parts
          ae_point_2d* pt_zero  = new ae_point_2d( get_x( 0.0, point, next_point ), 0.0 );

          tmp = (point->y + pt_zero->y) * (pt_zero->x - point->x) / 2.0;              // Negative if y1 is negative
          tmp2 =  (pt_zero->y + next_point->y) * (next_point->x - pt_zero->x) / 2.0;  // Negative if y2 is negative

          delete pt_zero;

          area += ( tmp  >= 0 )? tmp  : -tmp;
          area += ( tmp2 >= 0 )? tmp2 : -tmp2;
        }
        else
        {
          // The segment is entirely above or underneath the abscissa axis
          tmp = (point->y + next_point->y) * (next_point->x - point->x) / 2.0;

          area += ( tmp >= 0 )? tmp : -tmp;
        }
      }
    }
    
    point       = next_point;
    point_node  = next_point_node;
    next_point_node = next_point_node->get_next();
  }
  
  return area;
}


void ae_fuzzy_set::add_upper_bound( double upper_bound )
{
  ae_list_node* point_node = _points->get_first();
  ae_point_2d*  point;

  ae_list_node* prev_node;
  ae_list_node* next_node;
  ae_point_2d*  prev_point;
  ae_point_2d*  next_point;

  while ( point_node != NULL )
  {
    point = (ae_point_2d*)point_node->get_obj();

    if ( point->y > upper_bound )
    {
      prev_node = point_node->get_prev();
      next_node = point_node->get_next();

      if ( prev_node != NULL )
      {
        prev_point = (ae_point_2d*)prev_node->get_obj();

        if ( prev_point->y < upper_bound ) // In fact it can only be < or ==
        {
          // Create a point at the intersection of the segment and the ceiling defined by the upper bound
          _points->add_after( new ae_point_2d( get_x( upper_bound, prev_point, point ), upper_bound ), prev_node );
        }
      }

      if ( next_node != NULL )
      {
        next_point = (ae_point_2d*)next_node->get_obj();

        if ( next_point->y < upper_bound )
        {
          // Create a point at the intersection of the segment and the ceiling defined by the upper bound
          _points->add_after( new ae_point_2d( get_x( upper_bound, point, next_point ), upper_bound ), point_node );
        }
      }

      // Lower current point down to upper_bound
      point->y = upper_bound;
    }

    point_node = point_node->get_next();
  }
  
  #ifdef DEBUG
    _assert_order();
  #endif
}

void ae_fuzzy_set::add_lower_bound( double lower_bound )
{
  ae_list_node* point_node = _points->get_first();
  ae_point_2d*  point;

  ae_list_node* prev_node;
  ae_list_node* next_node;
  ae_point_2d*  prev_point;
  ae_point_2d*  next_point;

  while ( point_node != NULL )
  {
    point = (ae_point_2d*)point_node->get_obj();

    if ( point->y < lower_bound )
    {
      prev_node = point_node->get_prev();
      next_node = point_node->get_next();

      if ( prev_node != NULL )
      {
        prev_point = (ae_point_2d*)prev_node->get_obj();

        if ( prev_point->y > lower_bound ) // In fact it can only be > or ==
        {
          // Create a point at the intersection of the segment and the floor defined by the lower bound
          _points->add_after( new ae_point_2d( get_x( lower_bound, prev_point, point ), lower_bound ), prev_node );
            
          #ifdef DEBUG
            _assert_order();
          #endif
        }
        // else nothing to do
      }

      if ( next_node != NULL )
      {
        next_point = (ae_point_2d*)next_node->get_obj();

        if ( next_point->y > lower_bound )
        {
          // Create a point at the intersection of the segment and the floor defined by the lower bound
          _points->add_after( new ae_point_2d( get_x( lower_bound, point, next_point ), lower_bound ), point_node );
            
          #ifdef DEBUG
            _assert_order();
          #endif
        }
        // else nothing to do
      }

      // Raise current point up to lower_bound
      point->y = lower_bound;
    }

    point_node = point_node->get_next();
  }
  
  #ifdef DEBUG
    _assert_order();
  #endif
}


bool ae_fuzzy_set::is_identical_to( const ae_fuzzy_set * other) const
{
  if ( _points->get_nb_elts() != other->_points->get_nb_elts())
  {
    return false;
  }

  ae_list_node* point_node = _points->get_first();
  ae_point_2d*  point = NULL;

  ae_list_node* other_node = other->_points->get_first();
  ae_point_2d*  other_point = NULL;

  bool ok = true;

  while ( ok && (point_node != NULL) )
  {
    point = (ae_point_2d*)point_node->get_obj();
    other_point = (ae_point_2d*)other_node->get_obj();

    if ((point->x != other_point->x) || (point->y != other_point->y))
      {ok = false;}

    point_node = point_node->get_next();
    other_node = other_node->get_next();
  }
  
  return ok;
}


void ae_fuzzy_set::save( gzFile backup_file ) const
{
  int16_t nb_points = (_points == NULL) ? 0 : _points->get_nb_elts();
  gzwrite( backup_file, &nb_points, sizeof(nb_points) );
  
  if ( _points != NULL )
  {
    ae_list_node* point_node = _points->get_first();
    ae_point_2d*  point;
    for ( int16_t i = 0 ; i < nb_points ; i++ )
    {
      point = ( ae_point_2d* ) point_node->get_obj();
      
      point->save( backup_file );
      
      point_node = point_node->get_next();
    }
  }
}


void ae_fuzzy_set::load( gzFile backup_file )
{
  int16_t nb_points;
  gzread( backup_file, &nb_points, sizeof(nb_points) );
  
  if ( nb_points > 0 ) _points = new ae_list();
  
  for ( int16_t i = 0 ; i < nb_points ; i++ )
  {
    _points->add( new ae_point_2d( backup_file ) );
  }
}


// =================================================================
//                           Protected Methods
// =================================================================
double ae_fuzzy_set::get_y( double x, ae_list_node* list_entry ) const
{
  if ( list_entry == NULL ) list_entry = _points->get_first();

  assert( x >= X_MIN && x <= X_MAX );
  assert( list_entry != NULL );
  assert( ((ae_point_2d*)list_entry->get_obj())->x <= x );
  assert( _points->get_nb_elts() >= 2 );

  // Look for the first point of the list having point->x >= x
  ae_list_node* point_node  = list_entry;
  ae_point_2d*  point       = (ae_point_2d*) point_node->get_obj();

  while ( point->x < x )
  {
    point_node = point_node->get_next();
    point      = (ae_point_2d*) point_node->get_obj();
  }

  // If a point with abscissa x exists, return it's y, otherwise compute the needed y by interpolation
  // from the 2 flanking points
  if ( point->x == x ) return point->y;
  else
  {
    ae_point_2d* prev_point = (ae_point_2d*) point_node->get_prev()->get_obj();
    return prev_point->y + (point->y - prev_point->y) * (x - prev_point->x) / (point->x - prev_point->x);
  }
}

ae_list_node* ae_fuzzy_set::create_interpolated_point( double x, ae_list_node* list_entry )
{
  if ( list_entry == NULL ) list_entry = _points->get_first();

  assert( x >= X_MIN && x <= X_MAX );
  assert( ((ae_point_2d*)list_entry->get_obj())->x <= x );
  assert( _points->get_nb_elts() >= 2 );

  // Look for the first point of the list having point->x >= x
  ae_list_node* point_node  = list_entry;
  ae_point_2d*  point       = (ae_point_2d*) point_node->get_obj();

  while ( point->x < x )
  {
    point_node = point_node->get_next();
    point      = (ae_point_2d*) point_node->get_obj();
  }

  if ( point->x == x )
  {
    // The point already exists, nothing to do
    return point_node;
  }
  else
  {
    point_node = _points->add_before( new ae_point_2d( x, get_y( x ) ), point_node );
    return point_node;
  }
  
  #ifdef DEBUG
    _assert_order();
  #endif
}

void ae_fuzzy_set::_assert_order( void )
{
  assert( ((ae_point_2d*)_points->get_first()->get_obj())->x == X_MIN );
  assert( ((ae_point_2d*)_points->get_last()->get_obj())->x == X_MAX );
  
  ae_list_node* point_node  = _points->get_first();
  ae_list_node* next_point_node = point_node->get_next();
  ae_point_2d*  point;
  ae_point_2d*  next_point;

  while ( next_point_node != NULL )
  {
    point = (ae_point_2d*) point_node->get_obj();
    next_point = (ae_point_2d*) next_point_node->get_obj();
    
    assert( point->x <= next_point->x );
    
    point_node = point_node->get_next();
    next_point_node = next_point_node->get_next();
  }
}

