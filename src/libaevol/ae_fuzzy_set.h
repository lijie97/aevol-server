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


/** \class
 *  \brief This class provides management tools for "fuzzy sets"
 *         abscissa are bound between MIN_X and MAX_X (defined elsewhere)
 *         A "fuzzy set" should always have at least two points of abscissa
 *         MIN_X and MAX_X
 */


#ifndef __AE_FUZZY_SET_H__
#define __AE_FUZZY_SET_H__


#include <fuzzy.h>
typedef aevol::fuzzy ae_fuzzy_set;
#if 0
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ae_list.h>
#include <ae_point_2d.h>
#include <ae_macros.h>




// =================================================================
//                          Class declarations
// =================================================================
class exp_manager;






class ae_fuzzy_set
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_fuzzy_set( void );
    ae_fuzzy_set( const ae_fuzzy_set &model );
    ae_fuzzy_set( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_fuzzy_set( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_list<ae_point_2d*>* get_points( void );

    // =================================================================
    //                            Public Methods
    // =================================================================
    inline void initialize( void );
    void simplify( void );
    void print_points( void ) const;
    void add_triangle( double mean, double width, double height );
    void add( ae_fuzzy_set* to_add );
    void sub( ae_fuzzy_set* to_sub );
    double get_geometric_area( void ) const;
    double get_geometric_area( double start_segment, double end_segment ) const;
    inline double get_y( double x ) const;
    inline double get_x( double y, ae_point_2d* left_point, ae_point_2d* right_point ) const;
    inline ae_list_node<ae_point_2d*>* create_interpolated_point( double x );

    void add_upper_bound( double upper_bound );
    void add_lower_bound( double lower_bound );

    bool is_identical_to( const ae_fuzzy_set * other, double tolerance) const;

    void save( gzFile backup_file ) const;
    void load( gzFile backup_file );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*
    ae_fuzzy_set( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    */
    /*
    ae_fuzzy_set( const ae_fuzzy_set &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    */

    // =================================================================
    //                           Protected Methods
    // =================================================================
    void _assert_order( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    exp_manager*            _exp_m;
    ae_list<ae_point_2d*>*  _points;





  private :

    // WARNING :  The 2 following functions are private because they are dangerous.
    //            They are optimized versions of their omonyms, beginning their search at the given list entry.
    //            The following conditions are MANDATORY :
    //              list_entry must be a node of the point list (_points)
    //              the corresponding point's abscissa must be <= x
    double get_y( double x, ae_list_node<ae_point_2d*>* list_entry ) const;
    ae_list_node<ae_point_2d*>* create_interpolated_point( double x, ae_list_node<ae_point_2d*>* list_entry );
};


// =====================================================================
//                          Accessors definitions
// =====================================================================
ae_list<ae_point_2d*>* ae_fuzzy_set::get_points( void )
{
  return _points;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
void ae_fuzzy_set::initialize( void )
{
  _points->erase( true );
  _points->add( new ae_point_2d( X_MIN, 0.0 ) );
  _points->add( new ae_point_2d( X_MAX, 0.0 ) );
}

double ae_fuzzy_set::get_y( double x ) const
{
  return get_y( x, _points->get_first() );
}

double ae_fuzzy_set::get_x( double y, ae_point_2d* left_point, ae_point_2d* right_point ) const
{
  assert( left_point != NULL );
  assert( right_point != NULL );
  assert( left_point->y != right_point->y );
  assert( ( y >= right_point->y && y <= left_point->y ) || ( y >= left_point->y && y <= right_point->y ) );

  double tmp = left_point->x + (right_point->x - left_point->x) * (y - left_point->y) / (right_point->y - left_point->y);

  if ( tmp < left_point->x ) return left_point->x;
  else if ( tmp > right_point->x ) return right_point->x;
  else return tmp;
}

ae_list_node<ae_point_2d*>* ae_fuzzy_set::create_interpolated_point( double x )
{
  return create_interpolated_point( x, _points->get_first() );
}


#endif // 0
#endif // __AE_FUZZY_SET_H__
