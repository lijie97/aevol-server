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
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_fuzzy_set_X11.h>




//##############################################################################
//                                                                             #
//                           Class ae_fuzzy_set_X11                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_fuzzy_set_X11::ae_fuzzy_set_X11( void ) : ae_fuzzy_set()
{
}

// =================================================================
//                             Destructors
// =================================================================
ae_fuzzy_set_X11::~ae_fuzzy_set_X11( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_fuzzy_set_X11::display( ae_X11_window* win, color_map color, bool fill /*= false*/, bool bold /*= false*/ )
{
  double y_min = Y_MIN - 0.1 * Y_MAX; // Yields a bottom margin
  double y_max = Y_MAX * 1.1;         // Yields a top margin
  double delta_x = X_MAX - X_MIN;
  double delta_y = y_max - y_min;
  
  ae_list_node* node        = _points->get_first();
  ae_point_2d*  cur_point   = NULL;
  ae_point_2d*  next_point  = NULL;
  int16_t cur_x;
  int16_t cur_y;
  int16_t next_x;
  int16_t next_y;
  
  
  while ( (node != NULL) && (node->get_next() != NULL) )
  {
    cur_point   = (ae_point_2d*)node->get_obj();
    next_point  = (ae_point_2d*)node->get_next()->get_obj();
    
    // Display segment [cur_point, next_point]
    cur_x   = (      (cur_point->x -  X_MIN) / delta_x  ) * win->get_width();
    cur_y   = ( 1 - ((cur_point->y -  y_min) / delta_y) ) * win->get_height();
    next_x  = (      (next_point->x - X_MIN) / delta_x  ) * win->get_width();
    next_y  = ( 1 - ((next_point->y - y_min) / delta_y) ) * win->get_height();
    
    if ( fill )
    {
      char* fill_color;
      for ( int16_t i = cur_x ; i < next_x ; i++ )
      {
        fill_color = ae_X11_window::get_color( ((double)i / win->get_width()) * (X_MAX - X_MIN) );
        win->draw_line( i, ( 1 - ((0 -  y_min) / delta_y) ) * win->get_height(),
                        i, cur_y + (((i - cur_x) * (next_y - cur_y)) / (next_x - cur_x)) , fill_color );
        delete [] fill_color;
      }
    }
    
    win->draw_line( cur_x, cur_y, next_x, next_y, color, bold );
    
    node = node->get_next();
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
