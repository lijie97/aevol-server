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

#include <stdio.h>
#include <list>

#include <fuzzy_X11.h>
#include "ae_macros.h"

using std::list;

namespace aevol {

void Fuzzy_X11::display(ae_X11_window* win,
                               color_map      color,
                               bool           fill /*= false*/,
                               bool           bold /*= false*/ ) {
  double y_min = Y_MIN - 0.1 * Y_MAX; // Yields a bottom margin
  double y_max = Y_MAX * 1.1;         // Yields a top margin
  double delta_x = X_MAX - X_MIN;
  double delta_y = y_max - y_min;
  
  int16_t cur_x;
  int16_t cur_y;
  int16_t next_x;
  int16_t next_y;
  
  for (list<Point>::iterator p = points.begin() ; p != prev(points.end()) ; ++p) {
    list<Point>::iterator q = next(p);
    
    // Display segment [p, q]
    cur_x   = (      (p->x -  X_MIN) / delta_x  ) * win->get_width();
    cur_y   = ( 1 - ((p->y -  y_min) / delta_y) ) * win->get_height();
    next_x  = (      (q->x - X_MIN) / delta_x  ) * win->get_width();
    next_y  = ( 1 - ((q->y - y_min) / delta_y) ) * win->get_height();
    
    if (fill) {
      char* fill_color;
      for ( int16_t i = cur_x ; i < next_x ; i++ ) {
        fill_color = ae_X11_window::get_color( ((double)i / win->get_width()) * (X_MAX - X_MIN) );
        win->draw_line( i, ( 1 - ((0 -  y_min) / delta_y) ) * win->get_height(),
                        i, cur_y + (((i - cur_x) * (next_y - cur_y)) / (next_x - cur_x)) , fill_color );
        delete [] fill_color;
      }
    }
    win->draw_line( cur_x, cur_y, next_x, next_y, color, bold );
  }
}
} // namespace aevol
