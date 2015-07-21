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
//*****************************************************************************




// =================================================================
//                              Libraries
// =================================================================
#include <stdio.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "ExpManager_X11.h"
#include "World.h"
#include "Point.h"
#include "Individual_X11.h"
#include "X11Window.h"
#include "Fuzzy.h"

namespace aevol {


// XCheckMaskEvent() doesn't get ClientMessage Events so use XCheckIfEvent()
// with this Predicate function as a work-around (ClientMessage events
// are needed in order to catch "WM_DELETE_WINDOW")
static Bool AlwaysTruePredicate (Display*, XEvent*, char*) { return True; }



//##############################################################################
//                                                                             #
//                          Class ExpManager_X11                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
#define NB_WIN INT32_C(7) // Number of windows that can be showed => CDS, RNA, phenotype, ...

// =================================================================
//                             Constructors
// =================================================================
ExpManager_X11::ExpManager_X11(void) : ExpManager()
{
  // Basic initializations
  _win      = NULL;
  _win_size = NULL;
  _win_pos  = NULL;
  _win_name = NULL;

  _display_on             = false;
  _handle_display_on_off  = false;

  // Initialize XLib stuff
  _display  = XOpenDisplay(NULL);
  if (_display == NULL)
  {
    printf("ERROR:\tCould not open connection to X server.\n");
    printf("\tIf you are using aevol through SSH, you may use ssh -X.\n");
    exit(EXIT_FAILURE);
  }
  _screen   = XDefaultScreen(_display);
  _atoms    = new Atom[2];
  _atoms[0] = XInternAtom(_display, "WM_DELETE_WINDOW", False);
  _atoms[1] = XInternAtom(_display, "WM_PROTOCOLS", False);
  set_codes();
}

// =================================================================
//                             Destructor
// =================================================================
ExpManager_X11::~ExpManager_X11(void)
{
  delete [] _key_codes;
  delete [] _atoms;

  for (int8_t i = 0 ; i < NB_WIN ; i++)
  {
    if (_win != NULL)
    {
      if (_win[i] != NULL) delete _win[i];
    }
    if (_win_size != NULL)
    {
      if (_win_size[i] != NULL) delete [] _win_size[i];
    }
    if (_win_pos != NULL)
    {
      if (_win_pos[i] != NULL) delete [] _win_pos[i];
    }
  }
  if (_win != NULL) delete [] _win;

  XCloseDisplay(_display);

  if (_win_name != NULL) delete [] _win_name;
  if (_win_size != NULL) delete [] _win_size;
  if (_win_pos != NULL) delete [] _win_pos;

}

// =================================================================
//                            Public Methods
// =================================================================


bool ExpManager_X11::quit_signal_received(void)
{
  return _quit_signal_received;
}

void ExpManager_X11::display(void)
{
  // ---------------------
  // 1) Handle user events
  // ---------------------
  if (_display_on) handle_events();


  // --------------------------------------------------
  // 2) Handle signal that toggle the display on or off
  // --------------------------------------------------
  if (_handle_display_on_off)
  {
    _handle_display_on_off = false;

    if (_display_on)
    // Display was "on", close all windows
    // (after saving their current size and position)
    {
      for (int8_t num_win = 0 ; num_win < NB_WIN ; num_win++)
      {
        if (_win[num_win] != NULL)
        {
          // 1) Save current window position and size
          Window aWindow; // Unused
          int x_return, y_return;
          int dest_x_return, dest_y_return;
          unsigned int border_width_return, depth_return; // Unused
          XGetGeometry(_display, _win[num_win]->get_window(), &aWindow,
              &x_return, &y_return,
              &_win_size[num_win][0], &_win_size[num_win][1],
              &border_width_return, &depth_return);
          XTranslateCoordinates(_display, _win[num_win]->get_window(),
              DefaultRootWindow(_display), 0, 0, &dest_x_return, &dest_y_return, &aWindow);

          _win_pos[num_win][0] = dest_x_return - x_return;
          _win_pos[num_win][1] = dest_y_return - y_return;

          // 2) Delete window
          delete _win[num_win];
          _win[num_win] = NULL;
        }
      }

      XFlush(_display);
      delete _win;
      _win = NULL;

      _display_on = false;
    }
    else // Display was "off", open windows
    {
      _show_window |= 1;
      _new_show_window = _show_window;

      // If it's the first time the display is switched on, initialize it.
      if (_win == NULL) initialize();

      for (int8_t i = 0 ; i < NB_WIN ; i++)
      {
        if (get_show_window(i))
        {
          _win[i] = new X11Window(_display, _screen, _atoms,
              _win_pos[i][0], _win_pos[i][1],
              _win_size[i][0], _win_size[i][1], _win_name[i]);
        }
      }

      _display_on = true;
    }
  }


  // ----------
  // 3) Display
  // ----------
  if (_display_on)
  {
    for (int8_t i = 0 ; i < NB_WIN ; i++)
    {
      if (get_show_window(i))
      {
        if (get_new_show_window(i))
        {
          draw_window(i);
        }
        else
        {
          refresh_window(i);
        }
      }
    }

    // Refresh all windows
    XFlush(_display);
  }
}

void ExpManager_X11::handle_events(void)
{
  XEvent event;
  int8_t win_number;


  // XCheckMaskEvent() doesn't get ClientMessage Events so use XCheckIfEvent()
  // with the custom AlwaysTruePredicate function as a work-around (ClientMessage
  // events are needed in order to catch "WM_DELETE_WINDOW")
  int iCurrEvent    = 0;
  int iIgnoreNoise  = 0;
  while(XCheckIfEvent(_display, &event, AlwaysTruePredicate, 0))
  {
    iCurrEvent ++;
    win_number = identify_window(event.xany.window);

    if(win_number == -1) continue;
    // We discard this event because it occurred on a destroyed window
    // (e.g. the user pressed F3 and then moved or resized the window,
    // before the window was destroyed)

    switch(event.type)
    {
      case ConfigureNotify :
      {
        _win[win_number]->resize(event.xconfigure.width, event.xconfigure.height);
        //~ _win[win_number]->repos(event.xconfigure.x, event.xconfigure.y);

        // Mark window as having to be entirely redrawn
        _new_show_window |= 1 << win_number;
        break;
      }
      case MapNotify :
      {
        draw_window(win_number);
        break;
      }
      case Expose:
      {
        if(iCurrEvent > iIgnoreNoise)
        {
          draw_window(win_number);
          iIgnoreNoise = iCurrEvent + XQLength(_display);
        }
        break;
      }
      case ClientMessage :
      {
        if (((Atom) event.xclient.data.l[0]) == _atoms[0]) // The user closed the window by clicking on the cross
        {
          // 1) Save current window position and size
          Window aWindow; // Unused
          int x_return, y_return;
          int dest_x_return, dest_y_return;
          unsigned int border_width_return, depth_return; // Unused
          XGetGeometry(_display, _win[win_number]->get_window(), &aWindow, &x_return, &y_return,
                        &_win_size[win_number][0], &_win_size[win_number][1], &border_width_return, &depth_return);
          XTranslateCoordinates(_display, _win[win_number]->get_window(), DefaultRootWindow(_display), 0, 0, &dest_x_return, &dest_y_return, &aWindow);

          _win_pos[win_number][0] = dest_x_return - x_return;
          _win_pos[win_number][1] = dest_y_return - y_return;

          // 2) Delete window and mark as "not to be shown"
          delete _win[win_number];
          _win[win_number] = NULL;
          _show_window &= ~(1 << win_number);

          // 3) If it was the main that was closed, turn display off.
          if (win_number == 0)
          {
            _handle_display_on_off = true;
          }
        }
        break;
      }
      case KeyPress :
      {
        // Not sure a switch would work on any platform => use ifs instead
        if (event.xkey.keycode == _key_codes[KEY_F1] ||
             event.xkey.keycode == _key_codes[KEY_F2] ||
             event.xkey.keycode == _key_codes[KEY_F3] ||
             event.xkey.keycode == _key_codes[KEY_F4] ||
             event.xkey.keycode == _key_codes[KEY_F5] ||
             event.xkey.keycode == _key_codes[KEY_F6])
        {
          int8_t num_win;

          // Not sure a switch would work on any platform => use ifs instead
          if (event.xkey.keycode == _key_codes[KEY_F1]) num_win = 1;
          else if (event.xkey.keycode == _key_codes[KEY_F2]) num_win = 2;
          else if (event.xkey.keycode == _key_codes[KEY_F3]) num_win = 3;
          else if (event.xkey.keycode == _key_codes[KEY_F4]) num_win = 4;
          else if (event.xkey.keycode == _key_codes[KEY_F5]) num_win = 5;
          else  num_win = 6; // case where (event.xkey.keycode == _key_codes[KEY_F6])

          if (get_show_window(num_win))
          {
            // 1) Save current window position and size
            Window aWindow; // Unused
            int x_return, y_return;
            int dest_x_return, dest_y_return;
            unsigned int border_width_return, depth_return; // Unused
            XGetGeometry(_display, _win[num_win]->get_window(), &aWindow, &x_return, &y_return,
                          &_win_size[num_win][0], &_win_size[num_win][1], &border_width_return, &depth_return);
            XTranslateCoordinates(_display, _win[num_win]->get_window(), DefaultRootWindow(_display), 0, 0, &dest_x_return, &dest_y_return, &aWindow);

            _win_pos[num_win][0] = dest_x_return - x_return;
            _win_pos[num_win][1] = dest_y_return - y_return;

            // 2) Delete window and mark as "not to be shown"
            delete _win[num_win];
            _win[num_win] = NULL;
            _show_window &= ~(1 << num_win);
          }
          else
          {
            _win[num_win] = new X11Window( _display, _screen, _atoms, _win_pos[num_win][0], _win_pos[num_win][1],
                                                _win_size[num_win][0], _win_size[num_win][1], _win_name[num_win]);
            _new_show_window |= 1 << num_win;
            _show_window |= _new_show_window;
            draw_window(num_win);
          }
        }
        else if (event.xkey.keycode == _key_codes[KEY_P])
        {
          printf(" P A U S E D \n");
          bool pause_key  = false;
          while (! pause_key)
          {
            if (XCheckIfEvent(_display, &event, AlwaysTruePredicate, 0))
            {
              if (event.xkey.keycode == _key_codes[KEY_P])
              {
                pause_key = true;
              }
            }
          }
        }
        else if ((event.xkey.state & ControlMask) && (event.xkey.keycode == _key_codes[KEY_Q]))
        {
          printf(" Q U I T   R E Q U E S T E D\n");
          _quit_signal_received = true;
          // We do not exit here, because it is cleaner to let the main program
          // interpret the signal and call the destructor of the simulation.
          // This ensures that the connection to the X server will be closed
          // in a clean way, which is much better if we want to run other
          // X applications afterwards!
        }

        break;
      }
      default :
      {
        // printf("event not handled \n");
        break;
      }
    }
  }
}

void ExpManager_X11::toggle_display_on_off(void)
{
  // Mark action to be done
  _handle_display_on_off = true;
}



void ExpManager_X11::display(X11Window * win,
    const Fuzzy& fuzzy,
    color_map color,
    bool fill /*= false*/,
    bool bold /*= false*/)
{
  double y_min = Y_MIN - 0.1 * Y_MAX; // Yields a bottom margin
  double y_max = Y_MAX * 1.1;         // Yields a top margin
  double delta_x = X_MAX - X_MIN;
  double delta_y = y_max - y_min;
  
  int16_t cur_x;
  int16_t cur_y;
  int16_t next_x;
  int16_t next_y;
  
  for (list<Point>::const_iterator p = fuzzy.get_points().begin() ; p != prev(fuzzy.get_points().end()) ; ++p) {
    list<Point>::const_iterator q = next(p);
    
    // Display segment [p, q]
    cur_x   = (     (p->x -  X_MIN) / delta_x ) * win->get_width();
    cur_y   = (1 - ((p->y -  y_min) / delta_y)) * win->get_height();
    next_x  = (     (q->x - X_MIN) / delta_x ) * win->get_width();
    next_y  = (1 - ((q->y - y_min) / delta_y)) * win->get_height();
    
    if (fill) {
      char* fill_color;
      for (int16_t i = cur_x ; i < next_x ; i++) {
        fill_color = X11Window::get_color(((double)i / win->get_width()) * (X_MAX - X_MIN));
        win->draw_line(i, (1 - ((0 -  y_min) / delta_y)) * win->get_height(),
                        i, cur_y + (((i - cur_x) * (next_y - cur_y)) / (next_x - cur_x)) , fill_color);
        delete [] fill_color;
      }
    }
    win->draw_line(cur_x, cur_y, next_x, next_y, color, bold);
  }
}

// Display a grid of values
void ExpManager_X11::display_grid(X11Window * win, double** cell_grid)
{
  // printf("display grid\n");
  char t[40];
  int nb_colors = 50; 
  
  sprintf(t, "Generation = %" PRId64, AeTime::get_time());
  win->draw_string(15, 15, t);
  
  
  const int grid_width  = get_grid_width();
  const int grid_height = get_grid_height();

  int nb_slots_in_a_row = (int) grid_height;
  int slot_width = 200/nb_slots_in_a_row;
  int x1 = 50 + 50 + slot_width/2;
  int y1 = 75 + 50 + slot_width/2;

  // create the colormap colors to be used for grid plotting
  int cell_size = 5;

  // draw the color scale for fitness
  int y_step_size = grid_height*cell_size/nb_colors;
  for (int i = 0; i  < nb_colors; i++)
  {
    win->fill_rectangle(x1 - 30, y1 - 80 + y_step_size * i,
                         cell_size * 5, y_step_size,
                         _col_map[nb_colors-1-i]);
  }

  // find min/max of the matrix
  double grid_max = 0;
  double grid_min = 1000000;
  for (int x = 0 ; x < grid_width ; x++)
  {
    for (int y = 0 ; y < grid_height ; y++)
    {
       if (cell_grid[x][y] > grid_max) {grid_max = cell_grid[x][y];}
       if (cell_grid[x][y] < grid_min) {grid_min = cell_grid[x][y];}
     }
  }
  double col_sec_interval = (grid_max - grid_min)/49;

  char scale_txt[40];
  sprintf(scale_txt,"%.2e", grid_max);
  win->draw_string(x1-80, y1-80,scale_txt);
  sprintf(scale_txt,"%.2e", grid_min);
  win->draw_string(x1-80, y1-80+grid_height*cell_size,scale_txt);

  for (int x = 0; x < grid_width; x++)
  {
    for (int y = 0; y < grid_height; y++)
    {
      char * col_string;
      // calculate the color
      int new_col;
      if (col_sec_interval==0)
      {
        new_col = 0;
      }
      else
      {
        new_col = (int) floor((cell_grid[x][y] - grid_min) / col_sec_interval);
      }
      col_string = _col_map[new_col];

      // draw a colored rectangle for each cell
      win->fill_rectangle(x1 + 50 + x*cell_size, y1 - 80 + y*cell_size, cell_size, cell_size, col_string);
    }
  }
}


// =================================================================
//                           Protected Methods
// =================================================================
void ExpManager_X11::initialize(bool with_grid /*= false*/, bool with_plasmids /*= false*/)
{
  // Initialize window structures
  _win      = new X11Window * [NB_WIN];
  _win_size = new unsigned int* [NB_WIN];
  _win_pos  = new int* [NB_WIN];

  for (int8_t i = 0 ; i < NB_WIN ; i++)
  {
    _win[i] = NULL;

    // Default values
    _win_size[i] = new unsigned int[2];
    _win_size[i][0] = 300;
    _win_size[i][1] = 300;
    _win_pos[i] = new int[2];
    _win_pos[i][0]  = 0;
    _win_pos[i][1]  = 0;
  }

  // Set phenotype window width
  _win_size[1][0] = 600;

  // Set CDS and RNA window width
  if (with_plasmids)
  {
    _win_size[2][0] = 600;
    _win_size[3][0] = 600;
  }

  // Set initial positions if screen is large enough
  if (with_plasmids && with_grid)
  {
    //if (XDisplayWidth(_display, _screen) >= 900 && XDisplayHeight(_display, _screen) >= 650)
    {
      _win_pos[0][0]  = 0;
      _win_pos[0][1]  = 0;
      _win_pos[1][0]  = 300;
      _win_pos[1][1]  = 0;
      _win_pos[2][0]  = 0;
      _win_pos[2][1]  = 350;
      _win_pos[3][0]  = 0;
      _win_pos[3][1]  = 700;
    }
  }
  else if (with_plasmids)
  {
    //if (XDisplayWidth(_display, _screen) >= 900 && XDisplayHeight(_display, _screen) >= 650)
    {
      _win_pos[0][0]  = 0;
      _win_pos[0][1]  = 0;
      _win_pos[1][0]  = 300;
      _win_pos[1][1]  = 0;
      _win_pos[2][0]  = 0;
      _win_pos[2][1]  = 350;
      _win_pos[3][0]  = 0;
      _win_pos[3][1]  = 700;
    }
  }
  else if (with_grid)
  {
    //if (XDisplayWidth(_display, _screen) >= 900 && XDisplayHeight(_display, _screen) >= 650)
    {
      _win_pos[0][0]  = 0;
      _win_pos[0][1]  = 0;
      _win_pos[1][0]  = 300;
      _win_pos[1][1]  = 0;
      _win_pos[2][0]  = 0;
      _win_pos[2][1]  = 350;
      _win_pos[3][0]  = 300;
      _win_pos[3][1]  = 350;
    }
  }
  else // (! with_plasmids && ! with_grid)
  {
    //if (XDisplayWidth(_display, _screen) >= 900 && XDisplayHeight(_display, _screen) >= 650)
    {
      _win_pos[0][0]  = 0;
      _win_pos[0][1]  = 0;
      _win_pos[1][0]  = 300;
      _win_pos[1][1]  = 0;
      _win_pos[2][0]  = 0;
      _win_pos[2][1]  = 350;
      _win_pos[3][0]  = 300;
      _win_pos[3][1]  = 350;
    }
  }


  // Visible windows at the beginning of the run
  if (with_grid)
  {
    _show_window  = 0x007F; // hex for bin 1111111  => show first 7 windows
  }
  else
  {
    _show_window  = 0x000F;   // hex for bin 1111   => show first 4 windows
  }
  _new_show_window = _show_window;


  _win_name = new char*[NB_WIN];
  if (with_grid)  _win_name[0] = (char*) "Population grid";
  else              _win_name[0] = (char*) "Population";

  _win_name[1] = (char*) "Phenotypic profile";
  _win_name[2] = (char*) "Genes";
  _win_name[3] = (char*) "RNAs";
  _win_name[4] = (char*) "Secreted compound present";
  _win_name[5] = (char*) "Metabolic fitness";
  _win_name[6] = (char*) "Current secretion";

  compute_colormap();
}

int8_t ExpManager_X11::identify_window(Window winID)
{
  for (int8_t i = 0 ; i < NB_WIN ; i++)
  {
    if (_win[i] != NULL)
    {
      if (_win[i]->get_window() == winID) return i;
    }
  }

  return -1;
}

void ExpManager_X11::draw_window(int8_t win_number)
{
  if (_win[win_number] == NULL)
  {
    fprintf(stderr, "Error: cannot draw this window, it doesn't exist.\n");
    return;
  }

  X11Window * cur_win = _win[win_number];

  switch (win_number)
  {
    case 0:
    {
      break;
    }

    case 1:
    {
      cur_win->blacken();

      // Display colour bar
      char* color;
      for (int16_t i = 0 ; i < cur_win->get_width() ; i++)
      {
        color = X11Window::get_color(((double)i / cur_win->get_width()) * (X_MAX - X_MIN));
        //~ cur_win->draw_line(i, 0, i, cur_win->get_height() / 20, color);
        cur_win->draw_line(i, cur_win->get_height() * 19 / 20, i, cur_win->get_height(), color);
        delete [] color;
      }

      break;
    }

    case 2:
    {
      break;
    }

    case 3:
    {
      break;
    }

    // Amount of secreted compound present at each grid location
    case 4:
    {
      break;
    }

    // Metabolic fitness  grid
    case 5:
    {
      break;
    }

    // Current secretion (how much is secreted by each organism)
    case 6:
    {
      break;
    }
  }

  refresh_window(win_number);
  _new_show_window &= ~(1 << win_number);

  XFlush(_display);
}

void ExpManager_X11::refresh_window(int8_t win_number) {
  if (_win[win_number] == NULL)
  {
    fprintf(stderr, "Error: cannot draw this window, it doesn't exist.\n");
    return;
  }

  X11Window * cur_win = _win[win_number];

  switch (win_number)
  {
    // Main window (population)
    case 0 :
    {
      cur_win->blacken();

      double** grid = world()->get_total_fitness_grid();
      display_grid(cur_win, grid);

      // Has been allocated in ae_spatial_structure::get_total_fitness_grid()
      for (int16_t x = 0 ; x < get_grid_width() ; x++)
      {
        delete [] grid[x];
      }
      delete [] grid;
      break;
    }

    // Display phenotypes and phenotypic target
    case 1 :
    {
      // Blacken all the window except the colour bar
      cur_win->fill_rectangle(0, 0, cur_win->get_width(), cur_win->get_height() * 19 / 20, BLACK);

      // Get phenotypic target shorthand
      const PhenotypicTarget& phenotypic_target = get_best_indiv()->phenotypic_target();

      // Mark all the non-metabolic segments (paint them in grey)
      if (phenotypic_target.nb_segments() > 1)
      {
        PhenotypicSegment ** segments = phenotypic_target.segments();

        for (size_t i = 0 ; i < static_cast<size_t>(phenotypic_target.nb_segments()) ; i++) {
          if (segments[i]->feature != METABOLISM)
          {
            if (segments[i]->feature == NEUTRAL)
            {
              cur_win->fill_rectangle( cur_win->get_width() * segments[i]->start / (X_MAX-X_MIN), 0.0,
                                        cur_win->get_width() * (segments[i]->stop - segments[i]->start) / (X_MAX-X_MIN),
                                        cur_win->get_height() * 19 / 20, DARKER_GREY);
            }
            else
            {
              cur_win->fill_rectangle( cur_win->get_width() * segments[i]->start / (X_MAX-X_MIN), 0.0,
                                        cur_win->get_width() * (segments[i]->stop - segments[i]->start) / (X_MAX-X_MIN),
                                        cur_win->get_height() * 19 / 20, GREY);
            }
          }
        }
      }

      // Display all the phenotypes (blue)
//      for (const auto& indiv: get_indivs())
//      {
//        display(cur_win, *(indiv->get_phenotype()), BLUE);
//        if (indiv->get_allow_plasmids())
//        {
//          display(cur_win, *(indiv->get_genetic_unit(0).get_phenotypic_contribution()), YELLOW);
//          display(cur_win, *(indiv->get_genetic_unit(1).get_phenotypic_contribution()), GREEN);
//        }
//      }

      // Display best indiv's phenotype (white)
      display(cur_win, *(get_best_indiv()->get_phenotype()), WHITE, true);

      // Display phenotypic target (red)
      display(cur_win, phenotypic_target, RED, false, true);
    }
    break;

    // Display genes
    case 2 :
    {
      cur_win->blacken();

      Individual_X11 * indiv = dynamic_cast<Individual_X11 *>(get_best_indiv());
      indiv->display_cdss(cur_win);
    }
    break;

    // Display RNAs
    case 3 :
    {
      cur_win->blacken();

      Individual_X11 * indiv = dynamic_cast<Individual_X11 *>(get_best_indiv());
      indiv->display_rnas(cur_win);
    }
    break;

    // Display the amount of secreted compound present at each location
    case 4 :
    {
      cur_win->blacken();

      display_grid(cur_win, world()->get_secretion_present_grid());
    }
    break;

    // Display the metabolic fitness grid
    case 5 :
    {
      cur_win->blacken();

      display_grid(cur_win, world()->get_metabolic_fitness_grid());
    }
    break;

    // display current secretion (how much is secreted by each organism)
    case 6:
    {
      cur_win->blacken();

      display_grid(cur_win, world()->get_secreted_amount_grid());
    }
    break;
  }

  XFlush(_display);
}





void ExpManager_X11::set_codes(void) {
  _key_codes = new KeyCode[50];
  assert(_key_codes);

  _key_codes[KEY_ESCAPE]  = XKeysymToKeycode(_display, XK_Escape);
  _key_codes[KEY_F1]      = XKeysymToKeycode(_display, XK_F1);
  _key_codes[KEY_F2]      = XKeysymToKeycode(_display, XK_F2);
  _key_codes[KEY_F3]      = XKeysymToKeycode(_display, XK_F3);
  _key_codes[KEY_F4]      = XKeysymToKeycode(_display, XK_F4);
  _key_codes[KEY_F5]      = XKeysymToKeycode(_display, XK_F5);
  _key_codes[KEY_F6]      = XKeysymToKeycode(_display, XK_F6);
  _key_codes[KEY_F7]      = XKeysymToKeycode(_display, XK_F7);
  _key_codes[KEY_F8]      = XKeysymToKeycode(_display, XK_F8);
  _key_codes[KEY_F9]      = XKeysymToKeycode(_display, XK_F9);
  _key_codes[KEY_F10]     = XKeysymToKeycode(_display, XK_F10);
  _key_codes[KEY_F11]     = XKeysymToKeycode(_display, XK_F11);
  _key_codes[KEY_F12]     = XKeysymToKeycode(_display, XK_F12);
  _key_codes[KEY_A]       = XKeysymToKeycode(_display, XK_A);
  _key_codes[KEY_Q]       = XKeysymToKeycode(_display, XK_Q);
  _key_codes[KEY_W]       = XKeysymToKeycode(_display, XK_W);
  _key_codes[KEY_Z]       = XKeysymToKeycode(_display, XK_Z);
  _key_codes[KEY_S]       = XKeysymToKeycode(_display, XK_S);
  _key_codes[KEY_X]       = XKeysymToKeycode(_display, XK_X);
  _key_codes[KEY_E]       = XKeysymToKeycode(_display, XK_E);
  _key_codes[KEY_D]       = XKeysymToKeycode(_display, XK_D);
  _key_codes[KEY_C]       = XKeysymToKeycode(_display, XK_C);
  _key_codes[KEY_R]       = XKeysymToKeycode(_display, XK_R);
  _key_codes[KEY_F]       = XKeysymToKeycode(_display, XK_F);
  _key_codes[KEY_V]       = XKeysymToKeycode(_display, XK_V);
  _key_codes[KEY_T]       = XKeysymToKeycode(_display, XK_T);
  _key_codes[KEY_G]       = XKeysymToKeycode(_display, XK_G);
  _key_codes[KEY_B]       = XKeysymToKeycode(_display, XK_B);
  _key_codes[KEY_Y]       = XKeysymToKeycode(_display, XK_Y);
  _key_codes[KEY_H]       = XKeysymToKeycode(_display, XK_H);
  _key_codes[KEY_N]       = XKeysymToKeycode(_display, XK_N);
  _key_codes[KEY_U]       = XKeysymToKeycode(_display, XK_U);
  _key_codes[KEY_J]       = XKeysymToKeycode(_display, XK_J);
  _key_codes[KEY_I]       = XKeysymToKeycode(_display, XK_I);
  _key_codes[KEY_K]       = XKeysymToKeycode(_display, XK_K);
  _key_codes[KEY_O]       = XKeysymToKeycode(_display, XK_O);
  _key_codes[KEY_L]       = XKeysymToKeycode(_display, XK_L);
  _key_codes[KEY_P]       = XKeysymToKeycode(_display, XK_P);
  _key_codes[KEY_M]       = XKeysymToKeycode(_display, XK_M);
  _key_codes[KEY_1]       = XKeysymToKeycode(_display, XK_1);
  _key_codes[KEY_2]       = XKeysymToKeycode(_display, XK_2);
  _key_codes[KEY_3]       = XKeysymToKeycode(_display, XK_3);
  _key_codes[KEY_4]       = XKeysymToKeycode(_display, XK_4);
  _key_codes[KEY_5]       = XKeysymToKeycode(_display, XK_5);
  _key_codes[KEY_6]       = XKeysymToKeycode(_display, XK_6);
  _key_codes[KEY_7]       = XKeysymToKeycode(_display, XK_7);
  _key_codes[KEY_8]       = XKeysymToKeycode(_display, XK_8);
  _key_codes[KEY_9]       = XKeysymToKeycode(_display, XK_9);
}


void ExpManager_X11::compute_colormap(void) {
  _col_map = {
    (char*)"RGBi:1.0/0.0/0.0",
    (char*)"RGBi:1.0/0.1/0.0",   
    (char*)"RGBi:1.0/0.2/0.0",
    (char*)"RGBi:1.0/0.3/0.0",
    (char*)"RGBi:1.0/0.4/0.0",
    (char*)"RGBi:1.0/0.5/0.0",
    (char*)"RGBi:1.0/0.6/0.0",
    (char*)"RGBi:1.0/0.7/0.0",
    (char*)"RGBi:1.0/0.8/0.0",
    (char*)"RGBi:1.0/0.9/0.0",

    (char*)"RGBi:0.9/1.0/0.0",
    (char*)"RGBi:0.8/1.0/0.0",
    (char*)"RGBi:0.7/1.0/0.0",
    (char*)"RGBi:0.6/1.0/0.0",
    (char*)"RGBi:0.5/1.0/0.0",
    (char*)"RGBi:0.4/1.0/0.0",
    (char*)"RGBi:0.3/1.0/0.0",
    (char*)"RGBi:0.2/1.0/0.0",
    (char*)"RGBi:0.1/1.0/0.0",
    (char*)"RGBi:0.0/1.0/0.0",

    (char*)"RGBi:0.0/1.0/0.1",
    (char*)"RGBi:0.0/1.0/0.2",
    (char*)"RGBi:0.0/1.0/0.3",
    (char*)"RGBi:0.0/1.0/0.4",
    (char*)"RGBi:0.0/1.0/0.5",
    (char*)"RGBi:0.0/1.0/0.6",
    (char*)"RGBi:0.0/1.0/0.7",
    (char*)"RGBi:0.0/1.0/0.8",
    (char*)"RGBi:0.0/1.0/0.9",
    (char*)"RGBi:0.0/1.0/1.0",

    (char*)"RGBi:0.0/0.9/1.0",
    (char*)"RGBi:0.0/0.8/1.0",
    (char*)"RGBi:0.0/0.7/1.0",
    (char*)"RGBi:0.0/0.6/1.0",
    (char*)"RGBi:0.0/0.5/1.0",
    (char*)"RGBi:0.0/0.4/1.0",
    (char*)"RGBi:0.0/0.3/1.0",
    (char*)"RGBi:0.0/0.2/1.0",
    (char*)"RGBi:0.0/0.1/1.0",
    (char*)"RGBi:0.0/0.0/1.0",

    (char*)"RGBi:0.1/0.0/1.0",
    (char*)"RGBi:0.2/0.0/1.0",
    (char*)"RGBi:0.3/0.0/1.0",
    (char*)"RGBi:0.4/0.0/1.0",
    (char*)"RGBi:0.5/0.0/1.0",
    (char*)"RGBi:0.6/0.0/1.0",
    (char*)"RGBi:0.7/0.0/1.0",
    (char*)"RGBi:0.8/0.0/1.0",
    (char*)"RGBi:0.9/0.0/1.0",
    (char*)"RGBi:1.0/0.0/1.0",
  };

}
} // namespace aevol
