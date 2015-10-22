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
ExpManager_X11::ExpManager_X11() : ExpManager()
{
  // Basic initializations
  win_ = NULL;
  win_size_ = NULL;
  win_pos_ = NULL;
  win_name_ = NULL;

  display_on_ = false;
  handle_display_on_off_ = false;

  // Initialize XLib stuff
  display_ = XOpenDisplay(NULL);
  if (display_ == NULL)
  {
    printf("ERROR:\tCould not open connection to X server.\n");
    printf("\tIf you are using aevol through SSH, you may use ssh -X.\n");
    exit(EXIT_FAILURE);
  }
  screen_ = XDefaultScreen(display_);
  atoms_ = new Atom[2];
  atoms_[0] = XInternAtom(display_, "WM_DELETE_WINDOW", False);
  atoms_[1] = XInternAtom(display_, "WM_PROTOCOLS", False);
  set_codes();
}

// =================================================================
//                             Destructor
// =================================================================
ExpManager_X11::~ExpManager_X11()
{
  delete [] key_codes_;
  delete [] atoms_;

  for (int8_t i = 0 ; i < NB_WIN ; i++)
  {
    if (win_ != NULL)
    {
      if (win_[i] != NULL) delete win_[i];
    }
    if (win_size_ != NULL)
    {
      if (win_size_[i] != NULL) delete [] win_size_[i];
    }
    if (win_pos_ != NULL)
    {
      if (win_pos_[i] != NULL) delete [] win_pos_[i];
    }
  }
  if (win_ != NULL) delete [] win_;

  XCloseDisplay(display_);

  if (win_name_ != NULL) delete [] win_name_;
  if (win_size_ != NULL) delete [] win_size_;
  if (win_pos_ != NULL) delete [] win_pos_;

}

// =================================================================
//                            Public Methods
// =================================================================
void ExpManager_X11::display()
{
  // ---------------------
  // 1) Handle user events
  // ---------------------
  if (display_on_) handle_events();


  // --------------------------------------------------
  // 2) Handle signal that toggle the display on or off
  // --------------------------------------------------
  if (handle_display_on_off_)
  {
    handle_display_on_off_ = false;

    if (display_on_)
    // Display was "on", close all windows
    // (after saving their current size and position)
    {
      for (int8_t num_win = 0 ; num_win < NB_WIN ; num_win++)
      {
        if (win_[num_win] != NULL)
        {
          // 1) Save current window position and size
          Window aWindow; // Unused
          int x_return, y_return;
          int dest_x_return, dest_y_return;
          unsigned int border_width_return, depth_return; // Unused
          XGetGeometry(display_, win_[num_win]->get_window(), &aWindow,
              &x_return, &y_return,
              &win_size_[num_win][0], &win_size_[num_win][1],
              &border_width_return, &depth_return);
          XTranslateCoordinates(display_, win_[num_win]->get_window(),
              DefaultRootWindow(display_), 0, 0, &dest_x_return, &dest_y_return, &aWindow);

          win_pos_[num_win][0] = dest_x_return - x_return;
          win_pos_[num_win][1] = dest_y_return - y_return;

          // 2) Delete window
          delete win_[num_win];
          win_[num_win] = NULL;
        }
      }

      XFlush(display_);
      delete win_;
      win_ = NULL;

      display_on_ = false;
    }
    else // Display was "off", open windows
    {
      show_window_ |= 1;
      new_show_window_ = show_window_;

      // If it's the first time the display is switched on, initialize it.
      if (win_ == NULL) initialize();

      for (int8_t i = 0 ; i < NB_WIN ; i++)
      {
        if (show_window(i))
        {
          win_[i] = new X11Window(display_, screen_, atoms_,
              win_pos_[i][0], win_pos_[i][1],
              win_size_[i][0], win_size_[i][1], win_name_[i]);
        }
      }

      display_on_ = true;
    }
  }


  // ----------
  // 3) Display
  // ----------
  if (display_on_)
  {
    for (int8_t i = 0 ; i < NB_WIN ; i++)
    {
      if (show_window(i))
      {
        if (new_show_window(i))
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
    XFlush(display_);
  }
}

void ExpManager_X11::handle_events()
{
  XEvent event;
  int8_t win_number;


  // XCheckMaskEvent() doesn't get ClientMessage Events so use XCheckIfEvent()
  // with the custom AlwaysTruePredicate function as a work-around (ClientMessage
  // events are needed in order to catch "WM_DELETE_WINDOW")
  int iCurrEvent    = 0;
  int iIgnoreNoise  = 0;
  while(XCheckIfEvent(display_, &event, AlwaysTruePredicate, 0))
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
        win_[win_number]->resize(event.xconfigure.width, event.xconfigure.height);
        //~ _win[win_number]->repos(event.xconfigure.x, event.xconfigure.y);

        // Mark window as having to be entirely redrawn
        new_show_window_ |= 1 << win_number;
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
          iIgnoreNoise = iCurrEvent + XQLength(display_);
        }
        break;
      }
      case ClientMessage :
      {
        if (((Atom) event.xclient.data.l[0]) == atoms_[0]) // The user closed the window by clicking on the cross
        {
          // 1) Save current window position and size
          Window aWindow; // Unused
          int x_return, y_return;
          int dest_x_return, dest_y_return;
          unsigned int border_width_return, depth_return; // Unused
          XGetGeometry(display_, win_[win_number]->get_window(), &aWindow, &x_return, &y_return,
                        &win_size_[win_number][0], &win_size_[win_number][1], &border_width_return, &depth_return);
          XTranslateCoordinates(display_, win_[win_number]->get_window(), DefaultRootWindow(display_), 0, 0, &dest_x_return, &dest_y_return, &aWindow);

          win_pos_[win_number][0] = dest_x_return - x_return;
          win_pos_[win_number][1] = dest_y_return - y_return;

          // 2) Delete window and mark as "not to be shown"
          delete win_[win_number];
          win_[win_number] = NULL;
          show_window_ &= ~(1 << win_number);

          // 3) If it was the main that was closed, turn display off.
          if (win_number == 0)
          {
            handle_display_on_off_ = true;
          }
        }
        break;
      }
      case KeyPress :
      {
        // Not sure a switch would work on any platform => use ifs instead
        if (event.xkey.keycode == key_codes_[KEY_F1] ||
             event.xkey.keycode == key_codes_[KEY_F2] ||
             event.xkey.keycode == key_codes_[KEY_F3] ||
             event.xkey.keycode == key_codes_[KEY_F4] ||
             event.xkey.keycode == key_codes_[KEY_F5] ||
             event.xkey.keycode == key_codes_[KEY_F6])
        {
          int8_t num_win;

          // Not sure a switch would work on any platform => use ifs instead
          if (event.xkey.keycode == key_codes_[KEY_F1]) num_win = 1;
          else if (event.xkey.keycode == key_codes_[KEY_F2]) num_win = 2;
          else if (event.xkey.keycode == key_codes_[KEY_F3]) num_win = 3;
          else if (event.xkey.keycode == key_codes_[KEY_F4]) num_win = 4;
          else if (event.xkey.keycode == key_codes_[KEY_F5]) num_win = 5;
          else  num_win = 6; // case where (event.xkey.keycode == key_codes_[KEY_F6])

          if (show_window(num_win))
          {
            // 1) Save current window position and size
            Window aWindow; // Unused
            int x_return, y_return;
            int dest_x_return, dest_y_return;
            unsigned int border_width_return, depth_return; // Unused
            XGetGeometry(display_, win_[num_win]->get_window(), &aWindow, &x_return, &y_return,
                          &win_size_[num_win][0], &win_size_[num_win][1], &border_width_return, &depth_return);
            XTranslateCoordinates(display_, win_[num_win]->get_window(), DefaultRootWindow(display_), 0, 0, &dest_x_return, &dest_y_return, &aWindow);

            win_pos_[num_win][0] = dest_x_return - x_return;
            win_pos_[num_win][1] = dest_y_return - y_return;

            // 2) Delete window and mark as "not to be shown"
            delete win_[num_win];
            win_[num_win] = NULL;
            show_window_ &= ~(1 << num_win);
          }
          else
          {
            win_[num_win] = new X11Window(display_, screen_, atoms_, win_pos_[num_win][0], win_pos_[num_win][1],
                                                win_size_[num_win][0], win_size_[num_win][1], win_name_[num_win]);
            new_show_window_ |= 1 << num_win;
            show_window_ |= new_show_window_;
            draw_window(num_win);
          }
        }
        else if (event.xkey.keycode == key_codes_[KEY_P])
        {
          printf(" P A U S E D \n");
          bool pause_key  = false;
          while (! pause_key)
          {
            if (XCheckIfEvent(display_, &event, AlwaysTruePredicate, 0))
            {
              if (event.xkey.keycode == key_codes_[KEY_P])
              {
                pause_key = true;
              }
            }
          }
        }
        else if ((event.xkey.state & ControlMask) && (event.xkey.keycode == key_codes_[KEY_Q]))
        {
          printf(" Q U I T   R E Q U E S T E D\n");
          quit_signal_received_ = true;
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

void ExpManager_X11::toggle_display_on_off()
{
  // Mark action to be done
  handle_display_on_off_ = true;
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
  
  
  const int grid_width_ = grid_width();
  const int grid_height_ = grid_height();

  int nb_slots_in_a_row = (int) grid_height_;
  int slot_width = 200/nb_slots_in_a_row;
  int x1 = 50 + 50 + slot_width/2;
  int y1 = 75 + 50 + slot_width/2;

  // create the colormap colors to be used for grid plotting
  int cell_size = 5;

  // draw the color scale for fitness
  int y_step_size = grid_height_ *cell_size/nb_colors;
  for (int i = 0; i  < nb_colors; i++)
  {
    win->fill_rectangle(x1 - 30, y1 - 80 + y_step_size * i,
                         cell_size * 5, y_step_size,
                         _col_map[nb_colors-1-i]);
  }

  // find min/max of the matrix
  double grid_max = 0;
  double grid_min = 1000000;
  for (int x = 0 ; x < grid_width_; x++)
  {
    for (int y = 0 ; y < grid_height_; y++)
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
  win->draw_string(x1-80, y1-80+ grid_height_ *cell_size,scale_txt);

  for (int x = 0; x < grid_width_; x++)
  {
    for (int y = 0; y < grid_height_; y++)
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
  win_ = new X11Window * [NB_WIN];
  win_size_ = new unsigned int* [NB_WIN];
  win_pos_ = new int* [NB_WIN];

  for (int8_t i = 0 ; i < NB_WIN ; i++)
  {
    win_[i] = NULL;

    // Default values
    win_size_[i] = new unsigned int[2];
    win_size_[i][0] = 300;
    win_size_[i][1] = 300;
    win_pos_[i] = new int[2];
    win_pos_[i][0]  = 0;
    win_pos_[i][1]  = 0;
  }

  // Set phenotype window width
  win_size_[1][0] = 600;

  // Set CDS and RNA window width
  if (with_plasmids)
  {
    win_size_[2][0] = 600;
    win_size_[3][0] = 600;
  }

  // Set initial positions if screen is large enough
  if (with_plasmids && with_grid)
  {
    //if (XDisplayWidth(display_, screen_) >= 900 && XDisplayHeight(display_, screen_) >= 650)
    {
      win_pos_[0][0]  = 0;
      win_pos_[0][1]  = 0;
      win_pos_[1][0]  = 300;
      win_pos_[1][1]  = 0;
      win_pos_[2][0]  = 0;
      win_pos_[2][1]  = 350;
      win_pos_[3][0]  = 0;
      win_pos_[3][1]  = 700;
    }
  }
  else if (with_plasmids)
  {
    //if (XDisplayWidth(display_, screen_) >= 900 && XDisplayHeight(display_, screen_) >= 650)
    {
      win_pos_[0][0]  = 0;
      win_pos_[0][1]  = 0;
      win_pos_[1][0]  = 300;
      win_pos_[1][1]  = 0;
      win_pos_[2][0]  = 0;
      win_pos_[2][1]  = 350;
      win_pos_[3][0]  = 0;
      win_pos_[3][1]  = 700;
    }
  }
  else if (with_grid)
  {
    //if (XDisplayWidth(display_, screen_) >= 900 && XDisplayHeight(display_, screen_) >= 650)
    {
      win_pos_[0][0]  = 0;
      win_pos_[0][1]  = 0;
      win_pos_[1][0]  = 300;
      win_pos_[1][1]  = 0;
      win_pos_[2][0]  = 0;
      win_pos_[2][1]  = 350;
      win_pos_[3][0]  = 300;
      win_pos_[3][1]  = 350;
    }
  }
  else // (! with_plasmids && ! with_grid)
  {
    //if (XDisplayWidth(display_, screen_) >= 900 && XDisplayHeight(display_, screen_) >= 650)
    {
      win_pos_[0][0]  = 0;
      win_pos_[0][1]  = 0;
      win_pos_[1][0]  = 300;
      win_pos_[1][1]  = 0;
      win_pos_[2][0]  = 0;
      win_pos_[2][1]  = 350;
      win_pos_[3][0]  = 300;
      win_pos_[3][1]  = 350;
    }
  }


  // Visible windows at the beginning of the run
  if (with_grid)
  {
    show_window_ = 0x007F; // hex for bin 1111111  => show first 7 windows
  }
  else
  {
    show_window_ = 0x000F;   // hex for bin 1111   => show first 4 windows
  }
  new_show_window_ = show_window_;


  win_name_ = new char*[NB_WIN];
  if (with_grid)  win_name_[0] = (char*) "Population grid";
  else              win_name_[0] = (char*) "Population";

  win_name_[1] = (char*) "Phenotypic profile";
  win_name_[2] = (char*) "Genes";
  win_name_[3] = (char*) "RNAs";
  win_name_[4] = (char*) "Secreted compound present";
  win_name_[5] = (char*) "Metabolic fitness";
  win_name_[6] = (char*) "Current secretion";

  compute_colormap();
}

int8_t ExpManager_X11::identify_window(Window winID)
{
  for (int8_t i = 0 ; i < NB_WIN ; i++)
  {
    if (win_[i] != NULL)
    {
      if (win_[i]->get_window() == winID) return i;
    }
  }

  return -1;
}

void ExpManager_X11::draw_window(int8_t win_number)
{
  if (win_[win_number] == NULL)
  {
    fprintf(stderr, "Error: cannot draw this window, it doesn't exist.\n");
    return;
  }

  X11Window * cur_win = win_[win_number];

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
  new_show_window_ &= ~(1 << win_number);

  XFlush(display_);
}

void ExpManager_X11::refresh_window(int8_t win_number) {
  if (win_[win_number] == NULL)
  {
    fprintf(stderr, "Error: cannot draw this window, it doesn't exist.\n");
    return;
  }

  X11Window * cur_win = win_[win_number];

  switch (win_number)
  {
    // Main window (population)
    case 0 :
    {
      cur_win->blacken();

      double** grid = world()->get_total_fitness_grid();
      display_grid(cur_win, grid);

      // Has been allocated in ae_spatial_structure::get_total_fitness_grid()
      for (int16_t x = 0 ; x < grid_width() ; x++)
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
      const PhenotypicTarget& phenotypic_target = best_indiv()->phenotypic_target();

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
//      for (const auto& indiv: indivs())
//      {
//        display(cur_win, *(indiv->get_phenotype()), BLUE);
//        if (indiv->get_allow_plasmids())
//        {
//          display(cur_win, *(indiv->get_genetic_unit(0).get_phenotypic_contribution()), YELLOW);
//          display(cur_win, *(indiv->get_genetic_unit(1).get_phenotypic_contribution()), GREEN);
//        }
//      }

      // Display best indiv's phenotype (white)
      display(cur_win, *(best_indiv()->get_phenotype()), WHITE, true);

      // Display phenotypic target (red)
      display(cur_win, phenotypic_target, RED, false, true);
    }
    break;

    // Display genes
    case 2 :
    {
      cur_win->blacken();

      Individual_X11 * indiv = dynamic_cast<Individual_X11 *>(best_indiv());
      indiv->display_cdss(cur_win);
    }
    break;

    // Display RNAs
    case 3 :
    {
      cur_win->blacken();

      Individual_X11 * indiv = dynamic_cast<Individual_X11 *>(best_indiv());
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

  XFlush(display_);
}

void ExpManager_X11::set_codes() {
  key_codes_ = new KeyCode[50];
  assert(key_codes_);

  key_codes_[KEY_ESCAPE]  = XKeysymToKeycode(display_, XK_Escape);
  key_codes_[KEY_F1]      = XKeysymToKeycode(display_, XK_F1);
  key_codes_[KEY_F2]      = XKeysymToKeycode(display_, XK_F2);
  key_codes_[KEY_F3]      = XKeysymToKeycode(display_, XK_F3);
  key_codes_[KEY_F4]      = XKeysymToKeycode(display_, XK_F4);
  key_codes_[KEY_F5]      = XKeysymToKeycode(display_, XK_F5);
  key_codes_[KEY_F6]      = XKeysymToKeycode(display_, XK_F6);
  key_codes_[KEY_F7]      = XKeysymToKeycode(display_, XK_F7);
  key_codes_[KEY_F8]      = XKeysymToKeycode(display_, XK_F8);
  key_codes_[KEY_F9]      = XKeysymToKeycode(display_, XK_F9);
  key_codes_[KEY_F10]     = XKeysymToKeycode(display_, XK_F10);
  key_codes_[KEY_F11]     = XKeysymToKeycode(display_, XK_F11);
  key_codes_[KEY_F12]     = XKeysymToKeycode(display_, XK_F12);
  key_codes_[KEY_A]       = XKeysymToKeycode(display_, XK_A);
  key_codes_[KEY_Q]       = XKeysymToKeycode(display_, XK_Q);
  key_codes_[KEY_W]       = XKeysymToKeycode(display_, XK_W);
  key_codes_[KEY_Z]       = XKeysymToKeycode(display_, XK_Z);
  key_codes_[KEY_S]       = XKeysymToKeycode(display_, XK_S);
  key_codes_[KEY_X]       = XKeysymToKeycode(display_, XK_X);
  key_codes_[KEY_E]       = XKeysymToKeycode(display_, XK_E);
  key_codes_[KEY_D]       = XKeysymToKeycode(display_, XK_D);
  key_codes_[KEY_C]       = XKeysymToKeycode(display_, XK_C);
  key_codes_[KEY_R]       = XKeysymToKeycode(display_, XK_R);
  key_codes_[KEY_F]       = XKeysymToKeycode(display_, XK_F);
  key_codes_[KEY_V]       = XKeysymToKeycode(display_, XK_V);
  key_codes_[KEY_T]       = XKeysymToKeycode(display_, XK_T);
  key_codes_[KEY_G]       = XKeysymToKeycode(display_, XK_G);
  key_codes_[KEY_B]       = XKeysymToKeycode(display_, XK_B);
  key_codes_[KEY_Y]       = XKeysymToKeycode(display_, XK_Y);
  key_codes_[KEY_H]       = XKeysymToKeycode(display_, XK_H);
  key_codes_[KEY_N]       = XKeysymToKeycode(display_, XK_N);
  key_codes_[KEY_U]       = XKeysymToKeycode(display_, XK_U);
  key_codes_[KEY_J]       = XKeysymToKeycode(display_, XK_J);
  key_codes_[KEY_I]       = XKeysymToKeycode(display_, XK_I);
  key_codes_[KEY_K]       = XKeysymToKeycode(display_, XK_K);
  key_codes_[KEY_O]       = XKeysymToKeycode(display_, XK_O);
  key_codes_[KEY_L]       = XKeysymToKeycode(display_, XK_L);
  key_codes_[KEY_P]       = XKeysymToKeycode(display_, XK_P);
  key_codes_[KEY_M]       = XKeysymToKeycode(display_, XK_M);
  key_codes_[KEY_1]       = XKeysymToKeycode(display_, XK_1);
  key_codes_[KEY_2]       = XKeysymToKeycode(display_, XK_2);
  key_codes_[KEY_3]       = XKeysymToKeycode(display_, XK_3);
  key_codes_[KEY_4]       = XKeysymToKeycode(display_, XK_4);
  key_codes_[KEY_5]       = XKeysymToKeycode(display_, XK_5);
  key_codes_[KEY_6]       = XKeysymToKeycode(display_, XK_6);
  key_codes_[KEY_7]       = XKeysymToKeycode(display_, XK_7);
  key_codes_[KEY_8]       = XKeysymToKeycode(display_, XK_8);
  key_codes_[KEY_9]       = XKeysymToKeycode(display_, XK_9);
}

void ExpManager_X11::compute_colormap() {
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
